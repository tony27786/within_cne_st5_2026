#!/bin/bash
#SBATCH -J g0_X1_case
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH -o X1_case_h99_check.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$X1_OUT"/{genomicsdb,by_contig,logs,tmp,gt_source,report} "$RUNINFO_DIR"
cd "$X1_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner

require_tools gatk bcftools samtools python3 || exit 1
log_versions X1_case_h99_check gatk bcftools samtools python3
for f in "$REF_FA" "$REF_DICT" "$MASK_BED" "$ANALYZABLE_BED" "$GENOME_TXT" "$CONTIG_LIST" "$JOINT_OUT/main/st5_joint.cna3.raw.vcf.gz" "$PIPELINE_DIR/gt_source.py" "$PIPELINE_DIR/gt_mask.py" "$PIPELINE_DIR/case_h99_relationship.py"; do
  [[ -s "$f" ]] || { echo "[ERROR] missing $f" >&2; exit 1; }
done
load_samples || exit 1

# Samples: H99 controls + case-report isolates
X1_SAMPLES=()
for i in "${!S_NAME[@]}"; do
  case "${S_ROLE[$i]}" in technical_replicate|case_report_ST2_to_ST5) X1_SAMPLES+=( "${S_NAME[$i]}" ) ;; esac
done
MAP_TSV="$X1_OUT/sample_map.tsv"; : > "$MAP_TSV"
for s in "${X1_SAMPLES[@]}"; do
  for f in "$GVCF_OUT/${s}.g.vcf.gz" "$CALL_OUT/callable/${s}.cna3.bed.gz" "$MAP_OUT/${s}.md.bam"; do
    [[ -s "$f" ]] || { echo "[ERROR] missing $f" >&2; exit 1; }
  done
  printf '%s\t%s\n' "$s" "$GVCF_OUT/${s}.g.vcf.gz" >> "$MAP_TSV"
done
echo "[$(date '+%F %T')] >>> X1 samples: $(printf '%s ' "${X1_SAMPLES[@]}")"

FIFO_W() { fifo_width "" "$THREADS_TOTAL" 1 "$1"; }
open_fifo() { SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"; local i; for ((i=0;i<$1;i++)); do echo >&3; done; }
close_fifo() { wait || true; exec 3>&- 3<&- || true; }
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

# 1) Joint genotyping of the 7 samples, per contig
mapfile -t CONTIGS < "$CONTIG_LIST"
joint_contig() {
  local c="$1" ws="$X1_OUT/genomicsdb/$1" out="$X1_OUT/by_contig/$1.vcf.gz" log="$X1_OUT/logs/$1.log" tmp="$X1_OUT/tmp/$1"
  mkdir -p "$tmp"; : > "$log"
  local jopts; jopts="$(gatk_java_opts "$JOINT_JAVA_XMX" "$tmp")"
  rm -f "$out" "${out}.tbi" "$X1_OUT/by_contig/$1.failed"
  if ! gatk --java-options "$jopts" GenomicsDBImport --genomicsdb-workspace-path "$ws" --overwrite-existing-genomicsdb-workspace true --sample-name-map "$MAP_TSV" -L "$INTERVAL_DIR/${c}.intervals" --batch-size "$GDB_BATCH_SIZE" --reader-threads "$JOINT_READER_THREADS" --genomicsdb-shared-posixfs-optimizations true --tmp-dir "$tmp" >>"$log" 2>&1 \
     || ! gatk --java-options "$jopts" GenotypeGVCFs -R "$REF_FA" -V "gendb://${ws}" -L "$INTERVAL_DIR/${c}.intervals" -O "$out" --genomicsdb-shared-posixfs-optimizations true --tmp-dir "$tmp" >>"$log" 2>&1; then
    echo "[ERROR] joint genotyping failed: $c (see $log)"; touch "$X1_OUT/by_contig/$1.failed"; return 1
  fi
  rm -rf "$tmp"
  echo "[$(date '+%F %T')] >>> joint done: $c"
}
echo "[$(date '+%F %T')] >>> 1) joint genotyping, ${#CONTIGS[@]} contigs, $(FIFO_W "${#CONTIGS[@]}") in parallel"
open_fifo "$(FIFO_W "${#CONTIGS[@]}")"
for c in "${CONTIGS[@]}"; do run_parallel joint_contig "$c"; done
close_fifo
PARTS=()
for c in "${CONTIGS[@]}"; do
  [[ -s "$X1_OUT/by_contig/${c}.vcf.gz" && ! -e "$X1_OUT/by_contig/${c}.failed" ]] || { echo "[ERROR] contig $c failed; stopping" >&2; exit 1; }
  PARTS+=( "$X1_OUT/by_contig/${c}.vcf.gz" )
done
RAW="$X1_OUT/x1.cna3.raw.vcf.gz"
bcftools concat --threads "$BCF_THREADS" -Oz -o "$RAW" "${PARTS[@]}"
bcftools index --threads "$BCF_THREADS" -f -t "$RAW"

# 2) Same site rules as 06
JOPTS="$(gatk_java_opts "$FILTER_JAVA_XMX" "$X1_OUT/tmp")"
mapfile -t VF_ARGS < <(vf_filter_args)
SNP_RAW="$X1_OUT/x1.cna3.snp.raw.vcf.gz"; SITEFILT="$X1_OUT/x1.cna3.snp.sitefilt.vcf.gz"; PASSMASK="$X1_OUT/x1.cna3.snp.pass_unmasked.vcf.gz"; X1_VCF="$X1_OUT/x1.cna3.snp.filtered.vcf.gz"
gatk --java-options "$JOPTS" SelectVariants -R "$REF_FA" -V "$RAW" --select-type-to-include SNP -O "$SNP_RAW" >>"$X1_OUT/logs/filter.log" 2>&1
gatk --java-options "$JOPTS" VariantFiltration -R "$REF_FA" -V "$SNP_RAW" -O "$SITEFILT" "${VF_ARGS[@]}" >>"$X1_OUT/logs/filter.log" 2>&1
bcftools view -f PASS -T "^${MASK_BED}" --threads "$BCF_THREADS" -Oz -o "$PASSMASK" "$SITEFILT"
bcftools index --threads "$BCF_THREADS" -f -t "$PASSMASK"
bcftools query -f '%CHROM\t%POS\n' "$PASSMASK" > "$X1_OUT/sites.tsv"
echo "[$(date '+%F %T')] >>> 2) site filters: $(bcftools view -H "$SNP_RAW" | awk 'END{print NR}') SNPs -> $(awk 'END{print NR}' "$X1_OUT/sites.tsv") PASS outside the static mask"

# 3) Genotype source, genotype-level N
src_job() {
  python3 "$PIPELINE_DIR/gt_source.py" --gvcf "$GVCF_OUT/$1.g.vcf.gz" --sites "$X1_OUT/sites.tsv" --callable "$CALL_OUT/callable/$1.cna3.bed.gz" --genome "$GENOME_TXT" --out "$X1_OUT/gt_source/$1.tsv.gz" > "$X1_OUT/logs/gt_source.$1.log" 2>&1 || { echo "[ERROR] gt_source.py failed: $1"; touch "$X1_OUT/gt_source/$1.failed"; return 1; }
}
rm -f "$X1_OUT"/gt_source/*
open_fifo "$(FIFO_W "${#X1_SAMPLES[@]}")"
for s in "${X1_SAMPLES[@]}"; do run_parallel src_job "$s"; done
close_fifo
for s in "${X1_SAMPLES[@]}"; do
  [[ -s "$X1_OUT/gt_source/${s}.tsv.gz" ]] || { echo "[ERROR] genotype source missing for $s" >&2; exit 1; }
done
bcftools view -Ov "$PASSMASK" \
  | python3 "$PIPELINE_DIR/gt_mask.py" --gq-min "$GT_GQ_MIN" --dp-min "$GT_DP_MIN" --ad-frac-min "$GT_AD_FRAC_MIN" --src-dir "$X1_OUT/gt_source" --kept-out "$X1_OUT/refblock_kept.tsv.gz" --stats "$X1_OUT/gt_mask_stats.tsv" --site-stats "$X1_OUT/gt_mask_sites.tsv" \
  | bcftools view --threads "$BCF_THREADS" -Oz -o "$X1_VCF"
bcftools index --threads "$BCF_THREADS" -f -t "$X1_VCF"

# 4) Relationship report
python3 "$PIPELINE_DIR/case_h99_relationship.py" --vcf "$X1_VCF" --callable-dir "$CALL_OUT/callable" --analyzable "$ANALYZABLE_BED" --mask "$MASK_BED" --main-raw "$JOINT_OUT/main/st5_joint.cna3.raw.vcf.gz" --samples-tsv "$SAMPLES_TSV" --raw-root "$RAW_ROOT" --out-dir "$X1_OUT/report"

echo "=============================================================="
echo "[$(date '+%F %T')] >>> X1 report: $X1_OUT/report/ (summary.txt, pairwise_distance.tsv, h99_allele_carriers.tsv, site_patterns.tsv, sequencing_metadata.tsv)"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
