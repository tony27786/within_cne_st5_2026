#!/bin/bash
#SBATCH -J g1_01_prep
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=80
#SBATCH -o 01_prep.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
R_DIR="$PREP_OUT/r"; X_DIR="$PREP_OUT/xmask"; SETS_DIR="$PREP_OUT/sets"; FAIL_DIR="$PREP_OUT/.failed"; LOG_DIR="$PREP_OUT/logs"
mkdir -p "$R_DIR" "$X_DIR" "$SETS_DIR" "$FAIL_DIR" "$LOG_DIR" "$PREP_OUT/tmp" "$RUNINFO_DIR"
rm -f "$FAIL_DIR"/* "$PREP_OUT/$STAGE_DONE" 2>/dev/null || true

cd "$PREP_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_G1"
print_config_banner
require_tools bcftools bedtools samtools python3 || exit 1
log_versions 01_prep bcftools bedtools samtools python3

# 0) inputs
load_samples || exit 1
write_sample_table "$PREP_OUT/samples_callset.tsv"
need=( "$REF_FA" "${REF_FA}.fai" "$GENOME_TXT" "$GENOME_SORTED_TXT" "$ANALYZABLE_STATIC_BED" "$ANALYZABLE_UNION_BED" "$QC_GATE_TSV" "$METADATA_TSV" )
for cs in main h99; do need+=( "$(base1_vcf "$cs")" "$(base1_vcf "$cs").tbi" "$(sitefilt_vcf "$cs")" "$(indel_vcf "$cs")" ); done
for s in "${S_NAME[@]}"; do need+=( "$(base2_bed "$s")" ); done
require_files "${need[@]}" || exit 1
echo "[$(date '+%F %T')] >>> 0) inputs OK: G0 run ${G0_RUN_TAG}"

# 1) R list
EXPR="FMT/ADV=\"M\" & sSUM(FMT/AD)>=${R_MIN_AD_SUM} & sMAX(FMT/AD)/sSUM(FMT/AD)>=${R_IAF_LO} & sMAX(FMT/AD)/sSUM(FMT/AD)<=${R_IAF_HI}"
echo "[$(date '+%F %T')] >>> 1) R list: per-genotype expression: $EXPR"
bcftools query -i "$EXPR" -f '[%CHROM\t%POS\t%SAMPLE\n]' "$(base1_vcf main)" > "$R_DIR/new_iaf_genotypes.tsv"
python3 "$PIPELINE_DIR/r_list.py" --new-geno "$R_DIR/new_iaf_genotypes.tsv" --min-strains "$R_MIN_STRAINS" --out-dir "$R_DIR" || exit 1
md5_of "$R_DIR/R.bed" > "$R_DIR/R.bed.md5"

# 2) site-exclusion regions
echo "[$(date '+%F %T')] >>> 2) site-exclusion regions X_c"
for cs in main h99; do
  python3 "$PIPELINE_DIR/xmask.py" --name "$cs" --sitefilt "$(sitefilt_vcf "$cs")" --nonsnp "$(indel_vcf "$cs")" \
    --out-bed "$X_DIR/${cs}.bed" --out-summary "$X_DIR/${cs}.summary.tsv" || exit 1
done

# 3) common masks and A' per set
sort_bed() { LC_ALL=C sort -k1,1 -k2,2n "$@"; }
for set in "${PSEUDO_SETS[@]}"; do
  sd="$SETS_DIR/$set"; mkdir -p "$sd"/{masks,fasta,usable,stats,dump}
  set_samples "$set" > "$sd/samples.txt"
  region="$(set_region "$set")"; xname="$(set_xname "$set")"
  {
    sort_bed "$region" | bedtools complement -i - -g "$GENOME_SORTED_TXT"
    cat "$X_DIR/${xname}.bed"
    if set_uses_r "$set"; then cat "$R_DIR/R.bed"; fi
  } | sort_bed | bedtools merge -i - > "$sd/common_mask.bed"
  bedtools complement -i "$sd/common_mask.bed" -g "$GENOME_SORTED_TXT" > "$sd/aprime.bed"
  echo "[$(date '+%F %T')] >>> 3) set ${set}: $(awk 'END{print NR}' "$sd/samples.txt") samples; region $(basename "$region"); X ${xname}; R $(set_uses_r "$set" && echo yes || echo no); A' $(bed_bp "$sd/aprime.bed") bp"
done

# 4) pseudo-genomes
JOBS=()
for set in "${PSEUDO_SETS[@]}"; do while read -r s; do JOBS+=( "${set}:${s}" ); done < "$SETS_DIR/$set/samples.txt"; done
TOTAL_JOBS="${#JOBS[@]}"
MAX_JOBS="$(fifo_width "${PREP_MAX_JOBS:-}" "$THREADS_TOTAL" 1 "$TOTAL_JOBS")"
echo "[$(date '+%F %T')] >>> 4) ${TOTAL_JOBS} pseudo-genomes, ${MAX_JOBS} in parallel"
progress_init "$TOTAL_JOBS"

pseudo_fail() {   # <set> <sample> <reason>
  echo "[ERROR] $3: $1/$2" | tee -a "$LOG_DIR/${1}__${2}.log"
  echo "$3" > "$FAIL_DIR/${1}__${2}"
  rm -f "$SETS_DIR/$1/fasta/$2.fa" "$SETS_DIR/$1/fasta/$2.fa.fai" "$SETS_DIR/$1/usable/$2.usable.bed" "$SETS_DIR/$1/stats/$2.tsv"
  update_progress "$1/$2"
}
pseudo_job() {
  local set="${1%%:*}" s="${1#*:}"
  local sd="$SETS_DIR/$set" log="$LOG_DIR/${set}__${s}.log"; : > "$log"
  local mask="$sd/masks/$s.mask.bed" fa="$sd/fasta/$s.fa" vcf
  vcf="$(set_vcf "$set")"

  # To resume a partial run, uncomment to skip completed samples:
  # if [[ -s "$sd/usable/$s.usable.bed" && -s "$sd/stats/$s.tsv" ]]; then update_progress "$set/$s"; return 0; fi

  if ! { gzip -dc "$(base2_bed "$s")" | bedtools complement -i - -g "$GENOME_SORTED_TXT"; cat "$sd/common_mask.bed"; } \
       | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > "$mask" 2>>"$log"; then
    pseudo_fail "$set" "$s" "mask construction failed"; return 1
  fi
  if ! bcftools consensus -f "$REF_FA" -s "$s" -H 1 -M N -m "$mask" -o "$fa" "$vcf" >>"$log" 2>&1 || ! samtools faidx "$fa" 2>>"$log"; then
    pseudo_fail "$set" "$s" "bcftools consensus / faidx failed"; return 1
  fi
  if ! bcftools query -s "$s" -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%SRC]\n' "$vcf" > "$sd/dump/$s.tsv" 2>>"$log"; then
    pseudo_fail "$set" "$s" "genotype dump failed"; return 1
  fi
  if ! python3 "$PIPELINE_DIR/pseudo_check.py" --sample "$s" --pseudo "$fa" --dump "$sd/dump/$s.tsv" \
         --out-usable "$sd/usable/$s.usable.bed" --out-stats "$sd/stats/$s.tsv" >>"$log" 2>&1; then
    pseudo_fail "$set" "$s" "pseudo_check.py failed (see $log)"; return 1
  fi
  rm -f "$sd/dump/$s.tsv"
  update_progress "$set/$s"
}
open_fifo "$MAX_JOBS"
for j in "${JOBS[@]}"; do run_parallel pseudo_job "$j"; done
close_fifo
progress_done

# collect
n_ok=0
for j in "${JOBS[@]}"; do
  set="${j%%:*}"; s="${j#*:}"
  if [[ -s "$SETS_DIR/$set/usable/$s.usable.bed" && -s "$SETS_DIR/$set/stats/$s.tsv" ]]; then n_ok=$((n_ok + 1))
  elif [[ ! -e "$FAIL_DIR/${set}__${s}" ]]; then echo "no output and no error recorded (job killed?)" > "$FAIL_DIR/${set}__${s}"; fi
done
write_failed_list "$FAIL_DIR" "$PREP_OUT/FAILED_SAMPLES.txt"
for set in "${PSEUDO_SETS[@]}"; do
  sd="$SETS_DIR/$set"
  { printf 'sample\tusable_bp\tn_base1_sites\tn_gt_missing\n'; while read -r s; do tail -n +2 "$sd/stats/$s.tsv" 2>/dev/null || true; done < "$sd/samples.txt"; } > "$sd/stats.tsv"
done
echo "[$(date '+%F %T')] >>> pseudo-genomes complete: ${n_ok} / ${TOTAL_JOBS}"
FAILED=( "$FAIL_DIR"/* )
if (( ${#FAILED[@]} > 0 )); then
  echo "[ERROR] ${#FAILED[@]} pseudo-genome(s) failed (all are required downstream):" >&2
  for f in "${FAILED[@]}"; do printf '   %s  (%s)\n' "$(basename "$f")" "$(head -n 1 "$f")" >&2; done
  exit 1
fi
for set in "${PSEUDO_SETS[@]}"; do
  awk -F'\t' -v s="$set" 'NR > 1 {n++; u += $2; if (min == "" || $2 < min) min = $2; if ($2 > max) max = $2} END{printf "    %-11s %3d samples, usable bp min %d / mean %.0f / max %d\n", s, n, min, u / n, max}' "$SETS_DIR/$set/stats.tsv"
done

mark_stage_done "$PREP_OUT"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
