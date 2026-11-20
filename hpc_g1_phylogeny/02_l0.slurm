#!/bin/bash
#SBATCH -J g1_02_l0
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=80
#SBATCH -o 02_l0.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
SETS_DIR="$PREP_OUT/sets"
mkdir -p "$L0_OUT"/{pos,aln,dists,M,multi,main,sensR,logs} "$RUNINFO_DIR"
rm -f "$L0_OUT/$STAGE_DONE"
cd "$L0_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_G1"
print_config_banner
require_tools bcftools bedtools snp-dists python3 || exit 1
log_versions 02_l0 bcftools bedtools snp-dists python3
[[ -e "$PREP_OUT/$STAGE_DONE" ]] || { echo "[ERROR] stage 01 not complete ($PREP_OUT/$STAGE_DONE missing)" >&2; exit 1; }
load_samples || exit 1
SAMPLE_TABLE="$PREP_OUT/samples_callset.tsv"
SD_THREADS="${SNPDISTS_THREADS:-$THREADS_TOTAL}"
MAX_JOBS="$(fifo_width "${L0_MAX_JOBS:-}" "$THREADS_TOTAL" 1 200)"

# positions of each callset
bcftools query -f '%CHROM\t%POS\n' "$(base1_vcf main)" > pos/main.tsv
bcftools query -f '%CHROM\t%POS\n' "$(base1_vcf h99)" > pos/h99.tsv
set_positions() { case "$1" in main|sensR_main) echo pos/main.tsv ;; h99|sensR_h99) echo pos/h99.tsv ;; esac; }

# alignments + snp-dists, every set
for set in "${PSEUDO_SETS[@]}"; do
  python3 "$PIPELINE_DIR/l0_distances.py" extract --fasta-dir "$SETS_DIR/$set/fasta" --samples "$SETS_DIR/$set/samples.txt" \
    --positions "$(set_positions "$set")" --out "aln/${set}.fa" || exit 1
  snp-dists -j "$SD_THREADS" -m "aln/${set}.fa" > "dists/${set}.molten.tsv" 2> "logs/snp-dists.${set}.log" || { echo "[ERROR] snp-dists failed: $set" >&2; exit 1; }
  echo "[$(date '+%F %T')] >>> ${set}: snp-dists on $(awk 'END{print NR / 2}' "aln/${set}.fa") sequences"
done

# M_i = A' - U_i and multiinter, sets of the two analyses
m_job() {   # <set> <sample>
  bedtools subtract -a "$SETS_DIR/$1/aprime.bed" -b "$SETS_DIR/$1/usable/$2.usable.bed" > "M/$1/$2.bed" 2>>"logs/M.$1.log" \
    || { echo "[ERROR] bedtools subtract failed: $1/$2"; touch "M/$1/$2.failed"; return 1; }
}
for set in main h99 sensR_main sensR_h99; do
  mkdir -p "M/$set"; rm -f "M/$set"/*
  open_fifo "$MAX_JOBS"
  while read -r s; do run_parallel m_job "$set" "$s"; done < "$SETS_DIR/$set/samples.txt"
  close_fifo
  files=(); names=()
  while read -r s; do
    [[ ! -e "M/$set/$s.failed" && -e "M/$set/$s.bed" ]] || { echo "[ERROR] M missing for $set/$s" >&2; exit 1; }
    if [[ -s "M/$set/$s.bed" ]]; then files+=( "M/$set/$s.bed" ); names+=( "$s" ); fi
  done < "$SETS_DIR/$set/samples.txt"
  if (( ${#files[@]} >= 2 )); then
    bedtools multiinter -i "${files[@]}" -names "${names[@]}" > "multi/${set}.tsv"
  elif (( ${#files[@]} == 1 )); then
    awk -v n="${names[0]}" 'BEGIN{OFS="\t"} {print $1, $2, $3, 1, n}' "${files[0]}" > "multi/${set}.tsv"
  else
    : > "multi/${set}.tsv"
  fi
  echo "[$(date '+%F %T')] >>> ${set}: ${#files[@]} samples with unusable bases inside A'; multiinter $(awk 'END{print NR}' "multi/${set}.tsv") intervals"
done

# per analysis: pairwise table, rulers, SNP address
for an in main sensR; do
  if [[ "$an" == main ]]; then sm=main; sh=h99; else sm=sensR_main; sh=sensR_h99; fi
  od="$L0_OUT/$an"
  python3 "$PIPELINE_DIR/l0_distances.py" assemble --sample-table "$SAMPLE_TABLE" --metadata "$METADATA_TSV" \
    --set-samples "main=$SETS_DIR/$sm/samples.txt" "h99=$SETS_DIR/$sh/samples.txt" \
    --dists "main=dists/${sm}.molten.tsv" "h99=dists/${sh}.molten.tsv" \
    --multi "main=multi/${sm}.tsv" "h99=multi/${sh}.tsv" \
    --aprime "main=$SETS_DIR/$sm/aprime.bed" "h99=$SETS_DIR/$sh/aprime.bed" \
    --stats "main=$SETS_DIR/$sm/stats.tsv" "h99=$SETS_DIR/$sh/stats.tsv" \
    --gate-frac "$PAIR_GATE_FRAC" --out-dir "$od" || exit 1
  python3 "$PIPELINE_DIR/rulers_address.py" --pairs "$od/pairwise_L0.tsv" --sample-table "$SAMPLE_TABLE" --metadata "$METADATA_TSV" \
    --qc-gate "$QC_GATE_TSV" --usable-stats "$SETS_DIR/$sm/stats.tsv" --ladder "$LADDER" \
    --t-tech-percentile "$T_TECH_PERCENTILE" --t-event-round "$T_EVENT_ROUND" --low-usable-frac "$LOW_USABLE_FRAC" \
    --label "$an" --out-dir "$od" || exit 1
done

echo "=============================================================="
echo "[$(date '+%F %T')] >>> L0 done: $L0_OUT/{main,sensR}/ (pairwise_L0.tsv, ruler_*.tsv, anchors.tsv, snp_address*.tsv)"
echo "=============================================================="
mark_stage_done "$L0_OUT"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
