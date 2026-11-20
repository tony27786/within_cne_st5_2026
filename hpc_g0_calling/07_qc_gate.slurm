#!/bin/bash
#SBATCH -J g0_07_qc
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=80
#SBATCH -o 07_qc_gate.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
WORK_DIR="$QC_OUT/work"
mkdir -p "$QC_OUT" "$WORK_DIR" "$RUNINFO_DIR"

cd "$QC_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner

require_tools bedtools python3 || exit 1
log_versions 07_qc_gate bedtools python3
for f in "$MASK_SORTED_BED" "$GENOME_SORTED_TXT" "$NUCLEAR_BED" "$ANALYZABLE_BED" "$CONTIG_LIST" "$PIPELINE_DIR/qc_gate_table.py"; do
  [[ -s "$f" ]] || { echo "[ERROR] missing $f" >&2; exit 1; }
done
load_samples || exit 1
# The gate needs all three base-1 VCFs (afterok guarantees this when submitted by submit_all.sh; this guards manual runs).
for cs in "${CALLSETS[@]}"; do
  b1="$FILT_OUT/$cs/$(callset_prefix "$cs").cna3.snp.filtered.vcf.gz"
  [[ -s "$b1" ]] || { echo "[ERROR] base 1 missing for callset ${cs}: $b1 (run 06_filter first)" >&2; exit 1; }
done
status=0

# Sample table for the Python step: callset membership comes from load_samples (single source of truth for the role -> callset rule).
{
  printf 'sample\tgroup\trole\tcallset\n'
  for i in "${!S_NAME[@]}"; do
    printf '%s\t%s\t%s\t%s\n' "${S_NAME[$i]}" "${S_GROUP[$i]}" "${S_ROLE[$i]}" "${S_CALLSET[$i]}"
  done
} > "$WORK_DIR/samples_callset.tsv"

bed_bp()     { awk '{s += $3 - $2} END{print s + 0}' "$1"; }
bed_bp_nuc() { awk -v m="$MITO_CONTIG" '$1 != m {s += $3 - $2} END{print s + 0}' "$1"; }

# 1) Empirical inaccessible region + union with the static mask
REP_BEDS=(); MISSING_REPS=()
for i in "${!S_NAME[@]}"; do
  if [[ "${S_ROLE[$i]}" == "$EMPIRICAL_ROLE" ]]; then
    b="$CALL_OUT/callable/${S_NAME[$i]}.cna3.bed.gz"
    if [[ -s "$b" ]]; then REP_BEDS+=( "$b" ); else MISSING_REPS+=( "${S_NAME[$i]}" ); fi
  fi
done
EMP_BED="$QC_OUT/empirical_inaccessible.bed"
UNION_BED="$QC_OUT/mask_union.bed"
ANALYZABLE_UNION_BED="$QC_OUT/analyzable_union.bed"
SUMMARY_TSV="$QC_OUT/mask_union_summary.tsv"
rm -f "$EMP_BED" "$UNION_BED" "$ANALYZABLE_UNION_BED" "$SUMMARY_TSV"

if (( ${#MISSING_REPS[@]} > 0 )); then
  echo "[ERROR] empirical mask needs a base-2 mask for every representative; found ${#REP_BEDS[@]}, missing ${#MISSING_REPS[@]}" >&2
  if (( ${#MISSING_REPS[@]} > 0 )); then printf '   %s\n' "${MISSING_REPS[@]}" >&2; fi
  echo "        -> skipping the empirical mask; QC table is still written; stage will exit 1." >&2
  status=1
else
  n_rep=${#REP_BEDS[@]}
  k_min=$(( (EMPIRICAL_MIN_CALLABLE_PCT * n_rep + 99) / 100 ))    # 144 of 160
  echo "[$(date '+%F %T')] >>> empirical mask: accessible = callable in >= ${k_min} of ${n_rep} representatives (< ${EMPIRICAL_MIN_CALLABLE_PCT}% -> inaccessible)"
  if ! bedtools multiinter -i "${REP_BEDS[@]}" \
       | awk -v k="$k_min" 'BEGIN{OFS="\t"} $4 + 0 >= k {print $1, $2, $3}' \
       | bedtools merge -i - > "$WORK_DIR/accessible_ge${k_min}of${n_rep}.bed"; then
    echo "[ERROR] bedtools multiinter/merge failed" >&2; exit 1
  fi
  bedtools complement -i "$WORK_DIR/accessible_ge${k_min}of${n_rep}.bed" -g "$GENOME_SORTED_TXT" > "$EMP_BED"
  cat "$MASK_SORTED_BED" "$EMP_BED" | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > "$UNION_BED"
  bedtools subtract -a "$NUCLEAR_BED" -b "$UNION_BED" | LC_ALL=C sort -k1,1 -k2,2n > "$ANALYZABLE_UNION_BED"

  static_bp=$(bed_bp "$MASK_SORTED_BED");  static_nuc=$(bed_bp_nuc "$MASK_SORTED_BED")
  emp_bp=$(bed_bp "$EMP_BED");             emp_nuc=$(bed_bp_nuc "$EMP_BED")
  union_bp=$(bed_bp "$UNION_BED");         union_nuc=$(bed_bp_nuc "$UNION_BED")
  overlap_bp=$(( static_bp + emp_bp - union_bp ))
  emp_only_bp=$(( union_bp - static_bp ))
  analyzable_union_bp=$(bed_bp "$ANALYZABLE_UNION_BED")
  nuclear_bp=$(bed_bp "$NUCLEAR_BED")
  {
    printf 'metric\tvalue\n'
    printf 'representatives_used\t%s\n' "$n_rep"
    printf 'accessible_min_callable_samples\t%s\n' "$k_min"
    printf 'static_mask_intervals\t%s\n' "$(awk 'END{print NR}' "$MASK_SORTED_BED")"
    printf 'static_mask_bp\t%s\n' "$static_bp"
    printf 'static_mask_nuclear_bp\t%s\n' "$static_nuc"
    printf 'empirical_inaccessible_intervals\t%s\n' "$(awk 'END{print NR}' "$EMP_BED")"
    printf 'empirical_inaccessible_bp\t%s\n' "$emp_bp"
    printf 'empirical_inaccessible_nuclear_bp\t%s\n' "$emp_nuc"
    printf 'overlap_static_empirical_bp\t%s\n' "$overlap_bp"
    printf 'empirical_only_bp\t%s\n' "$emp_only_bp"
    printf 'union_intervals\t%s\n' "$(awk 'END{print NR}' "$UNION_BED")"
    printf 'union_bp\t%s\n' "$union_bp"
    printf 'union_nuclear_bp\t%s\n' "$union_nuc"
    awk -v u="$union_nuc" -v n="$nuclear_bp" 'BEGIN{printf "union_nuclear_pct_of_nuclear_genome\t%.4f\n", 100 * u / n}'
    awk -v s="$static_nuc" -v n="$nuclear_bp" 'BEGIN{printf "static_nuclear_pct_of_nuclear_genome\t%.4f\n", 100 * s / n}'
    printf 'analyzable_union_bp\t%s\n' "$analyzable_union_bp"
  } > "$SUMMARY_TSV"
  echo "[$(date '+%F %T')] >>> static U empirical:"
  sed 's/^/    /' "$SUMMARY_TSV"
fi

# 2) Per-sample callable bp inside the union-masked region (informational; FIFO parallel)
UNION_BP_DIR="$WORK_DIR/callable_union_bp"
rm -rf "$UNION_BP_DIR"; mkdir -p "$UNION_BP_DIR"
union_bp_job() {   # <sample>
  local s="$1" b="$CALL_OUT/callable/${1}.cna3.bed.gz" v="NA"
  if [[ -s "$b" && -s "$ANALYZABLE_UNION_BED" ]]; then
    v="$(bedtools intersect -a "$b" -b "$ANALYZABLE_UNION_BED" | awk '{s += $3 - $2} END{print s + 0}')" || v="NA"
  fi
  printf '%s\t%s\n' "$s" "$v" > "$UNION_BP_DIR/${s}.tsv"
}
MAX_JOBS="$(fifo_width "${QC_MAX_JOBS:-}" "$THREADS_TOTAL" 1 "${#S_NAME[@]}")"   # bedtools intersect is single-threaded
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }
for s in "${S_NAME[@]}"; do
  run_parallel union_bp_job "$s"
done
wait || true
exec 3>&- 3<&- || true
{
  printf 'sample\tcallable_bp_analyzable_union\n'
  for s in "${S_NAME[@]}"; do
    if [[ -s "$UNION_BP_DIR/${s}.tsv" ]]; then cat "$UNION_BP_DIR/${s}.tsv"; else printf '%s\tNA\n' "$s"; fi
  done
} > "$WORK_DIR/callable_union_bp.tsv"
echo "[$(date '+%F %T')] >>> callable bp inside the union-masked region: ${#S_NAME[@]} samples ($MAX_JOBS in parallel)"

# 3) QC gate table
BASE1_ARGS=(); GTMASK_ARGS=()
for cs in "${CALLSETS[@]}"; do
  prefix="$(callset_prefix "$cs")"
  BASE1_ARGS+=( --base1 "${cs}=$FILT_OUT/$cs/${prefix}.cna3.snp.filtered.vcf.gz" )
  GTMASK_ARGS+=( --gtmask "${cs}=$FILT_OUT/$cs/gt_mask_stats.tsv" )
done
MASK_SUMMARY_ARGS=()
if [[ -s "$SUMMARY_TSV" ]]; then MASK_SUMMARY_ARGS=( --mask-summary "$SUMMARY_TSV" ); fi

if ! python3 "$PIPELINE_DIR/qc_gate_table.py" \
      --samples "$WORK_DIR/samples_callset.tsv" \
      --analyzable-bp "$(bed_bp "$ANALYZABLE_BED")" \
      --callable-dir "$CALL_OUT/summary" \
      --union-callable "$WORK_DIR/callable_union_bp.tsv" \
      --flagstat-dir "$MAP_OUT/qc_metrics_map/flagstat" \
      --contigs "$CONTIG_LIST" --mito "$MITO_CONTIG" \
      "${BASE1_ARGS[@]}" "${GTMASK_ARGS[@]}" ${MASK_SUMMARY_ARGS[@]+"${MASK_SUMMARY_ARGS[@]}"} \
      --ref-role "$QC_REF_ROLE" \
      --callable-ratio-min "$QC_CALLABLE_RATIO_MIN" \
      --depth-median-min "$QC_DEPTH_MEDIAN_MIN" \
      --iaf-lo "$QC_IAF_LO" --iaf-hi "$QC_IAF_HI" \
      --iaf-min-ad "$QC_IAF_MIN_AD_SUM" \
      --refdist-z-max "$QC_REFDIST_ROBUST_Z_MAX" \
      --out-table "$QC_OUT/qc_gate.tsv" \
      --out-summary "$QC_OUT/qc_gate_summary.txt" \
      --out-depth "$QC_OUT/depth_by_chrom.tsv" \
      --out-iaf-sites "$QC_OUT/intermediate_af_sites_h99_case.tsv"; then
  echo "[ERROR] qc_gate_table.py failed" >&2
  exit 1
fi

echo "=============================================================="
echo "[$(date '+%F %T')] >>> G0 QC gate written:"
echo "    $QC_OUT/qc_gate.tsv          (one row per sample)"
echo "    $QC_OUT/qc_gate_summary.txt"
echo "    $QC_OUT/depth_by_chrom.tsv"
echo "    $QC_OUT/intermediate_af_sites_h99_case.tsv"
if [[ -s "$SUMMARY_TSV" ]]; then
  echo "    $UNION_BED  (static U empirical; sizes in $SUMMARY_TSV)"
fi
echo "=============================================================="
if (( status != 0 )); then
  echo "[ERROR] empirical inaccessible mask was NOT produced (missing representative masks)." >&2
  exit 1
fi
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
