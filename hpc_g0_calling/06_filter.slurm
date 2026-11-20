#!/bin/bash
#SBATCH -J g0_06_filter
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=80
#SBATCH -o 06_filter.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
FAIL_DIR="$FILT_OUT/.failed"
mkdir -p "$FILT_OUT" "$FAIL_DIR" "$RUNINFO_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$FILT_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner

require_tools gatk bcftools samtools python3 || exit 1
log_versions 06_filter gatk bcftools samtools python3
for f in "$REF_FA" "$REF_DICT" "$MASK_BED" "$GENOME_TXT" "$PIPELINE_DIR/gt_mask.py" "$PIPELINE_DIR/gt_source.py"; do
  [[ -s "$f" ]] || { echo "[ERROR] missing $f" >&2; exit 1; }
done
load_samples || exit 1

mapfile -t FILTER_NAMES < <(filter_names)
mapfile -t VF_ARGS < <(vf_filter_args)
echo "[$(date '+%F %T')] >>> VariantFiltration arguments:"
printf '    %s\n' "${VF_ARGS[@]}"
echo "[$(date '+%F %T')] >>> genotype N if: GQ<${GT_GQ_MIN} | DP<${GT_DP_MIN} | AD[GT]/sum(AD)<${GT_AD_FRAC_MIN} | GT='*'; rule A for reference-block REF calls"

# FIFO semaphore (re-created per phase)
open_fifo() {   # <width>
  SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
  local i; for ((i=0;i<$1;i++)); do echo >&3; done
}
close_fifo() { wait || true; exec 3>&- 3<&- || true; }
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

# Progress bar (reset per phase)
PROGRESS_WIDTH=45
COUNTER_FILE="$(mktemp)"; LOCK_FILE="$(mktemp)"
progress_reset() { TOTAL_JOBS="$1"; echo 0 > "$COUNTER_FILE"; }
update_progress() {
  local sample="$1"
  {
    flock -x 9
    local done_jobs; done_jobs="$(<"$COUNTER_FILE")"; done_jobs=$((done_jobs + 1))
    echo "$done_jobs" > "$COUNTER_FILE"
    local percent=$(( done_jobs * 100 / TOTAL_JOBS ))
    local bar_len=$(( done_jobs * PROGRESS_WIDTH / TOTAL_JOBS )) bar="" i
    for ((i=0; i<bar_len; i++)); do bar+="#"; done
    for ((i=bar_len; i<PROGRESS_WIDTH; i++)); do bar+="."; done
    printf "[PROGRESS] [%s] %3d%% (%d/%d) latest: %s\n" \
      "$bar" "$percent" "$done_jobs" "$TOTAL_JOBS" "$sample"
  } 9>>"$LOCK_FILE"
}

count_records() { bcftools view -H "$1" | awk 'END{print NR}'; }

filt_fail() {   # <callset> <reason>
  local cs="$1" reason="$2" prefix
  prefix="$(callset_prefix "$cs")"
  echo "[ERROR] ${reason}: callset ${cs}" | tee -a "$FILT_OUT/$cs/logs/${prefix}.filter.log"
  echo "$reason" > "$FAIL_DIR/$cs"
  rm -f "$FILT_OUT/$cs/${prefix}.cna3.snp.filtered.tmp.vcf.gz" "$FILT_OUT/$cs/${prefix}.cna3.snp.filtered.tmp.vcf.gz.tbi"
  rm -f "$FILT_OUT/$cs/${prefix}.cna3.snp.filtered.vcf.gz" "$FILT_OUT/$cs/${prefix}.cna3.snp.filtered.vcf.gz.tbi"
  update_progress "$cs"
}

# Phase 1: site filters + static mask, one callset
prep_callset() {
  local cs="$1" prefix od log jopts
  prefix="$(callset_prefix "$cs")"
  od="$FILT_OUT/$cs"; mkdir -p "$od/logs" "$od/tmp"
  log="$od/logs/${prefix}.filter.log"; : > "$log"
  jopts="$(gatk_java_opts "$FILTER_JAVA_XMX" "$od/tmp")"
  local in="$JOINT_OUT/$cs/${prefix}.cna3.raw.vcf.gz"
  local snp_raw="$od/${prefix}.cna3.snp.raw.vcf.gz"
  local indel_raw="$od/${prefix}.cna3.indel.raw.vcf.gz"
  local sitefilt="$od/${prefix}.cna3.snp.sitefilt.vcf.gz"
  local passmask="$od/${prefix}.cna3.snp.pass_unmasked.vcf.gz"
  local n_raw n_snp n_nonsnp

  echo "[$(date '+%F %T')] >>> START phase 1: ${cs} (${in})" | tee -a "$log"
  rm -f "$od/.phase1.ok"
  if [[ ! -s "$in" || ! -s "${in}.tbi" ]]; then
    filt_fail "$cs" "raw joint VCF missing (stage 05): $in"; return 1
  fi

  # 1) SNP / non-SNP split. REF=A ALT=G,* stays SNP (htsjdk ignores '*' when typing); a sample called '*' becomes N in phase 3.
  if ! gatk --java-options "$jopts" SelectVariants -R "$REF_FA" -V "$in" --select-type-to-include SNP -O "$snp_raw" >>"$log" 2>&1; then
    filt_fail "$cs" "SelectVariants (SNP) failed"; return 1
  fi
  if ! gatk --java-options "$jopts" SelectVariants -R "$REF_FA" -V "$in" --select-type-to-exclude SNP -O "$indel_raw" >>"$log" 2>&1; then
    filt_fail "$cs" "SelectVariants (non-SNP) failed"; return 1
  fi
  if ! n_raw=$(count_records "$in") || ! n_snp=$(count_records "$snp_raw") || ! n_nonsnp=$(count_records "$indel_raw"); then
    filt_fail "$cs" "record counting failed"; return 1
  fi
  if (( n_raw != n_snp + n_nonsnp )); then
    filt_fail "$cs" "SNP/non-SNP split is not a partition: raw=${n_raw} snp=${n_snp} nonsnp=${n_nonsnp}"; return 1
  fi
  printf 'n_raw_records\t%s\nn_snp_records\t%s\nn_nonsnp_records_saved_raw\t%s\n' "$n_raw" "$n_snp" "$n_nonsnp" > "$od/.phase1_counts.tsv"

  # 2) Annotation missing rates over all SNP records (before filtering).
  if ! bcftools query -f '%QUAL\t%INFO/MQ\t%INFO/QD\t%INFO/FS\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\t%INFO/SOR\n' "$snp_raw" 2>>"$log" \
     | awk -F'\t' -v n="$n_snp" '
         BEGIN { split("QUAL MQ QD FS MQRankSum ReadPosRankSum SOR", h, " ") }
         { for (i = 1; i <= 7; i++) if ($i == "." || $i == "") miss[i]++ }
         END {
           printf "annotation\tn_snp_records\tn_missing\tmissing_rate\n"
           for (i = 1; i <= 7; i++) printf "%s\t%d\t%d\t%.6f\n", h[i], n, miss[i] + 0, (n > 0 ? (miss[i] + 0) / n : 0)
         }' > "$od/annotation_missing.tsv"; then
    filt_fail "$cs" "annotation missing-rate report failed"; return 1
  fi

  # 3) Seven independent named site filters (fail-if expressions; floating-point literals).
  if ! gatk --java-options "$jopts" VariantFiltration -R "$REF_FA" -V "$snp_raw" -O "$sitefilt" "${VF_ARGS[@]}" >>"$log" 2>&1; then
    filt_fail "$cs" "VariantFiltration failed"; return 1
  fi

  # 4) Per-filter counts: n_fail_any is non-exclusive; n_fail_only_this = failed this filter and nothing else.
  if ! bcftools query -f '%FILTER\n' "$sitefilt" 2>>"$log" \
     | awk -F';' -v names="$(printf '%s,' "${FILTER_NAMES[@]}")" '
         { n++; if ($0 == "PASS") { p++; next }
           for (i = 1; i <= NF; i++) c[$i]++
           if (NF == 1) o[$1]++ }
         END {
           k = split(names, nm, ",")
           printf "filter\tn_fail_any\tn_fail_only_this\n"
           for (i = 1; i <= k; i++) if (nm[i] != "") printf "%s\t%d\t%d\n", nm[i], c[nm[i]] + 0, o[nm[i]] + 0
           printf "PASS\t%d\t-\n", p + 0
           printf "TOTAL_SNP_RECORDS\t%d\t-\n", n + 0
         }' > "$od/filter_counts.tsv"; then
    filt_fail "$cs" "filter count report failed"; return 1
  fi

  # 5) PASS -> drop static mask -> pass_unmasked + site list for gt_source.py.
  rm -f "$passmask" "${passmask}.tbi"
  if ! bcftools view -f PASS -T "^${MASK_BED}" --threads "$BCF_THREADS" -Oz -o "$passmask" "$sitefilt" 2>>"$log" || ! bcftools index --threads "$BCF_THREADS" -t "$passmask" 2>>"$log"; then
    filt_fail "$cs" "PASS / static-mask step failed"; return 1
  fi
  if ! bcftools query -f '%CHROM\t%POS\n' "$passmask" > "$od/sites.tsv" 2>>"$log"; then
    filt_fail "$cs" "site list failed"; return 1
  fi
  touch "$od/.phase1.ok"
  echo "[$(date '+%F %T')] >>> DONE  phase 1: ${cs}  raw=${n_raw} snp=${n_snp} pass_unmasked=$(awk 'END{print NR}' "$od/sites.tsv")" | tee -a "$log"
  update_progress "$cs"
}

# Phase 2: genotype source per sample
source_job() {   # <callset> <sample>
  local cs="$1" s="$2"
  local od="$FILT_OUT/$cs"
  local out="$od/gt_source/${s}.tsv.gz"
  local log="$od/logs/gt_source.${s}.log"
  local g="$GVCF_OUT/${s}.g.vcf.gz" cb="$CALL_OUT/callable/${s}.cna3.bed.gz"
  rm -f "$out" "${out}.tmp"
  if [[ ! -s "$g" || ! -s "$cb" ]]; then
    echo "[ERROR] gVCF or base 2 missing: $s" | tee -a "$log"; echo "gVCF or base 2 missing" > "$FAIL_DIR/src__${s}"; update_progress "$s"; return 1
  fi
  if ! python3 "$PIPELINE_DIR/gt_source.py" --gvcf "$g" --sites "$od/sites.tsv" --callable "$cb" --genome "$GENOME_TXT" --out "${out}.tmp" > "$log" 2>&1; then
    echo "[ERROR] gt_source.py failed: $s (see $log)"; echo "gt_source.py failed" > "$FAIL_DIR/src__${s}"; rm -f "${out}.tmp"; update_progress "$s"; return 1
  fi
  mv -f "${out}.tmp" "$out"
  update_progress "$s"
}

# Phase 3: genotype-level N, one callset
finish_callset() {
  local cs="$1" prefix od log
  prefix="$(callset_prefix "$cs")"
  od="$FILT_OUT/$cs"
  log="$od/logs/${prefix}.filter.log"
  local passmask="$od/${prefix}.cna3.snp.pass_unmasked.vcf.gz"
  local sitefilt="$od/${prefix}.cna3.snp.sitefilt.vcf.gz"
  local base1="$od/${prefix}.cna3.snp.filtered.vcf.gz"
  local base1_tmp="$od/${prefix}.cna3.snp.filtered.tmp.vcf.gz"
  local n_pass n_base1 got exp n_samples rc i s

  echo "[$(date '+%F %T')] >>> START phase 3: ${cs}" | tee -a "$log"
  for i in "${!S_NAME[@]}"; do
    s="${S_NAME[$i]}"
    if [[ "${S_CALLSET[$i]}" == "$cs" && ! -s "$od/gt_source/${s}.tsv.gz" ]]; then
      filt_fail "$cs" "genotype source missing for $s (phase 2)"; return 1
    fi
  done

  # Genotype-level N with rule A -> base 1.
  rm -f "$base1_tmp" "${base1_tmp}.tbi"
  if ! bcftools view -Ov "$passmask" 2>>"$log" \
     | python3 "$PIPELINE_DIR/gt_mask.py" --gq-min "$GT_GQ_MIN" --dp-min "$GT_DP_MIN" --ad-frac-min "$GT_AD_FRAC_MIN" --src-dir "$od/gt_source" --kept-out "$od/refblock_kept.tsv.gz" --stats "$od/gt_mask_stats.tsv" --site-stats "$od/gt_mask_sites.tsv" 2>>"$log" \
     | bcftools view --threads "$BCF_THREADS" -Oz -o "$base1_tmp" 2>>"$log"; then
    filt_fail "$cs" "genotype-N pipeline failed"; return 1
  fi
  if ! bcftools index --threads "$BCF_THREADS" -t "$base1_tmp" 2>>"$log"; then
    filt_fail "$cs" "bcftools index failed"; return 1
  fi
  mv -f "$base1_tmp" "$base1"; mv -f "${base1_tmp}.tbi" "${base1}.tbi"

  if ! n_pass=$(bcftools view -H -f PASS "$sitefilt" | awk 'END{print NR}') || ! n_base1=$(count_records "$base1"); then
    filt_fail "$cs" "record counting failed"; return 1
  fi
  got="$(bcftools query -l "$base1" 2>/dev/null | LC_ALL=C sort | tr '\n' ' ' || true)"
  exp="$(for i in "${!S_NAME[@]}"; do if [[ "${S_CALLSET[$i]}" == "$cs" ]]; then echo "${S_NAME[$i]}"; fi; done | LC_ALL=C sort | tr '\n' ' ')"
  if [[ "$got" != "$exp" ]]; then
    filt_fail "$cs" "sample list of base 1 differs from the roster"; return 1
  fi
  n_samples=$(printf '%s' "$got" | wc -w | tr -d ' ')

  {
    printf 'metric\tvalue\n'
    printf 'callset\t%s\n' "$cs"
    printf 'n_samples\t%s\n' "$n_samples"
    cat "$od/.phase1_counts.tsv"
    printf 'n_snp_pass_site_filters\t%s\n' "$n_pass"
    printf 'n_pass_removed_by_static_mask\t%s\n' "$(( n_pass - n_base1 ))"
    printf 'n_base1_sites\t%s\n' "$n_base1"
    awk -F'\t' 'NR > 1 {print "gt_mask_" $1 "\t" $2}' "$od/gt_mask_sites.tsv"
  } > "$od/filter_summary.tsv"

  echo "[$(date '+%F %T')] >>> DONE  phase 3: ${cs}  pass=${n_pass} base1=${n_base1} rule-A kept=$(awk -F'\t' '$1=="genotypes_n_ruleA_kept"{print $2}' "$od/gt_mask_sites.tsv") -> $base1" | tee -a "$log"
  echo "[$(date '+%F %T')] >>>       AD evidence (FORMAT/ADV, all genotypes): $(awk -F'\t' '$1 ~ /^genotypes_n_adv_/ {sub("genotypes_n_adv_", "", $1); printf "%s=%s ", $1, $2}' "$od/gt_mask_sites.tsv")" | tee -a "$log"
  n_adv_x=$(awk -F'\t' '$1=="genotypes_n_adv_X"{print $2}' "$od/gt_mask_sites.tsv")
  if [[ "${n_adv_x:-0}" != "0" ]]; then
    echo "[WARN] ${cs}: ${n_adv_x} genotypes whose AD differs from the sample gVCF variant record (ADV=X); they stay out of the intermediate-fraction denominator, see gt_mask_stats.tsv" | tee -a "$log"
  fi
  update_progress "$cs"
}

# Run phase 1
echo "[$(date '+%F %T')] >>> phase 1: site filters + static mask (${#CALLSETS[@]} callsets)"
progress_reset "${#CALLSETS[@]}"; open_fifo "${#CALLSETS[@]}"
for cs in "${CALLSETS[@]}"; do run_parallel prep_callset "$cs"; done
close_fifo
OK_CS=()
for cs in "${CALLSETS[@]}"; do
  if [[ -e "$FILT_OUT/$cs/.phase1.ok" ]]; then OK_CS+=( "$cs" ); fi
done

# Run phase 2
SRC_TASKS=()
for i in "${!S_NAME[@]}"; do
  for cs in ${OK_CS[@]+"${OK_CS[@]}"}; do
    if [[ "${S_CALLSET[$i]}" == "$cs" ]]; then SRC_TASKS+=( "${cs}:${S_NAME[$i]}" ); fi
  done
done
if (( ${#SRC_TASKS[@]} > 0 )); then
  for cs in "${OK_CS[@]}"; do mkdir -p "$FILT_OUT/$cs/gt_source"; done
  WIDTH="$(fifo_width "${FILTER_MAX_JOBS:-}" "$THREADS_TOTAL" 1 "${#SRC_TASKS[@]}")"
  echo "[$(date '+%F %T')] >>> phase 2: genotype source for ${#SRC_TASKS[@]} samples, ${WIDTH} in parallel"
  progress_reset "${#SRC_TASKS[@]}"; open_fifo "$WIDTH"
  for t in "${SRC_TASKS[@]}"; do run_parallel source_job "${t%%:*}" "${t#*:}"; done
  close_fifo
fi

# Run phase 3
if (( ${#OK_CS[@]} > 0 )); then
  echo "[$(date '+%F %T')] >>> phase 3: genotype-level N (rule A) (${#OK_CS[@]} callsets)"
  progress_reset "${#OK_CS[@]}"; open_fifo "${#OK_CS[@]}"
  for cs in "${OK_CS[@]}"; do run_parallel finish_callset "$cs"; done
  close_fifo
fi
rm -f "$COUNTER_FILE" "$LOCK_FILE"

# Failure report + summaries
n_ok=0
for cs in "${CALLSETS[@]}"; do
  prefix="$(callset_prefix "$cs")"
  if [[ -s "$FILT_OUT/$cs/${prefix}.cna3.snp.filtered.vcf.gz" && -s "$FILT_OUT/$cs/filter_summary.tsv" && ! -e "$FAIL_DIR/$cs" ]]; then
    n_ok=$((n_ok + 1))
    echo "---- ${cs}: filter_summary.tsv"; sed 's/^/    /' "$FILT_OUT/$cs/filter_summary.tsv"
    echo "---- ${cs}: annotation_missing.tsv"; sed 's/^/    /' "$FILT_OUT/$cs/annotation_missing.tsv"
    echo "---- ${cs}: filter_counts.tsv"; sed 's/^/    /' "$FILT_OUT/$cs/filter_counts.tsv"
  elif [[ ! -e "$FAIL_DIR/$cs" ]]; then
    echo "no output and no error recorded (job killed?)" > "$FAIL_DIR/$cs"
  fi
done
write_failed_list "$FAIL_DIR" "$FILT_OUT/FAILED_SAMPLES.txt"
if (( n_ok != ${#CALLSETS[@]} )); then
  echo "=============================================================="
  echo "[ERROR] $(( ${#CALLSETS[@]} - n_ok )) callset(s) FAILED filtering (list: $FILT_OUT/FAILED_SAMPLES.txt):" >&2
  for f in "$FAIL_DIR"/*; do printf '   %s  (%s)\n' "$(basename "$f")" "$(head -n 1 "$f")" >&2; done
  echo "=============================================================="
  exit 1
fi
echo "[$(date '+%F %T')] >>> filtering complete. Base 1: $FILT_OUT/main/st5_joint.cna3.snp.filtered.vcf.gz (+ h99, case)"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
