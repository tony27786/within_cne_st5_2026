#!/bin/bash
#SBATCH -J g0_03_callable
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=80
#SBATCH -o 03_callable.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
CALLABLE_DIR="$CALL_OUT/callable"
SUMMARY_DIR="$CALL_OUT/summary"
MOS_DIR="$CALL_OUT/mosdepth"
LOG_DIR="$CALL_OUT/logs"
FAIL_DIR="$CALL_OUT/.failed"
mkdir -p "$CALLABLE_DIR" "$SUMMARY_DIR" "$MOS_DIR" "$LOG_DIR" "$CALL_OUT/tmp" "$FAIL_DIR" "$RUNINFO_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$CALL_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner

require_tools mosdepth bedtools bgzip tabix || exit 1
log_versions 03_callable mosdepth bedtools bgzip tabix
[[ -s "$ANALYZABLE_BED" ]] || { echo "[ERROR] missing $ANALYZABLE_BED (00_ref_prep)" >&2; exit 1; }
ANALYZABLE_BP=$(awk '{s += $3 - $2} END{print s + 0}' "$ANALYZABLE_BED")
echo "[$(date '+%F %T')] >>> mosdepth -Q ${CALLABLE_MAPQ_MIN} -F ${MOSDEPTH_EXCLUDE_FLAGS} --quantize 0:${CALLABLE_DP_MIN}: (no --fast-mode)"

# Samples (static manifest)
load_samples || exit 1
SAMPLES=( "${S_NAME[@]}" )
printf '%s\n' "${SAMPLES[@]}" > "$RUNINFO_DIR/03_callable.samples.txt"
TOTAL_JOBS="${#SAMPLES[@]}"

# Thread budget
MAX_JOBS="$(fifo_width "${MOS_MAX_JOBS:-}" "$THREADS_TOTAL" "$MOS_JOB_CORES" "$TOTAL_JOBS")"
echo "[$(date '+%F %T')] >>> CPUs: $THREADS_TOTAL  parallel jobs: $MAX_JOBS  (mosdepth -t $MOS_THREADS)"

# FIFO semaphore
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

# Progress bar
PROGRESS_WIDTH=45
COUNTER_FILE="$(mktemp)"; LOCK_FILE="$(mktemp)"; echo 0 > "$COUNTER_FILE"
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

# Callable-mask task
mos_fail() {   # <sample> <reason>
  local sample="$1" reason="$2"
  echo "[ERROR] ${reason}: $sample" | tee -a "$LOG_DIR/${sample}.log"
  echo "$reason" > "$FAIL_DIR/$sample"
  rm -f "$CALLABLE_DIR/${sample}.cna3.bed.gz" "$CALLABLE_DIR/${sample}.cna3.bed.gz.tbi" "$CALLABLE_DIR/${sample}.cna3.bed.gz.tmp"
  rm -f "$SUMMARY_DIR/${sample}.callable.tsv" "$SUMMARY_DIR/${sample}.depth_by_chrom.tsv"
  rm -rf "$CALL_OUT/tmp/${sample}"
  update_progress "$sample"
}

mos_job() {
  local sample="$1"
  local bam="$MAP_OUT/${sample}.md.bam"
  local prefix="$MOS_DIR/${sample}"
  local bed_out="$CALLABLE_DIR/${sample}.cna3.bed.gz"
  local sum_out="$SUMMARY_DIR/${sample}.callable.tsv"
  local chr_out="$SUMMARY_DIR/${sample}.depth_by_chrom.tsv"
  local log="$LOG_DIR/${sample}.log"; : > "$log"
  local tmp_dir="$CALL_OUT/tmp/${sample}"; mkdir -p "$tmp_dir"
  local cb_genome cb_analyzable med_depth mean_depth

  echo "[$(date '+%F %T')] >>> START callable: $sample" | tee -a "$log"

  # To resume a partial run, uncomment to skip completed samples:
  # if [[ -s "$bed_out" && -s "${bed_out}.tbi" && -s "$sum_out" && -s "$chr_out" ]]; then
  #   echo "[$(date '+%F %T')] >>> skip (done): $sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 0
  # fi

  if [[ ! -s "$bam" || ! -s "${bam}.bai" ]]; then
    mos_fail "$sample" "BAM/BAI missing (stage 02)"; return 1
  fi
  if ! MOSDEPTH_Q0=LOW MOSDEPTH_Q1=CALLABLE TMPDIR="$tmp_dir" \
       mosdepth -t "$MOS_THREADS" -n \
                -Q "$CALLABLE_MAPQ_MIN" -F "$MOSDEPTH_EXCLUDE_FLAGS" \
                --quantize "0:${CALLABLE_DP_MIN}:" \
                --by "$ANALYZABLE_BED" \
                "$prefix" "$bam" >>"$log" 2>&1; then
    mos_fail "$sample" "mosdepth failed"; return 1
  fi
  if [[ ! -s "${prefix}.quantized.bed.gz" || ! -s "${prefix}.mosdepth.region.dist.txt" || ! -s "${prefix}.mosdepth.summary.txt" ]]; then
    mos_fail "$sample" "mosdepth outputs missing"; return 1
  fi

  # Base 2: CALLABLE runs only, lexicographically sorted, bgzip + tabix.
  if ! gzip -dc "${prefix}.quantized.bed.gz" \
       | awk 'BEGIN{OFS="\t"} $4 == "CALLABLE" {print $1, $2, $3}' \
       | LC_ALL=C sort -k1,1 -k2,2n \
       | bgzip -c > "${bed_out}.tmp" 2>>"$log"; then
    mos_fail "$sample" "callable BED extraction failed"; return 1
  fi
  mv -f "${bed_out}.tmp" "$bed_out"
  if ! tabix -f -p bed "$bed_out" 2>>"$log"; then
    mos_fail "$sample" "tabix failed"; return 1
  fi

  if ! cb_genome=$(gzip -dc "$bed_out" | awk '{s += $3 - $2} END{print s + 0}'); then
    mos_fail "$sample" "cannot read $bed_out"; return 1
  fi
  if ! cb_analyzable=$(bedtools intersect -a "$bed_out" -b "$ANALYZABLE_BED" 2>>"$log" | awk '{s += $3 - $2} END{print s + 0}'); then
    mos_fail "$sample" "bedtools intersect failed"; return 1
  fi
  # Median depth = largest d with P(depth >= d) >= 0.5 over the analyzable region (mosdepth prints that cumulative fraction with 2 decimals).
  med_depth=$(awk '$1 == "total" && $3 + 0 >= 0.5 { if ($2 + 0 > m) m = $2 + 0 } END{print m + 0}' "${prefix}.mosdepth.region.dist.txt")
  mean_depth=$(awk '$1 == "total_region" {print $4; f = 1} END{if (!f) print "NA"}' "${prefix}.mosdepth.summary.txt")
  {
    printf 'sample\tcallable_bp_genome\tcallable_bp_analyzable\tanalyzable_bp\tcallable_frac_analyzable\tmedian_depth_analyzable\tmean_depth_analyzable\n'
    awk -v s="$sample" -v g="$cb_genome" -v a="$cb_analyzable" -v t="$ANALYZABLE_BP" \
        -v md="$med_depth" -v mn="$mean_depth" \
        'BEGIN{printf "%s\t%d\t%d\t%d\t%.6f\t%s\t%s\n", s, g, a, t, (t > 0 ? a / t : 0), md, mn}'
  } > "$sum_out"
# Per-chromosome median depth over the analyzable part of each nuclear chromosome.
  {
    printf 'chrom\tmedian_depth\n'
    awk '$1 != "total" && $3 + 0 >= 0.5 { if (!($1 in m) || $2 + 0 > m[$1]) m[$1] = $2 + 0 } END{for (c in m) print c "\t" m[c]}' "${prefix}.mosdepth.region.dist.txt" | LC_ALL=C sort -k1,1
  } > "$chr_out"

  rm -rf "$tmp_dir"
  echo "[$(date '+%F %T')] >>> DONE  callable: $sample  callable_frac=$(awk -F'\t' 'NR==2{print $5}' "$sum_out")  median_depth=${med_depth}" | tee -a "$log"
  update_progress "$sample"
}

# Submit
for s in "${SAMPLES[@]}"; do
  run_parallel mos_job "$s"
done
wait || true
exec 3>&- 3<&- || true
rm -f "$COUNTER_FILE" "$LOCK_FILE"

# Failure report
n_ok=0
for s in "${SAMPLES[@]}"; do
  if [[ -s "$CALLABLE_DIR/${s}.cna3.bed.gz" && -s "$CALLABLE_DIR/${s}.cna3.bed.gz.tbi" && -s "$SUMMARY_DIR/${s}.callable.tsv" && -s "$SUMMARY_DIR/${s}.depth_by_chrom.tsv" ]]; then
    n_ok=$((n_ok + 1))
  elif [[ ! -e "$FAIL_DIR/$s" ]]; then
    echo "no output and no error recorded (job killed?)" > "$FAIL_DIR/$s"
  fi
done
write_failed_list "$FAIL_DIR" "$CALL_OUT/FAILED_SAMPLES.txt"
FAILED=( "$FAIL_DIR"/* )
if (( ${#FAILED[@]} > 0 )); then
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED the callable mask:"
  for f in "${FAILED[@]}"; do printf '   %s  (%s)\n' "$(basename "$f")" "$(head -n 1 "$f")"; done
  echo "       full list: $CALL_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
fi
echo "[$(date '+%F %T')] >>> callable masks complete: ${n_ok} / ${TOTAL_JOBS}  (base 2 under $CALLABLE_DIR/)"
(( n_ok == 0 )) && { echo "[ERROR] No callable mask produced; halting pipeline." >&2; exit 1; }

echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
