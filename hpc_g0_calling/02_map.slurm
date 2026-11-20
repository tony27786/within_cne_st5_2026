#!/bin/bash
#SBATCH -J g0_02_map
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=1000G
#SBATCH --cpus-per-task=80
#SBATCH -o 02_map.out

# Load shared config (PIPELINE_DIR = where sbatch was invoked)
PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Output layout
METRICS_OUT="$MAP_OUT/qc_metrics_map"
FLAGSTAT_DIR="$METRICS_OUT/flagstat"
STATS_DIR="$METRICS_OUT/stats"
IDXSTATS_DIR="$METRICS_OUT/idxstats"
MARKDUP_DIR="$METRICS_OUT/markdup"
LOG_DIR="$METRICS_OUT/logs"
FAIL_DIR="$MAP_OUT/.failed"
mkdir -p "$MAP_OUT" "$FLAGSTAT_DIR" "$STATS_DIR" "$IDXSTATS_DIR" "$MARKDUP_DIR" "$LOG_DIR" "$FAIL_DIR" "$RUNINFO_DIR"

cd "$MAP_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner
resolve_shard || exit 1

# Dependencies & reference (indexes are built by 00_ref_prep)
require_tools bwa-mem2 samtools || exit 1
log_versions "02_map.shard${SHARD_ID}" bwa-mem2 samtools
[[ -f "$REF_FA" ]] || { echo "[ERROR] Reference FASTA not found: $REF_FA" >&2; exit 1; }
for ext in fai 0123 amb ann bwt.2bit.64 pac; do
  [[ -s "${REF_FA}.${ext}" ]] || { echo "[ERROR] missing ${REF_FA}.${ext} (00_ref_prep builds it)" >&2; exit 1; }
done

# Samples of this shard (static manifest; never a directory scan)
load_samples || exit 1
mapfile -t IDX < <(shard_indices)
SAMPLES=()
for i in "${IDX[@]}"; do SAMPLES+=( "${S_NAME[$i]}" ); done
TOTAL_JOBS="${#SAMPLES[@]}"
(( TOTAL_JOBS == 0 )) && { echo "[ERROR] shard ${SHARD_ID} has no samples." >&2; exit 1; }
for s in "${SAMPLES[@]}"; do rm -f "$FAIL_DIR/$s"; done   # only this shard's markers
printf '%s\n' "${SAMPLES[@]}" > "$RUNINFO_DIR/02_map.shard${SHARD_ID}.samples.txt"
echo "[$(date '+%F %T')] >>> ${TOTAL_JOBS} samples for mapping in shard ${SHARD_ID}."

# Per-job thread budget
JOB_CPUS="$MAP_JOB_CORES"; (( JOB_CPUS > THREADS )) && JOB_CPUS="$THREADS"; (( JOB_CPUS < 1 )) && JOB_CPUS=1
MAX_JOBS="$(fifo_width "${MAP_MAX_JOBS:-}" "$THREADS" "$JOB_CPUS" "$TOTAL_JOBS")"
BWA_THREADS=$(( JOB_CPUS - 1 )); (( BWA_THREADS < 1 )) && BWA_THREADS=1   # 1 core left for fixmate + sort intake in the same pipe
POST_THREADS="$MAP_POST_THREADS"; (( POST_THREADS > JOB_CPUS )) && POST_THREADS="$JOB_CPUS"
THREADS_PER_JOB="$JOB_CPUS"
echo "[$(date '+%F %T')] >>> Total CPUs: $THREADS; parallel jobs: $MAX_JOBS; CPUs/job: $JOB_CPUS"
echo "[$(date '+%F %T')] >>> Per-job: bwa-mem2 -t $BWA_THREADS; samtools sort/markdup/index/stats -@ $POST_THREADS; sort -m $MAP_SORT_MEM per thread"

# FIFO semaphore
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() {
  read -u 3
  (
    trap 'echo >&3' EXIT
    export OMP_NUM_THREADS="$THREADS_PER_JOB" THREADS="$THREADS_PER_JOB"
    "$@"
  ) &
}

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

# Mapping task
map_fail() {   # <sample> <reason>
  local sample="$1" reason="$2"
  echo "[ERROR] ${reason}: $sample" | tee -a "$LOG_DIR/${sample}.map.log"
  echo "$reason" > "$FAIL_DIR/$sample"
  rm -f "$MAP_OUT/${sample}.md.tmp.bam" "$MAP_OUT/${sample}.md.bam" "$MAP_OUT/${sample}.md.bam.bai"
  rm -rf "$MAP_OUT/tmp/${sample}"
  update_progress "$sample"
}

map_job() {
  local sample="$1"
  local r1="$FASTP_OUT/${sample}/${sample}_clean_R1.fq.gz"
  local r2="$FASTP_OUT/${sample}/${sample}_clean_R2.fq.gz"
  local fastp_json="$FASTP_OUT/${sample}/${sample}.json"
  local bam_out="$MAP_OUT/${sample}.md.bam"
  local bam_tmp="$MAP_OUT/${sample}.md.tmp.bam"
  local flag_out="$FLAGSTAT_DIR/${sample}.flagstat.txt"
  local stats_out="$STATS_DIR/${sample}.stats.txt"
  local idxstats_out="$IDXSTATS_DIR/${sample}.idxstats.txt"
  local markdup_out="$MARKDUP_DIR/${sample}.markdup.txt"
  local log_out="$LOG_DIR/${sample}.map.log"; : > "$log_out"
  local rg="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA"
  local tmp_dir="$MAP_OUT/tmp/${sample}"; mkdir -p "$tmp_dir"

  echo "[$(date '+%F %T')] >>> START map: $sample" | tee -a "$log_out"

  # To resume a partial run, uncomment to skip completed samples:
  # if [[ -s "$bam_out" && -s "${bam_out}.bai" && -s "$flag_out" ]]; then
  #   echo "[$(date '+%F %T')] >>> skip (done): $sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 0
  # fi

  if [[ ! -s "$r1" || ! -s "$r2" || ! -s "$fastp_json" ]]; then
    map_fail "$sample" "clean reads missing (stage 01)"; return 1
  fi
  rm -f "$bam_out" "${bam_out}.bai" "$bam_tmp"

  if ! bwa-mem2 mem -t "$BWA_THREADS" -R "$rg" "$REF_FA" "$r1" "$r2" 2>>"$log_out" \
       | samtools fixmate -m -u - - 2>>"$log_out" \
       | samtools sort -u -@ "$POST_THREADS" -m "$MAP_SORT_MEM" -T "${tmp_dir}/${sample}.sort" - 2>>"$log_out" \
       | samtools markdup -@ "$POST_THREADS" -T "${tmp_dir}/${sample}.markdup" -f "$markdup_out" - "$bam_tmp" 2>>"$log_out"; then
    map_fail "$sample" "bwa-mem2 | fixmate | sort | markdup failed"; return 1
  fi
  mv -f "$bam_tmp" "$bam_out"
  if ! samtools index -@ "$POST_THREADS" "$bam_out" 2>>"$log_out"; then
    map_fail "$sample" "samtools index failed"; return 1
  fi
  if ! samtools quickcheck "$bam_out"; then
    map_fail "$sample" "samtools quickcheck failed"; return 1
  fi

  samtools flagstat -@ "$POST_THREADS" "$bam_out" > "$flag_out"     2>>"$log_out" || true
  samtools stats    -@ "$POST_THREADS" "$bam_out" > "$stats_out"    2>>"$log_out" || true
  samtools idxstats                    "$bam_out" > "$idxstats_out" 2>>"$log_out" || true

  rm -rf "$tmp_dir"
  echo "[$(date '+%F %T')] >>> DONE  map: $sample" | tee -a "$log_out"
  update_progress "$sample"
}

# Submit
for s in "${SAMPLES[@]}"; do
  run_parallel map_job "$s"
done
wait || true
exec 3>&- 3<&- || true
rm -f "$COUNTER_FILE" "$LOCK_FILE"

# Failure report (this shard)
n_ok=0
for s in "${SAMPLES[@]}"; do
  if [[ -s "$MAP_OUT/${s}.md.bam" && -s "$MAP_OUT/${s}.md.bam.bai" ]]; then
    n_ok=$((n_ok + 1))
  elif [[ ! -e "$FAIL_DIR/$s" ]]; then
    echo "no output and no error recorded (job killed?)" > "$FAIL_DIR/$s"
  fi
done
write_failed_list "$FAIL_DIR" "$MAP_OUT/FAILED_SAMPLES.txt"
FAILED=()
for s in "${SAMPLES[@]}"; do
  if [[ -e "$FAIL_DIR/$s" ]]; then FAILED+=( "$s" ); fi
done
if (( ${#FAILED[@]} > 0 )); then
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED mapping in shard ${SHARD_ID}:"
  for s in "${FAILED[@]}"; do printf '   %s  (%s)\n' "$s" "$(head -n 1 "$FAIL_DIR/$s")"; done
  echo "       full list (all shards): $MAP_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
fi
echo "[$(date '+%F %T')] >>> Mapping complete (shard ${SHARD_ID}). BAMs produced: ${n_ok} / ${TOTAL_JOBS}"
(( n_ok == 0 )) && { echo "[ERROR] No BAM produced in this shard; halting pipeline." >&2; exit 1; }

echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
