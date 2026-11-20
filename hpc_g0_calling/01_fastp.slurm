#!/bin/bash
#SBATCH -J g0_01_fastp
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=80
#SBATCH -o 01_fastp.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
FAIL_DIR="$FASTP_OUT/.failed"
mkdir -p "$FASTP_OUT" "$FAIL_DIR" "$RUNINFO_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$FASTP_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner

require_tools fastp || exit 1
log_versions 01_fastp fastp
[[ -d "$RAW_ROOT" ]] || { echo "[ERROR] RAW_ROOT not found: $RAW_ROOT" >&2; exit 1; }

# Samples (static manifest)
load_samples || exit 1
SAMPLES=(); R1_LIST=(); R2_LIST=()
for i in "${!S_NAME[@]}"; do
  SAMPLES+=( "${S_NAME[$i]}" )
  R1_LIST+=( "$RAW_ROOT/${S_R1[$i]}" )
  R2_LIST+=( "$RAW_ROOT/${S_R2[$i]}" )
done
printf '%s\n' "${SAMPLES[@]}" > "$RUNINFO_DIR/01_fastp.samples.txt"
TOTAL_JOBS="${#SAMPLES[@]}"
echo "[$(date '+%F %T')] >>> ${TOTAL_JOBS} samples for fastp."

# Per-job thread budget: -w worker threads + reader/writer threads (FASTP_IO_CORES)
THREADS_PER_JOB="$FASTP_THREADS"
(( THREADS_PER_JOB < 1 )) && THREADS_PER_JOB=1
# at most 16 worker threads per fastp job
(( THREADS_PER_JOB > 16 )) && THREADS_PER_JOB=16
MAX_JOBS="$(fifo_width "${FASTP_MAX_JOBS:-}" "$THREADS_TOTAL" "$(( THREADS_PER_JOB + FASTP_IO_CORES ))" "$TOTAL_JOBS")"
echo "[$(date '+%F %T')] >>> CPUs: $THREADS_TOTAL  parallel jobs: $MAX_JOBS  fastp -w $THREADS_PER_JOB (+ ${FASTP_IO_CORES} reader/writer cores per job)"

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

# fastp task
fastp_fail() {   # <sample> <reason>
  local sample="$1" reason="$2" odir="$FASTP_OUT/$1"
  echo "[ERROR] ${reason}: $sample"
  echo "$reason" > "$FAIL_DIR/$sample"
  # never leave truncated reads that stage 02 could pick up
  rm -f "$odir/${sample}_clean_R1.fq.gz" "$odir/${sample}_clean_R2.fq.gz" "$odir/${sample}.json"
  update_progress "$sample"
}

fastp_job() {
  local sample="$1" r1="$2" r2="$3"
  local odir="$FASTP_OUT/$sample"; mkdir -p "$odir"
  local log="$odir/${sample}.fastp.log"
  echo "[$(date '+%F %T')] >>> START fastp: $sample"

  # To resume a partial run, uncomment to skip completed samples:
  # if [[ -s "$odir/${sample}_clean_R1.fq.gz" && -s "$odir/${sample}_clean_R2.fq.gz" && -s "$odir/${sample}.json" ]]; then
  #   echo "[$(date '+%F %T')] >>> skip (done): $sample"; update_progress "$sample"; return 0
  # fi

  if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    fastp_fail "$sample" "raw FASTQ missing ($r1 | $r2)"; return 1
  fi
  if ! fastp \
        -i "$r1" -I "$r2" \
        -o "$odir/${sample}_clean_R1.fq.gz" \
        -O "$odir/${sample}_clean_R2.fq.gz" \
        -j "$odir/${sample}.json" \
        -h "$odir/${sample}.html" \
        -w "$THREADS_PER_JOB" >"$log" 2>&1; then
    fastp_fail "$sample" "fastp failed (see $log)"; return 1
  fi
  echo "[$(date '+%F %T')] >>> DONE  fastp: $sample"
  update_progress "$sample"
}

# Submit
for i in "${!SAMPLES[@]}"; do
  run_parallel fastp_job "${SAMPLES[$i]}" "${R1_LIST[$i]}" "${R2_LIST[$i]}"
done
wait || true
exec 3>&- 3<&- || true
rm -f "$COUNTER_FILE" "$LOCK_FILE"

# Failure report
# A sample is done only with both reads and the json (written last), so a killed job cannot pass truncated reads on.
n_ok=0
for s in "${SAMPLES[@]}"; do
  if [[ -s "$FASTP_OUT/$s/${s}_clean_R1.fq.gz" && -s "$FASTP_OUT/$s/${s}_clean_R2.fq.gz" && -s "$FASTP_OUT/$s/${s}.json" ]]; then
    n_ok=$((n_ok + 1))
  elif [[ ! -e "$FAIL_DIR/$s" ]]; then
    echo "no output and no error recorded (job killed?)" > "$FAIL_DIR/$s"
  fi
done
write_failed_list "$FAIL_DIR" "$FASTP_OUT/FAILED_SAMPLES.txt"
FAILED=( "$FAIL_DIR"/* )
if (( ${#FAILED[@]} > 0 )); then
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED fastp:"
  for f in "${FAILED[@]}"; do printf '   %s  (%s)\n' "$(basename "$f")" "$(head -n 1 "$f")"; done
  echo "       full list: $FASTP_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
fi
echo "[$(date '+%F %T')] >>> fastp complete. Cleaned pairs: ${n_ok} / ${TOTAL_JOBS}"
(( n_ok == 0 )) && { echo "[ERROR] No cleaned reads produced; halting pipeline." >&2; exit 1; }

echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
