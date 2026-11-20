#!/bin/bash
#SBATCH -J g0_04_hc
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=1200G
#SBATCH --cpus-per-task=80
#SBATCH -o 04_hc_gvcf.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
LOG_DIR="$GVCF_OUT/logs"
FAIL_DIR="$GVCF_OUT/.failed"
mkdir -p "$GVCF_OUT" "$LOG_DIR" "$GVCF_OUT/tmp" "$FAIL_DIR" "$RUNINFO_DIR"

cd "$GVCF_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner
resolve_shard || exit 1

require_tools gatk tabix || exit 1
log_versions "04_hc_gvcf.shard${SHARD_ID}" gatk tabix
for f in "$REF_FA" "${REF_FA}.fai" "$REF_DICT" "$GENOME_TXT"; do
  [[ -s "$f" ]] || { echo "[ERROR] missing $f (00_ref_prep)" >&2; exit 1; }
done
# Every reference contig must show up in each gVCF index (completeness check).
EXPECTED_CONTIG_SET="$(awk '{print $1}' "$GENOME_TXT" | LC_ALL=C sort | tr '\n' ' ')"

# Samples of this shard (static manifest)
load_samples || exit 1
mapfile -t IDX < <(shard_indices)
SAMPLES=()
for i in "${IDX[@]}"; do SAMPLES+=( "${S_NAME[$i]}" ); done
TOTAL_JOBS="${#SAMPLES[@]}"
(( TOTAL_JOBS == 0 )) && { echo "[ERROR] shard ${SHARD_ID} has no samples." >&2; exit 1; }
for s in "${SAMPLES[@]}"; do rm -f "$FAIL_DIR/$s"; done   # only this shard's markers
printf '%s\n' "${SAMPLES[@]}" > "$RUNINFO_DIR/04_hc_gvcf.shard${SHARD_ID}.samples.txt"
echo "[$(date '+%F %T')] >>> ${TOTAL_JOBS} samples for HaplotypeCaller in shard ${SHARD_ID}."

# Thread budget
MAX_JOBS="$(fifo_width "${HC_MAX_JOBS:-}" "$THREADS_TOTAL" "$HC_PAIRHMM_THREADS" "$TOTAL_JOBS")"
echo "[$(date '+%F %T')] >>> CPUs: $THREADS_TOTAL  parallel HC: $MAX_JOBS  pair-HMM threads: $HC_PAIRHMM_THREADS  heap: $HC_JAVA_XMX"

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

# HaplotypeCaller task
hc_fail() {   # <sample> <reason>
  local sample="$1" reason="$2"
  echo "[ERROR] ${reason}: $sample" | tee -a "$LOG_DIR/${sample}.hc.log"
  echo "$reason" > "$FAIL_DIR/$sample"
  rm -f "$GVCF_OUT/${sample}.tmp.g.vcf.gz" "$GVCF_OUT/${sample}.tmp.g.vcf.gz.tbi"
  rm -rf "$GVCF_OUT/tmp/${sample}"
  update_progress "$sample"
}

hc_job() {
  local sample="$1"
  local bam="$MAP_OUT/${sample}.md.bam"
  local out="$GVCF_OUT/${sample}.g.vcf.gz"
  local tmp_out="$GVCF_OUT/${sample}.tmp.g.vcf.gz"
  local log="$LOG_DIR/${sample}.hc.log"; : > "$log"
  local tmp_dir="$GVCF_OUT/tmp/${sample}"; mkdir -p "$tmp_dir"
  local seen

  echo "[$(date '+%F %T')] >>> START HC: $sample" | tee -a "$log"

  # To resume a partial run, uncomment to skip completed samples:
  # if [[ -s "$out" && -s "${out}.tbi" ]]; then
  #   echo "[$(date '+%F %T')] >>> skip (done): $sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 0
  # fi

  if [[ ! -s "$bam" || ! -s "${bam}.bai" ]]; then
    hc_fail "$sample" "BAM/BAI missing (stage 02)"; return 1
  fi
  rm -f "$tmp_out" "${tmp_out}.tbi"
  if ! gatk --java-options "$(gatk_java_opts "$HC_JAVA_XMX" "$tmp_dir")" HaplotypeCaller \
        -R "$REF_FA" \
        -I "$bam" \
        -O "$tmp_out" \
        -ERC GVCF \
        --sample-ploidy 1 \
        --pcr-indel-model "$HC_PCR_INDEL_MODEL" \
        --native-pair-hmm-threads "$HC_PAIRHMM_THREADS" \
        --tmp-dir "$tmp_dir" >>"$log" 2>&1; then
    hc_fail "$sample" "HaplotypeCaller failed (see $log)"; return 1
  fi
  if [[ ! -s "${tmp_out}.tbi" ]]; then
    hc_fail "$sample" "gVCF index missing"; return 1
  fi
  seen="$(tabix -l "$tmp_out" 2>>"$log" | LC_ALL=C sort | tr '\n' ' ' || true)"
  if [[ "$seen" != "$EXPECTED_CONTIG_SET" ]]; then
    hc_fail "$sample" "gVCF does not cover all reference contigs (got: ${seen:-none})"; return 1
  fi
  mv -f "$tmp_out" "$out"
  mv -f "${tmp_out}.tbi" "${out}.tbi"

  rm -rf "$tmp_dir"
  echo "[$(date '+%F %T')] >>> DONE  HC: $sample" | tee -a "$log"
  update_progress "$sample"
}

# Submit
for s in "${SAMPLES[@]}"; do
  run_parallel hc_job "$s"
done
wait || true
exec 3>&- 3<&- || true
rm -f "$COUNTER_FILE" "$LOCK_FILE"

# Failure report (this shard)
n_ok=0
for s in "${SAMPLES[@]}"; do
  if [[ -s "$GVCF_OUT/${s}.g.vcf.gz" && -s "$GVCF_OUT/${s}.g.vcf.gz.tbi" ]]; then
    n_ok=$((n_ok + 1))
  elif [[ ! -e "$FAIL_DIR/$s" ]]; then
    echo "no output and no error recorded (job killed?)" > "$FAIL_DIR/$s"
  fi
done
write_failed_list "$FAIL_DIR" "$GVCF_OUT/FAILED_SAMPLES.txt"
FAILED=()
for s in "${SAMPLES[@]}"; do
  if [[ -e "$FAIL_DIR/$s" ]]; then FAILED+=( "$s" ); fi
done
if (( ${#FAILED[@]} > 0 )); then
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED HaplotypeCaller in shard ${SHARD_ID}:"
  for s in "${FAILED[@]}"; do printf '   %s  (%s)\n' "$s" "$(head -n 1 "$FAIL_DIR/$s")"; done
  echo "       full list (all shards): $GVCF_OUT/FAILED_SAMPLES.txt"
  echo "       NOTE: 05_joint refuses to run until every gVCF exists."
  echo "=============================================================="
fi
echo "[$(date '+%F %T')] >>> HaplotypeCaller complete (shard ${SHARD_ID}). gVCFs: ${n_ok} / ${TOTAL_JOBS}"
(( n_ok == 0 )) && { echo "[ERROR] No gVCF produced in this shard; halting pipeline." >&2; exit 1; }

echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
