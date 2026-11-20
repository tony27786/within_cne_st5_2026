#!/bin/bash
#SBATCH -J g0_05_joint
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=1200G
#SBATCH --cpus-per-task=80
#SBATCH -o 05_joint.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
FAIL_DIR="$JOINT_OUT/.failed"
mkdir -p "$JOINT_OUT/tmp" "$FAIL_DIR" "$RUNINFO_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$JOINT_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner

require_tools gatk bcftools || exit 1
log_versions 05_joint gatk bcftools
for f in "$REF_FA" "${REF_FA}.fai" "$REF_DICT" "$GENOME_TXT" "$CONTIG_LIST"; do
  [[ -s "$f" ]] || { echo "[ERROR] missing $f (00_ref_prep)" >&2; exit 1; }
done
mapfile -t CONTIGS < "$CONTIG_LIST"                                            # .fai order (merge order)
mapfile -t CONTIGS_BY_LEN < <(sort -k2,2nr "$GENOME_TXT" | awk '{print $1}')   # longest first (scheduling)

# Samples + completeness gate
load_samples || exit 1
MISSING=()
for i in "${!S_NAME[@]}"; do
  s="${S_NAME[$i]}"; g="$GVCF_OUT/${s}.g.vcf.gz"
  if [[ ! -s "$g" || ! -s "${g}.tbi" ]]; then
    MISSING+=( "$s (no gVCF)" ); continue
  fi
  sm="$(bcftools query -l "$g" 2>/dev/null || true)"
  if [[ "$sm" != "$s" ]]; then MISSING+=( "$s (gVCF sample name: ${sm:-unreadable})" ); fi
done
if (( ${#MISSING[@]} > 0 )); then
  printf '%s\n' "${MISSING[@]}" > "$JOINT_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
  echo "[ERROR] completeness gate: ${#MISSING[@]} of ${#S_NAME[@]} gVCFs missing/invalid; no callset is joint-called:" >&2
  printf '   %s\n' "${MISSING[@]}" >&2
  echo "        list: $JOINT_OUT/FAILED_SAMPLES.txt" >&2
  echo "=============================================================="
  exit 1
fi
rm -f "$JOINT_OUT/FAILED_SAMPLES.txt"
echo "[$(date '+%F %T')] >>> completeness gate passed: all ${#S_NAME[@]} gVCFs present"

# Per-callset sample maps (sample<TAB>gVCF)
for cs in "${CALLSETS[@]}"; do
  mkdir -p "$JOINT_OUT/$cs/by_contig" "$JOINT_OUT/$cs/genomicsdb" "$JOINT_OUT/$cs/logs"
  map="$JOINT_OUT/$cs/sample_map.tsv"
  : > "$map"
  for i in "${!S_NAME[@]}"; do
    if [[ "${S_CALLSET[$i]}" == "$cs" ]]; then
      printf '%s\t%s\n' "${S_NAME[$i]}" "$GVCF_OUT/${S_NAME[$i]}.g.vcf.gz" >> "$map"
    fi
  done
  n_map=$(awk 'END{print NR}' "$map")
  cp -f "$map" "$RUNINFO_DIR/05_joint.${cs}.sample_map.tsv"
  echo "[$(date '+%F %T')] >>> callset ${cs}: ${n_map} samples -> $map"
done

# Task list: main first, longest contigs first
TASKS=()
for cs in "${CALLSETS[@]}"; do
  for c in "${CONTIGS_BY_LEN[@]}"; do TASKS+=( "${cs}:${c}" ); done
done
TOTAL_JOBS="${#TASKS[@]}"
MAX_JOBS="$(fifo_width "${JOINT_MAX_JOBS:-}" "$THREADS_TOTAL" 1 "$TOTAL_JOBS")"
echo "[$(date '+%F %T')] >>> ${TOTAL_JOBS} (callset, contig) tasks; parallel: $MAX_JOBS; GenomicsDBImport --reader-threads $JOINT_READER_THREADS; heap: $JOINT_JAVA_XMX"

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

# Joint-genotyping task: one (callset, contig)
joint_fail() {   # <task> <reason>
  local task="$1" reason="$2"
  local cs="${task%%:*}" contig="${task#*:}"
  echo "[ERROR] ${reason}: ${task}" | tee -a "$JOINT_OUT/$cs/logs/${contig}.log"
  echo "$reason" > "$FAIL_DIR/${cs}__${contig}"
  rm -f "$JOINT_OUT/$cs/by_contig/${contig}.tmp.vcf.gz" "$JOINT_OUT/$cs/by_contig/${contig}.tmp.vcf.gz.tbi"
  rm -rf "$JOINT_OUT/tmp/${cs}__${contig}"
  update_progress "$task"
}

joint_task() {
  local task="$1"
  local cs="${task%%:*}" contig="${task#*:}"
  local ws="$JOINT_OUT/$cs/genomicsdb/${contig}"
  local vcf="$JOINT_OUT/$cs/by_contig/${contig}.vcf.gz"
  local vcf_tmp="$JOINT_OUT/$cs/by_contig/${contig}.tmp.vcf.gz"
  local ivl="$INTERVAL_DIR/${contig}.intervals"
  local log="$JOINT_OUT/$cs/logs/${contig}.log"; : > "$log"
  local tmp_dir="$JOINT_OUT/tmp/${cs}__${contig}"; mkdir -p "$tmp_dir"
  local jopts; jopts="$(gatk_java_opts "$JOINT_JAVA_XMX" "$tmp_dir")"

  echo "[$(date '+%F %T')] >>> START joint: ${task}" | tee -a "$log"

  # To resume a partial run, uncomment to skip completed tasks:
  # if [[ -s "$vcf" && -s "${vcf}.tbi" ]]; then
  #   echo "[$(date '+%F %T')] >>> skip (done): ${task}"; rm -rf "$tmp_dir"; update_progress "$task"; return 0
  # fi

  if [[ ! -s "$ivl" ]]; then joint_fail "$task" "interval file missing: $ivl"; return 1; fi
  if ! gatk --java-options "$jopts" GenomicsDBImport \
        --genomicsdb-workspace-path "$ws" \
        --overwrite-existing-genomicsdb-workspace true \
        --sample-name-map "$JOINT_OUT/$cs/sample_map.tsv" \
        -L "$ivl" \
        --batch-size "$GDB_BATCH_SIZE" \
        --reader-threads "$JOINT_READER_THREADS" \
        --genomicsdb-shared-posixfs-optimizations true \
        --tmp-dir "$tmp_dir" >>"$log" 2>&1; then
    joint_fail "$task" "GenomicsDBImport failed (see $log)"; return 1
  fi
  rm -f "$vcf_tmp" "${vcf_tmp}.tbi"
  if ! gatk --java-options "$jopts" GenotypeGVCFs \
        -R "$REF_FA" \
        -V "gendb://${ws}" \
        -L "$ivl" \
        -O "$vcf_tmp" \
        --genomicsdb-shared-posixfs-optimizations true \
        --tmp-dir "$tmp_dir" >>"$log" 2>&1; then
    joint_fail "$task" "GenotypeGVCFs failed (see $log)"; return 1
  fi
  if [[ ! -s "$vcf_tmp" || ! -s "${vcf_tmp}.tbi" ]]; then
    joint_fail "$task" "GenotypeGVCFs output/index missing"; return 1
  fi
  mv -f "$vcf_tmp" "$vcf"
  mv -f "${vcf_tmp}.tbi" "${vcf}.tbi"

  rm -rf "$tmp_dir"
  echo "[$(date '+%F %T')] >>> DONE  joint: ${task}" | tee -a "$log"
  update_progress "$task"
}

# Submit
for t in "${TASKS[@]}"; do
  run_parallel joint_task "$t"
done
wait || true
exec 3>&- 3<&- || true
rm -f "$COUNTER_FILE" "$LOCK_FILE"

# Merge contigs per callset + roster check
cs_failed=0
for cs in "${CALLSETS[@]}"; do
  prefix="$(callset_prefix "$cs")"
  raw="$JOINT_OUT/$cs/${prefix}.cna3.raw.vcf.gz"
  raw_tmp="$JOINT_OUT/$cs/${prefix}.cna3.raw.tmp.vcf.gz"
  bad=()
  parts=()
  for c in "${CONTIGS[@]}"; do
    part="$JOINT_OUT/$cs/by_contig/${c}.vcf.gz"
    if [[ -e "$FAIL_DIR/${cs}__${c}" || ! -s "$part" ]]; then
      bad+=( "$c" )
      if [[ ! -e "$FAIL_DIR/${cs}__${c}" ]]; then
        echo "no output and no error recorded (job killed?)" > "$FAIL_DIR/${cs}__${c}"
      fi
    fi
    parts+=( "$part" )
  done
  if (( ${#bad[@]} > 0 )); then
    echo "[ERROR] callset ${cs}: ${#bad[@]} contig task(s) failed: $(printf '%s ' "${bad[@]}")" >&2
    cs_failed=1; continue
  fi
  rm -f "$raw" "${raw}.tbi" "$raw_tmp" "${raw_tmp}.tbi"
  if ! bcftools concat --threads "$BCF_THREADS" -Oz -o "$raw_tmp" "${parts[@]}" 2>>"$JOINT_OUT/$cs/logs/concat.log" || ! bcftools index --threads "$BCF_THREADS" -t "$raw_tmp" 2>>"$JOINT_OUT/$cs/logs/concat.log"; then
    echo "[ERROR] callset ${cs}: bcftools concat/index failed" >&2
    echo "concat failed" > "$FAIL_DIR/${cs}__concat"
    cs_failed=1; continue
  fi
  mv -f "$raw_tmp" "$raw"; mv -f "${raw_tmp}.tbi" "${raw}.tbi"

  got="$(bcftools query -l "$raw" 2>/dev/null | LC_ALL=C sort | tr '\n' ' ' || true)"
  exp="$(for i in "${!S_NAME[@]}"; do if [[ "${S_CALLSET[$i]}" == "$cs" ]]; then echo "${S_NAME[$i]}"; fi; done | LC_ALL=C sort | tr '\n' ' ')"
  n_got=$(printf '%s' "$got" | wc -w | tr -d ' ')
  if [[ "$got" != "$exp" ]]; then
    echo "[ERROR] callset ${cs}: sample list of $raw differs from the roster" >&2
    echo "roster mismatch" > "$FAIL_DIR/${cs}__roster"
    cs_failed=1; continue
  fi
  n_rec=$(bcftools view -H "$raw" | awk 'END{print NR}')
  echo "[$(date '+%F %T')] >>> callset ${cs}: ${n_got} samples, ${n_rec} records -> $raw"
done

# Failure report
FAILED=( "$FAIL_DIR"/* )
if (( ${#FAILED[@]} > 0 )); then
  for f in "${FAILED[@]}"; do printf '%s\t%s\n' "$(basename "$f")" "$(head -n 1 "$f")"; done > "$JOINT_OUT/FAILED_TASKS.txt"
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} joint-genotyping task(s) FAILED (callset__contig):"
  sed 's/^/   /' "$JOINT_OUT/FAILED_TASKS.txt"
  echo "       list: $JOINT_OUT/FAILED_TASKS.txt"
  echo "=============================================================="
else
  rm -f "$JOINT_OUT/FAILED_TASKS.txt"
fi
if (( cs_failed )); then
  echo "[ERROR] at least one callset is incomplete; stopping before 06_filter." >&2
  exit 1
fi
echo "[$(date '+%F %T')] >>> joint genotyping complete: main / h99 / case raw VCFs under $JOINT_OUT/<callset>/"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
