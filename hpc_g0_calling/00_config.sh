#!/usr/bin/env bash
# Shared configuration; sourced by every stage and submit_all.sh (defines variables and helpers only, no side effects).

# SLURM
SBATCH_PARTITION="<your_partition>"   # your cluster partition; keep in sync with '#SBATCH -p' in each stage
# Heavy stages take a whole node (--cpus-per-task=80); stages whose tools are single-threaded stay small.

# Conda: launcher and env prefixes (envs are activated by absolute path).
CONDA_LAUNCHER="${CONDA_LAUNCHER:-/path/to/your/conda}"   # path to your conda launcher
CONDA_ENV_MAIN="${CONDA_ENV_MAIN:-/path/to/your/envs/main}"   # env with fastp bwa-mem2 samtools gatk bcftools bedtools mosdepth python3

# Tested tool versions; a mismatch only warns (log_versions).
EXPECT_VERSIONS=( "fastp=1.3.7" "bwa-mem2=2.3" "samtools=1.24" "gatk=4.6.2.0" "bcftools=1.24" "bedtools=2.31.1" "mosdepth=0.3.14" )

# Roots
CLIN_ROOT="${CLIN_ROOT:-/path/to/your/work_dir}"   # output root; run directories are created below it
# Reference FASTA; bwa-mem2 index, .fai and GATK .dict are built next to it when missing.
REF_FA="${REF_FA:-/path/to/your/reference/GCF_000149245.1_CNA3_genomic.fna}"   # path to the reference FASTA
# Raw reads; fastq_R1 / fastq_R2 in samples.tsv are relative to this directory.
RAW_ROOT="${RAW_ROOT:-/path/to/your/raw_fastq}"   # path to the raw FASTQ directory

# Pipeline dir, sample manifest, run tag
PIPELINE_DIR="${PIPELINE_DIR:-$(pwd)}"
# data/samples.tsv columns: sample, group, role, fastq_R1, fastq_R2. Samples come from the manifest, never from scanning RAW_ROOT.
SAMPLES_TSV="${SAMPLES_TSV:-${PIPELINE_DIR}/../data/samples.tsv}"

# Priority: exported RUN_TAG env  >  run_tag.txt in pipeline dir  >  today
if [[ -z "${RUN_TAG:-}" ]]; then
  if [[ -f "${PIPELINE_DIR}/run_tag.txt" ]]; then
    RUN_TAG="$(tr -d '[:space:]' < "${PIPELINE_DIR}/run_tag.txt")"
  else
    RUN_TAG="$(date +%y%m%d)"
  fi
fi

# Derived run directories (all keyed on RUN_TAG)
RUN_ROOT="${CLIN_ROOT}/g0_${RUN_TAG}"
RUNINFO_DIR="$RUN_ROOT/run_info"     # per-stage tool versions + processed sample lists
REF_DIR="$RUN_ROOT/00_ref"
FASTP_OUT="$RUN_ROOT/01_fastp"
MAP_OUT="$RUN_ROOT/02_map"
CALL_OUT="$RUN_ROOT/03_callable"
GVCF_OUT="$RUN_ROOT/04_gvcf"
JOINT_OUT="$RUN_ROOT/05_joint"
FILT_OUT="$RUN_ROOT/06_filter"
QC_OUT="$RUN_ROOT/07_qc_gate"
X1_OUT="$RUN_ROOT/X1_case_h99_check"

# Reference-level files written by 00_ref_prep
REF_DICT="${REF_FA%.*}.dict"                        # GATK looks for <ref minus extension>.dict
MASK_BED="$REF_DIR/mask_cna3_v2.bed"                # static mask
MASK_SORTED_BED="$REF_DIR/mask_cna3_v2.sorted.bed"  # same intervals, lexicographic (bedtools)
GENOME_TXT="$REF_DIR/genome.txt"                    # contig<TAB>length, .fai order
GENOME_SORTED_TXT="$REF_DIR/genome.sorted.txt"      # same, lexicographic (bedtools -g)
NUCLEAR_BED="$REF_DIR/nuclear.bed"                  # 14 nuclear chromosomes
ANALYZABLE_BED="$REF_DIR/analyzable_static.bed"     # nuclear minus static mask
INTERVAL_DIR="$REF_DIR/intervals"                   # <contig>.intervals (GATK -L)
CONTIG_LIST="$INTERVAL_DIR/chroms.list"             # contigs, .fai order

# Mitochondrial contig of the CNA3 reference; the other contigs are the nuclear genome (14 chromosomes).
MITO_CONTIG="NC_018792.1"

# Frozen parameters: plain assignments, so an exported variable of the same name cannot override them.
# Site-level hard filters on SNPs: seven independent named filters (06_filter), GATK's recommended SNP thresholds.
FLT_QUAL_MIN=30                 # QUAL >= 30 (same as GenotypeGVCFs' default -stand-call-conf)
FLT_MQ_MIN=40                   # INFO/MQ (RMS mapping quality) >= 40
FLT_QD_MIN=2                    # QD >= 2
FLT_FS_MAX=60                   # FS <= 60
FLT_MQRANKSUM_MIN=-12.5         # MQRankSum >= -12.5
FLT_READPOSRANKSUM_MIN=-8       # ReadPosRankSum >= -8
FLT_SOR_MAX=3                   # SOR <= 3
# HaplotypeCaller PCR indel model: NONE for PCR-free libraries; changing it changes the gVCFs (rerun 04-07 and X1).
HC_PCR_INDEL_MODEL=NONE
# Genotype-level N (06_filter -> gt_mask.py).
GT_GQ_MIN=30                    # GQ < 30 -> N
GT_DP_MIN=10                    # DP < 10 -> N
GT_AD_FRAC_MIN=0.90             # major-allele AD fraction < 0.90 -> N
# Callable mask = base 2 (03_callable).
CALLABLE_DP_MIN=10              # depth >= 10 ...
CALLABLE_MAPQ_MIN=30            # ... counting only reads with MAPQ >= 30
# Empirical inaccessible region: positions callable in < EMPIRICAL_MIN_CALLABLE_PCT % of the representatives.
EMPIRICAL_MIN_CALLABLE_PCT=90
EMPIRICAL_ROLE="representative"
# Strain QC gate.
QC_CALLABLE_RATIO_MIN=0.90      # callable fraction >= 0.90 x cohort median
QC_DEPTH_MEDIAN_MIN=30          # median depth >= 30x
# Intermediate allele fraction is reported only (no threshold): sites whose major-allele AD fraction lies in [QC_IAF_LO, QC_IAF_HI].
QC_IAF_LO=0.2
QC_IAF_HI=0.8

# Implementation choices (hard-coded like the frozen block).
# Major allele = the called (GT) allele, share = AD[GT] / sum(AD); '*' calls and missing AD become N.
# Site filters treat a missing annotation as not filtered.
QC_REF_ROLE="representative"    # cohort medians for the gate = the representatives
QC_IAF_MIN_AD_SUM=10            # denominator: base-1 sites with sum(AD) >= 10 and AD confirmed against the gVCF (FORMAT/ADV = M)
QC_REFDIST_ROBUST_Z_MAX=3.5     # distance-to-H99 outlier: |robust z| > 3.5 on D_norm
MOSDEPTH_EXCLUDE_FLAGS=1796     # unmapped|secondary|qcfail|duplicate = HC's default read filters
GDB_BATCH_SIZE=50               # GenomicsDBImport --batch-size
# Rule A: reference-block REF calls with GQ >= 30 and DP = block MIN_DP are judged by base 2 instead of DP (SRC=BC); see gt_mask.py.

# Threads and memory (env-overridable; they never change results)
# Budgets are per node (80 CPUs, ~1.4 TB); a stage runs allocated CPUs / cores-per-process jobs in FIFO parallel.
# 01 fastp: -w FASTP_THREADS workers plus FASTP_IO_CORES reader/writer cores per job.
FASTP_THREADS="${FASTP_THREADS:-8}"
FASTP_IO_CORES="${FASTP_IO_CORES:-2}"
# 02 mapping: a job counts as MAP_JOB_CORES (bwa-mem2 -t MAP_JOB_CORES-1 plus 1 core for fixmate/sort); samtools uses -@ MAP_POST_THREADS.
MAP_SHARDS="${MAP_SHARDS:-2}"                 # submit_all.sh: --array=0-$((MAP_SHARDS-1)), one node per shard
MAP_JOB_CORES="${MAP_JOB_CORES:-8}"
MAP_POST_THREADS="${MAP_POST_THREADS:-4}"
# samtools sort -m is per thread (8G x 4 = 32 GB holds a whole sample, no temp files).
MAP_SORT_MEM="${MAP_SORT_MEM:-8G}"
# 03 mosdepth: -t only adds BAM decompression threads (no gain beyond ~4).
MOS_THREADS="${MOS_THREADS:-3}"
MOS_JOB_CORES="${MOS_JOB_CORES:-3}"
# 04 HaplotypeCaller: single-threaded except the pair-HMM; with 2 pair-HMM threads a sample keeps 2 cores busy.
HC_SHARDS="${HC_SHARDS:-2}"
HC_PAIRHMM_THREADS="${HC_PAIRHMM_THREADS:-2}"
HC_JAVA_XMX="${HC_JAVA_XMX:-16g}"
# 05 joint genotyping: GenomicsDBImport and GenotypeGVCFs are single-threaded, so each (callset, contig) task counts as 1 core.
JOINT_READER_THREADS="${JOINT_READER_THREADS:-4}"
JOINT_JAVA_XMX="${JOINT_JAVA_XMX:-32g}"
# 06 filtering: one single-threaded gt_source.py per sample in FIFO parallel; bcftools --threads only compresses -Oz output.
FILTER_JAVA_XMX="${FILTER_JAVA_XMX:-32g}"
BCF_THREADS="${BCF_THREADS:-4}"               # bcftools --threads for -Oz output (05 concat, 06 base 1) and index
# 07 QC gate: per-sample bedtools intersects run in FIFO parallel; other steps are single-threaded.
# All GATK calls: few GC threads per JVM so that many JVMs share a node.
GATK_GC_THREADS="${GATK_GC_THREADS:-2}"

# FIFO width = allocated CPUs / cores per process, between 1 and the number of tasks; an exported override wins.
fifo_width() {   # <override-or-empty> <cpus> <cores_per_process> <tasks>
  local override="$1" cpus="$2" per="$3" tasks="$4" n
  if [[ -n "$override" ]]; then n="$override"; else n=$(( cpus / per )); fi
  (( n < 1 )) && n=1
  (( n > tasks )) && n="$tasks"
  echo "$n"
}

# Callsets: main / h99 / case are joint-called separately and never mixed.
CALLSETS=(main h99 case)

callset_of_role() {
  case "$1" in
    representative|within_patient_extra|culture_replicate_of_Cn169) echo "main" ;;
    technical_replicate)    echo "h99" ;;
    case_report_ST2_to_ST5) echo "case" ;;
    *) return 1 ;;
  esac
}

callset_prefix() {
  case "$1" in
    main) echo "st5_joint" ;;    # -> st5_joint.cna3.snp.filtered.vcf.gz
    h99)  echo "h99_joint" ;;
    case) echo "case_joint" ;;
    *) return 1 ;;
  esac
}

# Helpers

# Activate a conda env while preserving the caller's 'set -u' state.
activate_conda() {
  local env_name="$1"
  local had_u=0
  case "$-" in *u*) had_u=1 ;; esac
  local conda_base
  conda_base="$("$CONDA_LAUNCHER" info --base)"
  # shellcheck source=/dev/null
  source "${conda_base}/etc/profile.d/conda.sh"
  set +u
  conda activate "$env_name"
  (( had_u )) && set -u
  echo "[$(date '+%F %T')] conda env = ${env_name} (CONDA_PREFIX=${CONDA_PREFIX:-NA})"
}

# Apply the standard strict mode + single-thread BLAS used by all stages.
apply_strict_mode() {
  set -euo pipefail
  trap 'code=$?; echo -e "\n[ERROR] Exit $code at line $LINENO" >&2' ERR
  IFS=$'\n\t'
  shopt -s nullglob
  export PYTHONUNBUFFERED=1
  export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
}

# Print the resolved configuration banner.
print_config_banner() {
  local disk
  disk="$(df -hP "$CLIN_ROOT" 2>/dev/null | awk 'NR==2{print $4" free of "$2}' || true)"
  echo "[$(date '+%F %T')] >>> RUN_TAG      = $RUN_TAG"
  echo "[$(date '+%F %T')] >>> PIPELINE_DIR = $PIPELINE_DIR"
  echo "[$(date '+%F %T')] >>> RUN_ROOT     = $RUN_ROOT"
  echo "[$(date '+%F %T')] >>> REF_FA       = $REF_FA"
  echo "[$(date '+%F %T')] >>> RAW_ROOT     = $RAW_ROOT"
  echo "[$(date '+%F %T')] >>> SAMPLES_TSV  = $SAMPLES_TSV"
  echo "[$(date '+%F %T')] >>> CONDA_ENV    = $CONDA_ENV_MAIN"
  echo "[$(date '+%F %T')] >>> MASK_BED     = $MASK_BED"
  echo "[$(date '+%F %T')] >>> HOST / JOB   = $(hostname) / ${SLURM_JOB_ID:-NA}${SLURM_ARRAY_TASK_ID:+ (array task $SLURM_ARRAY_TASK_ID)}  CPUs=${SLURM_CPUS_PER_TASK:-NA}"
  echo "[$(date '+%F %T')] >>> DISK         = ${disk:-unknown} ($CLIN_ROOT)"
}

# Fail (return 1) if any of the given commands is missing from PATH.
require_tools() {
  local t missing=0
  for t in "$@"; do
    if ! command -v "$t" >/dev/null 2>&1; then
      echo "[ERROR] $t not in PATH (CONDA_PREFIX=${CONDA_PREFIX:-none})" >&2
      missing=1
    fi
  done
  (( missing == 0 ))
}

# Record tool versions to $RUNINFO_DIR/<stage>.versions.txt; WARN on mismatch with EXPECT_VERSIONS (never fails).
log_versions() {
  local stage="$1"; shift
  local out="$RUNINFO_DIR/${stage}.versions.txt" t v e
  mkdir -p "$RUNINFO_DIR"
  {
    echo "# stage=${stage} run_tag=${RUN_TAG} host=$(hostname) date=$(date '+%F %T')"
    echo "# conda_prefix=${CONDA_PREFIX:-NA}"
    for t in "$@"; do
      case "$t" in
        fastp)       v="$(fastp --version 2>&1 | head -n 1 || true)" ;;
        bwa-mem2)    v="$(bwa-mem2 version 2>/dev/null | tail -n 1 || true)" ;;
        gatk)        v="$(gatk --version 2>/dev/null | grep -m 1 -i 'gatk' || true)" ;;
        python3)     v="$(python3 --version 2>&1 | head -n 1 || true)" ;;
        samtools|bcftools|bedtools|mosdepth|bgzip|tabix) v="$("$t" --version 2>&1 | head -n 1 || true)" ;;
        *)           v="$(command -v "$t" || true)" ;;
      esac
      printf '%s\t%s\n' "$t" "${v:-NA}"
      for e in "${EXPECT_VERSIONS[@]}"; do
        if [[ "${e%%=*}" == "$t" && "$v" != *"${e#*=}"* ]]; then
          echo "[WARN] $t version '${v:-NA}' does not contain expected '${e#*=}'" >&2
        fi
      done
    done
  } > "$out"
  echo "[$(date '+%F %T')] >>> tool versions -> $out"
  sed 's/^/    /' "$out"
}

# Load SAMPLES_TSV into index-aligned arrays S_NAME S_GROUP S_ROLE S_R1 S_R2 S_CALLSET (file order); returns 1 on any problem.
load_samples() {
  local f="$SAMPLES_TSV"
  [[ -s "$f" ]] || { echo "[ERROR] SAMPLES_TSV missing or empty: $f" >&2; return 1; }
  S_NAME=(); S_GROUP=(); S_ROLE=(); S_R1=(); S_R2=(); S_CALLSET=()
  local -A seen=() n_role=() n_cs=()
  local lineno=0 errs=0 name group role r1 r2 extra cs r n total
  while IFS=$'\t' read -r name group role r1 r2 extra || [[ -n "${name:-}" ]]; do
    lineno=$((lineno + 1))
    r2="${r2%$'\r'}"; extra="${extra%$'\r'}"
    if (( lineno == 1 )); then
      if [[ "${name}|${group}|${role}|${r1}|${r2}" != "sample|group|role|fastq_R1|fastq_R2" || -n "$extra" ]]; then
        echo "[ERROR] $f: header must be sample<TAB>group<TAB>role<TAB>fastq_R1<TAB>fastq_R2" >&2
        return 1
      fi
      continue
    fi
    if [[ -z "${name}${group}${role}${r1}${r2}" ]]; then continue; fi
    if [[ -n "$extra" || -z "$r2" ]]; then
      echo "[ERROR] $f line $lineno: expected exactly 5 tab-separated columns" >&2
      errs=$((errs + 1)); continue
    fi
    if [[ ! "$name" =~ ^[A-Za-z0-9._-]+$ ]]; then
      echo "[ERROR] $f line $lineno: unsafe sample name '$name'" >&2
      errs=$((errs + 1)); continue
    fi
    if [[ -n "${seen[$name]:-}" ]]; then
      echo "[ERROR] $f line $lineno: duplicate sample '$name'" >&2
      errs=$((errs + 1)); continue
    fi
    seen[$name]=1
    if ! cs="$(callset_of_role "$role")"; then
      echo "[ERROR] $f line $lineno: unknown role '$role' ($name)" >&2
      errs=$((errs + 1)); continue
    fi
    S_NAME+=( "$name" ); S_GROUP+=( "$group" ); S_ROLE+=( "$role" )
    S_R1+=( "$r1" ); S_R2+=( "$r2" ); S_CALLSET+=( "$cs" )
    n_role[$role]=$(( ${n_role[$role]:-0} + 1 ))
    n_cs[$cs]=$(( ${n_cs[$cs]:-0} + 1 ))
  done < "$f"

  total=${#S_NAME[@]}
  if (( errs > 0 )); then
    echo "[ERROR] $f failed validation ($errs problem(s))" >&2
    return 1
  fi
  echo "[$(date '+%F %T')] >>> samples.tsv OK: ${total} = main ${n_cs[main]:-0} + h99 ${n_cs[h99]:-0} + case ${n_cs[case]:-0}"
}

# Resolve this task's shard of an array job -> SHARD_ID / N_SHARDS (not an array job: one shard = all samples).
resolve_shard() {
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    SHARD_ID="$SLURM_ARRAY_TASK_ID"
    N_SHARDS="${SLURM_ARRAY_TASK_COUNT:-$(( ${SLURM_ARRAY_TASK_MAX:-0} + 1 ))}"
    if [[ "${SLURM_ARRAY_TASK_MIN:-0}" != "0" ]] || (( ${SLURM_ARRAY_TASK_MAX:-$((N_SHARDS - 1))} != N_SHARDS - 1 )); then
      echo "[ERROR] submit arrays as --array=0-<K-1> (min=${SLURM_ARRAY_TASK_MIN:-?} max=${SLURM_ARRAY_TASK_MAX:-?} count=${N_SHARDS})" >&2
      return 1
    fi
  else
    SHARD_ID=0
    N_SHARDS=1
  fi
  echo "[$(date '+%F %T')] >>> shard ${SHARD_ID} of ${N_SHARDS} (samples with row index % ${N_SHARDS} == ${SHARD_ID})"
}

# Print the S_* indices handled by this shard (call load_samples + resolve_shard first).
shard_indices() {
  local i
  for i in "${!S_NAME[@]}"; do
    if (( i % N_SHARDS == SHARD_ID )); then echo "$i"; fi
  done
}

# Rebuild <out_file> (FAILED_SAMPLES.txt) from the marker files in <fail_dir> under flock (array shards share it); remove it when none.
write_failed_list() {
  local fail_dir="$1" out_file="$2"
  (
    flock -x 9
    local f names=()
    for f in "$fail_dir"/*; do
      if [[ -f "$f" ]]; then names+=( "$(basename "$f")" ); fi
    done
    if (( ${#names[@]} > 0 )); then
      printf '%s\n' "${names[@]}" | LC_ALL=C sort > "$out_file"
    else
      rm -f "$out_file"
    fi
  ) 9>>"${fail_dir}/../.failed.lock"
}

# Write a threshold as a JEXL floating-point literal (GATK compares an integer literal as a long and fails on values such as MQ=60.00).
jexl_num() {
  if [[ "$1" == *.* ]]; then echo "$1"; else echo "$1.0"; fi
}

# VariantFiltration arguments for the seven site filters, one per line (used by 06 and X1); each filter is its own JEXL expression.
vf_filter_args() {
  printf '%s\n' --missing-values-evaluate-as-failing false \
    --filter-name "QUAL${FLT_QUAL_MIN}" --filter-expression "QUAL < $(jexl_num "$FLT_QUAL_MIN")" \
    --filter-name "MQ${FLT_MQ_MIN}" --filter-expression "MQ < $(jexl_num "$FLT_MQ_MIN")" \
    --filter-name "QD${FLT_QD_MIN}" --filter-expression "QD < $(jexl_num "$FLT_QD_MIN")" \
    --filter-name "FS${FLT_FS_MAX}" --filter-expression "FS > $(jexl_num "$FLT_FS_MAX")" \
    --filter-name "MQRankSum${FLT_MQRANKSUM_MIN}" --filter-expression "MQRankSum < $(jexl_num "$FLT_MQRANKSUM_MIN")" \
    --filter-name "ReadPosRankSum${FLT_READPOSRANKSUM_MIN}" --filter-expression "ReadPosRankSum < $(jexl_num "$FLT_READPOSRANKSUM_MIN")" \
    --filter-name "SOR${FLT_SOR_MAX}" --filter-expression "SOR > $(jexl_num "$FLT_SOR_MAX")"
}

# The filter names of vf_filter_args, one per line (06 per-filter counts, submit_all.sh base-1 header check).
filter_names() {
  vf_filter_args | awk 'prev == "--filter-name" {print} {prev = $0}'
}

# Java options for GATK: fixed heap, few GC threads (many JVMs per node), tmp dir.
gatk_java_opts() {
  local xmx="$1" tmp="$2"
  echo "-Xmx${xmx} -XX:ParallelGCThreads=${GATK_GC_THREADS} -XX:ConcGCThreads=1 -Djava.io.tmpdir=${tmp}"
}
