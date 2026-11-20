#!/usr/bin/env bash
# Shared configuration for G1/G2; sourced by every stage and submit_all.sh (defines variables and helpers only, no side effects).

# SLURM
SBATCH_PARTITION="<your_partition>"   # your cluster partition; keep in sync with '#SBATCH -p' in each stage
# Heavy stages take a whole node (--cpus-per-task=80); stages whose tools are single-threaded stay small (03, 07).

# Conda: launcher and env prefix
CONDA_LAUNCHER="${CONDA_LAUNCHER:-/path/to/your/conda}"   # path to your conda launcher
CONDA_ENV_G1="${CONDA_ENV_G1:-/path/to/your/envs/g1}"   # env with bcftools bedtools samtools snp-dists iqtree3 python3 (numpy scipy networkx dendropy)

# Tested tool versions; a mismatch only warns (log_versions).
EXPECT_VERSIONS=( "bcftools=1.24" "bedtools=2.31.1" "samtools=1.24" "snp-dists=1.2.0" "iqtree3=3.1.4" "python3=3.12" )

# Roots
CLIN_ROOT="${CLIN_ROOT:-/path/to/your/work_dir}"   # output root; the G0 run directory lives below it
REF_FA="${REF_FA:-/path/to/your/reference/GCF_000149245.1_CNA3_genomic.fna}"   # path to the reference FASTA

# G0 input run (frozen, read-only): the run directory is ${CLIN_ROOT}/g0_<tag>; set G0_RUN_TAG to the tag of your finished G0 run
G0_RUN_TAG="${G0_RUN_TAG:-<g0_run_tag>}"
G0_ROOT="${CLIN_ROOT}/g0_${G0_RUN_TAG}"
G0_REF_DIR="$G0_ROOT/00_ref"
G0_MAP="$G0_ROOT/02_map"
G0_CALL="$G0_ROOT/03_callable"
G0_GVCF="$G0_ROOT/04_gvcf"
G0_JOINT="$G0_ROOT/05_joint"
G0_FILT="$G0_ROOT/06_filter"
G0_QC="$G0_ROOT/07_qc_gate"
GENOME_TXT="$G0_REF_DIR/genome.txt"                 # contig<TAB>length, .fai order
GENOME_SORTED_TXT="$G0_REF_DIR/genome.sorted.txt"   # same, lexicographic (bedtools -g)
ANALYZABLE_STATIC_BED="$G0_REF_DIR/analyzable_static.bed"   # nuclear minus static mask
ANALYZABLE_UNION_BED="$G0_QC/analyzable_union.bed"          # nuclear minus (static U empirical) = the analyzable region
CONTIG_LIST="$G0_REF_DIR/intervals/chroms.list"     # contigs, .fai order
QC_GATE_TSV="$G0_QC/qc_gate.tsv"
MITO_CONTIG="NC_018792.1"

# Pipeline dir, manifests, run tag
PIPELINE_DIR="${PIPELINE_DIR:-$(pwd)}"
# Same manifest as G0: sample, group, role, fastq_R1, fastq_R2.
SAMPLES_TSV="${SAMPLES_TSV:-${PIPELINE_DIR}/../data/samples.tsv}"
# Strain metadata TSV: sample, patient_id, admission_no, sampling_date, specimen, hiv_status, st_mlst (no outcome column).
METADATA_TSV="${PIPELINE_DIR}/../data/strain_metadata.tsv"

# Priority: exported RUN_TAG env  >  run_tag.txt in pipeline dir  >  today
if [[ -z "${RUN_TAG:-}" ]]; then
  if [[ -f "${PIPELINE_DIR}/run_tag.txt" ]]; then
    RUN_TAG="$(tr -d '[:space:]' < "${PIPELINE_DIR}/run_tag.txt")"
  else
    RUN_TAG="$(date +%y%m%d)"
  fi
fi

# Derived run directories (all keyed on RUN_TAG)
RUN_ROOT="${CLIN_ROOT}/g1_${RUN_TAG}"
RUNINFO_DIR="$RUN_ROOT/run_info"
PREP_OUT="$RUN_ROOT/01_prep"
L0_OUT="$RUN_ROOT/02_l0"
L1_OUT="$RUN_ROOT/03_l1"
TREE_OUT="$RUN_ROOT/04_tree"
SENS_OUT="$RUN_ROOT/05_sens_l1"
DISPLAY_OUT="$RUN_ROOT/06_display"
ARCHIVE_OUT="$RUN_ROOT/07_archive"
STAGE_DONE=".stage_complete"   # written last by each stage on success; submit_all.sh --from checks it upstream

# Frozen parameters: plain assignments, so an exported variable of the same name is overwritten.
# R list: positions where >= R_MIN_STRAINS strains have a measured (ADV=M) intermediate AD fraction in [R_IAF_LO, R_IAF_HI].
R_IAF_LO=0.2
R_IAF_HI=0.8
R_MIN_AD_SUM=10
R_MIN_STRAINS=10
# L0 pairwise distances and SNP address.
PAIR_GATE_FRAC=0.85             # N_shared >= 0.85 x N_median, N_median = median over the representative pairs of the same analysis
LADDER="5 10 25 50 100"         # SNP address thresholds (D_raw); a D_norm version is reported as sensitivity
T_EVENT_ROUND=5                 # t_event = max D_raw of same-day different-site pairs, rounded up to a multiple of 5
T_TECH_PERCENTILE=95            # t_tech = 95th percentile of the 10 H99 D_raw, numpy.percentile(method="linear")
LOW_USABLE_FRAC=0.90            # ruler B flag: usable bp < 0.90 x median of the representatives (mirrors the G0 callable gate)
# L1 and display tree.
L1_ROLE="representative"
MAX_N_FRAC=0.10                 # a column enters L1 when at most floor(0.10 x n_taxa) taxa are N (16 of 160, 16 of 165)
# Tree.
IQTREE_MODEL="MFP"              # ModelFinder; -fconst for the ascertainment, never +ASC
IQTREE_UFBOOT=1000
IQTREE_ALRT=1000
IQTREE_SEED=20261011
SUPPORT_ALRT_MIN=80             # supported clade: SH-aLRT >= 80 and UFBoot >= 95
SUPPORT_UFBOOT_MIN=95
ASR_LOWPOST=0.95                # descriptive only: substitutions with a max posterior < 0.95 at either end
PROBE_MODEL="GTR+G4"            # substitution model of the IQ-TREE thread-count probe (04)

# Threads and memory (env-overridable; they never change results except where noted)
# 04 fixes the IQ-TREE thread count once (it can change the tree) and records it in 04_tree/iqtree_threads.txt; IQTREE_THREADS overrides.
IQTREE_THREADS="${IQTREE_THREADS:-}"
SNPDISTS_THREADS="${SNPDISTS_THREADS:-}"      # default: all allocated CPUs
BCF_THREADS="${BCF_THREADS:-4}"

# FIFO width = allocated CPUs / cores per process, between 1 and the number of tasks; an exported override wins.
fifo_width() {   # <override-or-empty> <cpus> <cores_per_process> <tasks>
  local override="$1" cpus="$2" per="$3" tasks="$4" n
  if [[ -n "$override" ]]; then n="$override"; else n=$(( cpus / per )); fi
  (( n < 1 )) && n=1
  (( n > tasks )) && n="$tasks"
  echo "$n"
}

# Callsets and roles (same rules as G0)
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
    main) echo "st5_joint" ;;
    h99)  echo "h99_joint" ;;
    case) echo "case_joint" ;;
    *) return 1 ;;
  esac
}

# G0 file paths used by G1.
base1_vcf()    { echo "$G0_FILT/$1/$(callset_prefix "$1").cna3.snp.filtered.vcf.gz"; }
sitefilt_vcf() { echo "$G0_FILT/$1/$(callset_prefix "$1").cna3.snp.sitefilt.vcf.gz"; }
indel_vcf()    { echo "$G0_FILT/$1/$(callset_prefix "$1").cna3.indel.raw.vcf.gz"; }
base2_bed()    { echo "$G0_CALL/callable/$1.cna3.bed.gz"; }

# Pseudo-genome sets
# set -> callset of its samples / VCF / analyzable region / site-exclusion region / whether R is masked as well.
PSEUDO_SETS=(main h99 sensR_main sensR_h99)
set_vcf() {
  case "$1" in
    main|sensR_main) base1_vcf main ;;
    h99|sensR_h99)   base1_vcf h99 ;;
    *) return 1 ;;
  esac
}
set_region() { echo "$ANALYZABLE_UNION_BED"; }   # analyzable region of the set (the same for every set)
set_xname() {    # site-exclusion region name (01_prep/xmask/<name>.bed)
  case "$1" in
    main|sensR_main) echo main ;;
    h99|sensR_h99)   echo h99 ;;
    *) return 1 ;;
  esac
}
set_uses_r() { case "$1" in sensR_main|sensR_h99) return 0 ;; *) return 1 ;; esac; }

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
  echo "[$(date '+%F %T')] >>> G0_ROOT      = $G0_ROOT (read-only)"
  echo "[$(date '+%F %T')] >>> REF_FA       = $REF_FA"
  echo "[$(date '+%F %T')] >>> SAMPLES_TSV  = $SAMPLES_TSV"
  echo "[$(date '+%F %T')] >>> METADATA_TSV = $METADATA_TSV"
  echo "[$(date '+%F %T')] >>> CONDA_ENV    = ${CONDA_PREFIX:-NA}"
  echo "[$(date '+%F %T')] >>> HOST / JOB   = $(hostname) / ${SLURM_JOB_ID:-NA}  CPUs=${SLURM_CPUS_PER_TASK:-NA}"
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

# Fail (return 1) if any of the given files is missing or empty.
require_files() {
  local f missing=0
  for f in "$@"; do
    if [[ ! -s "$f" ]]; then echo "[ERROR] missing or empty: $f" >&2; missing=1; fi
  done
  (( missing == 0 ))
}

# md5 of a file (Linux md5sum).
md5_of() { md5sum "$1" | awk '{print $1}'; }

# Record tool versions to $RUNINFO_DIR/<stage>.versions.txt; WARN on mismatch with EXPECT_VERSIONS (never fails).
log_versions() {
  local stage="$1"; shift
  local out="$RUNINFO_DIR/${stage}.versions.txt" t v e
  mkdir -p "$RUNINFO_DIR"
  {
    echo "# stage=${stage} run_tag=${RUN_TAG} g0_run_tag=${G0_RUN_TAG} host=$(hostname) date=$(date '+%F %T')"
    echo "# conda_prefix=${CONDA_PREFIX:-NA}"
    echo "# cpu=$(awk -F': ' '/model name/ {print $2; exit}' /proc/cpuinfo 2>/dev/null || echo NA)"
    for t in "$@"; do
      case "$t" in
        python3)   v="$(python3 --version 2>&1 | head -n 1 || true)"
                   v="$v; $(python3 -c 'import numpy, scipy, networkx, dendropy; print("numpy", numpy.__version__, "scipy", scipy.__version__, "networkx", networkx.__version__, "dendropy", dendropy.__version__)' 2>/dev/null || echo 'python modules NA')" ;;
        iqtree3)   v="$(iqtree3 --version 2>&1 | head -n 1 || true)" ;;
        snp-dists) v="$(snp-dists -v 2>&1 | head -n 1 || true)" ;;
        samtools|bcftools|bedtools|bgzip|tabix) v="$("$t" --version 2>&1 | head -n 1 || true)" ;;
        *)         v="$(command -v "$t" || true)" ;;
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

# Load SAMPLES_TSV into index-aligned arrays S_NAME S_GROUP S_ROLE S_CALLSET (file order); returns 1 on any problem.
load_samples() {
  local f="$SAMPLES_TSV"
  [[ -s "$f" ]] || { echo "[ERROR] SAMPLES_TSV missing or empty: $f" >&2; return 1; }
  S_NAME=(); S_GROUP=(); S_ROLE=(); S_CALLSET=()
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
    S_NAME+=( "$name" ); S_GROUP+=( "$group" ); S_ROLE+=( "$role" ); S_CALLSET+=( "$cs" )
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

# Samples of a pseudo-genome set, in samples.tsv order (call load_samples first).
set_samples() {
  local set="$1" i
  for i in "${!S_NAME[@]}"; do
    case "$set" in
      main|sensR_main) [[ "${S_CALLSET[$i]}" == main ]] && echo "${S_NAME[$i]}" ;;
      h99|sensR_h99)   [[ "${S_CALLSET[$i]}" == h99 ]] && echo "${S_NAME[$i]}" ;;
    esac
  done
  return 0
}

# Write the sample table used by the Python steps: sample, group, role, callset (file order).
write_sample_table() {   # <out.tsv>
  local i
  {
    printf 'sample\tgroup\trole\tcallset\n'
    for i in "${!S_NAME[@]}"; do
      printf '%s\t%s\t%s\t%s\n' "${S_NAME[$i]}" "${S_GROUP[$i]}" "${S_ROLE[$i]}" "${S_CALLSET[$i]}"
    done
  } > "$1"
}

# Rebuild <out_file> (FAILED_SAMPLES.txt) from the marker files in <fail_dir>; remove it when there are none.
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

# FIFO semaphore and progress bar used by the per-sample stages.
open_fifo() {   # <width>
  local i
  SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
  for ((i=0;i<$1;i++)); do echo >&3; done
}
close_fifo() { wait || true; exec 3>&- 3<&- || true; }
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }
PROGRESS_WIDTH=45
progress_init() {   # <total>
  PROGRESS_TOTAL="$1"
  COUNTER_FILE="$(mktemp)"; LOCK_FILE="$(mktemp)"; echo 0 > "$COUNTER_FILE"
}
update_progress() {
  local label="$1"
  {
    flock -x 9
    local done_jobs; done_jobs="$(<"$COUNTER_FILE")"; done_jobs=$((done_jobs + 1))
    echo "$done_jobs" > "$COUNTER_FILE"
    local percent=$(( done_jobs * 100 / PROGRESS_TOTAL ))
    local bar_len=$(( done_jobs * PROGRESS_WIDTH / PROGRESS_TOTAL )) bar="" i
    for ((i=0; i<bar_len; i++)); do bar+="#"; done
    for ((i=bar_len; i<PROGRESS_WIDTH; i++)); do bar+="."; done
    printf "[PROGRESS] [%s] %3d%% (%d/%d) latest: %s\n" "$bar" "$percent" "$done_jobs" "$PROGRESS_TOTAL" "$label"
  } 9>>"$LOCK_FILE"
}
progress_done() { rm -f "$COUNTER_FILE" "$LOCK_FILE"; }

# bp of a BED (plain or gz).
bed_bp() { if [[ "$1" == *.gz ]]; then gzip -dc "$1"; else cat "$1"; fi | awk '{s += $3 - $2} END{print s + 0}'; }

# 05 and 06 wait here until 04 has written 04_tree/iqtree_threads.txt; returns 1 if 04 is gone without it, or after WAIT_THREADS_MAX_S.
WAIT_THREADS_MAX_S="${WAIT_THREADS_MAX_S:-21600}"
wait_for_iqtree_threads() {
  local f="$TREE_OUT/iqtree_threads.txt" waited=0
  while [[ ! -s "$f" ]]; do
    if command -v squeue >/dev/null 2>&1 && [[ -z "$(squeue -h -u "${USER:-$(id -un)}" -n g1_04_tree -o '%i' 2>/dev/null)" ]]; then
      sleep 10
      if [[ -s "$f" ]]; then break; fi
      echo "[ERROR] 04_tree is neither queued nor running and $f does not exist (04 failed before fixing the thread count?)" >&2
      return 1
    fi
    if (( waited >= WAIT_THREADS_MAX_S )); then echo "[ERROR] no IQ-TREE thread count after ${waited} s ($f)" >&2; return 1; fi
    if (( waited % 600 == 0 )); then echo "[$(date '+%F %T')] >>> waiting for 04_tree to fix the IQ-TREE thread count ($f)"; fi
    sleep 30; waited=$((waited + 30))
  done
  echo "[$(date '+%F %T')] >>> IQ-TREE thread count fixed by 04: $(awk -F'\t' '$1 == "threads" {print $2}' "$f")"
}

# Mark a stage complete (last line of every stage on success).
mark_stage_done() { date '+%F %T' > "$1/$STAGE_DONE"; }
