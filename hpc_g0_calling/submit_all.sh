#!/usr/bin/env bash
# submit_all.sh -- submit the G0 pipeline (new run or resume) as a SLURM dependency chain.
#
#   00_ref_prep ──┐
#   01_fastp ─────┴─ 02_map[0,1] ── 04_hc_gvcf[0,1] ─┬─ 05_joint ─┬─ 06_filter ── 07_qc_gate
#                                                    └─ 03_callable┘
#
#   02 afterok:00:01   04 afterok:02   03 afterok:02,afterany:04   05 afterok:04   06 afterok:03:05   07 afterok:06
#
# Usage:
#   cd /path/to/hpc_g0_calling
#   bash submit_all.sh                          # new run, RUN_TAG = today (yymmdd)
#   bash submit_all.sh <run_tag>                # new run with a pinned RUN_TAG (refused if that run directory already has content)
#   bash submit_all.sh <run_tag> --resume       # continue an existing run: rerun only what is incomplete or out of date
#   bash submit_all.sh <run_tag> --from 06      # force 06 and everything downstream of it; upstream stages are still checked
#   add --dry-run to any of these to print the plan without submitting anything
#
# Optional: export HC_SHARDS=1 MAP_SHARDS=1 to use one node per stage instead of two.
set -euo pipefail
shopt -s nullglob

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_DIR"

usage() { sed -n '/^# Usage:/,/^#   add --dry-run/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }
die() { echo "[ERROR] $*" >&2; exit 1; }

# Arguments
MODE="new"; FROM=""; DRY_RUN=0; TAG_ARG=""
while (( $# > 0 )); do
  case "$1" in
    --resume)  [[ "$MODE" == new ]] || die "use either --resume or --from"; MODE="resume" ;;
    --from)    [[ "$MODE" == new ]] || die "use either --resume or --from"; MODE="from"; FROM="${2:-}"; shift
               [[ "$FROM" =~ ^0[0-7]$ ]] || die "--from needs a stage 00-07, got '${FROM}'" ;;
    --dry-run) DRY_RUN=1 ;;
    -h|--help) usage; exit 0 ;;
    -*)        usage >&2; die "unknown option $1" ;;
    *)         [[ -z "$TAG_ARG" ]] || die "more than one RUN_TAG given"; TAG_ARG="$1" ;;
  esac
  shift
done
if [[ "$MODE" != new && -z "$TAG_ARG" ]]; then die "--resume / --from need the RUN_TAG of the run to continue, e.g. bash submit_all.sh <run_tag> --resume"; fi
[[ -z "$TAG_ARG" || "$TAG_ARG" =~ ^[0-9A-Za-z_-]+$ ]] || die "unsafe RUN_TAG '$TAG_ARG'"

RUN_TAG="${TAG_ARG:-$(date +%y%m%d)}"
export PIPELINE_DIR RUN_TAG
source "${PIPELINE_DIR}/00_config.sh"

# Pipeline files and manifest
for f in 00_config.sh 00_ref_prep.slurm 01_fastp.slurm 02_map.slurm 03_callable.slurm 04_hc_gvcf.slurm 05_joint.slurm 06_filter.slurm 07_qc_gate.slurm gt_mask.py gt_source.py qc_gate_table.py; do
  [[ -f "$f" ]] || die "missing pipeline file: $PIPELINE_DIR/$f"
done
[[ -f "$SAMPLES_TSV" ]] || die "missing sample manifest: $SAMPLES_TSV"
load_samples
(( MAP_SHARDS >= 1 && HC_SHARDS >= 1 )) || die "MAP_SHARDS / HC_SHARDS must be >= 1"

# Nothing of this pipeline may still be queued or running
X1_ACTIVE=0
if command -v squeue >/dev/null 2>&1; then
  QUEUE="$(squeue -h -u "${USER:-$(id -un)}" -o '%i %j %T' 2>/dev/null || true)"
  ACTIVE="$(printf '%s\n' "$QUEUE" | awk '$2 ~ /^g0_0[0-7]_/' || true)"
  if [[ -n "$ACTIVE" ]]; then
    echo "[ERROR] G0 jobs are still queued or running; nothing submitted (wait, or scancel them first):" >&2
    printf '    %s\n' "$ACTIVE" >&2
    exit 1
  fi
  if printf '%s\n' "$QUEUE" | awk '$2 == "g0_X1_case" {f = 1} END {exit !f}'; then X1_ACTIVE=1; fi
fi

# Stage table
STAGES=(00 01 02 03 04 05 06 07)
declare -A SCRIPT=( [00]=00_ref_prep.slurm [01]=01_fastp.slurm [02]=02_map.slurm [03]=03_callable.slurm [04]=04_hc_gvcf.slurm [05]=05_joint.slurm [06]=06_filter.slurm [07]=07_qc_gate.slurm )
declare -A UP=( [00]="" [01]="" [02]="00 01" [03]="00 02" [04]="00 02" [05]="00 04" [06]="00 02 03 04 05" [07]="00 02 03 06" )   # stages whose products each stage reads (rerun propagation)
declare -A SCHED=( [00]="" [01]="" [02]="00 01" [03]="02" [04]="02" [05]="04" [06]="03 05" [07]="06" )                         # direct SLURM afterok edges; with rerun propagation they cover every stage in UP
declare -A OK=() WHY=() OLDEST=() NEWEST=() RUN=() REASON=() FORCED=() JID=()

roster_of() { local i; for i in "${!S_NAME[@]}"; do if [[ "${S_CALLSET[$i]}" == "$1" ]]; then echo "${S_NAME[$i]}"; fi; done | LC_ALL=C sort | tr '\n' ' '; }

# Product checks (login node: coreutils, gzip, awk only)
EOF_HEX="1f8b08040000000000ff0600424302001b0003000000000000000000"   # BGZF end-of-file block (BAM, bgzipped VCF / BED)
bgzf_ok() { [[ "$(tail -c 28 "$1" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')" == "$EOF_HEX" ]]; }
mtime() { date -r "$1" +%s; }
MTIME_TOL_S=300
vcf_header() { gzip -cd "$1" 2>/dev/null | sed '/^#CHROM/q' || true; }   # header up to and including #CHROM
P_BAD=0; P_WHY=""; P_OLD=""; P_NEW=""
p_reset() { P_BAD=0; P_WHY=""; P_OLD=""; P_NEW=""; }
p_bad() { P_BAD=$((P_BAD + 1)); if [[ -z "$P_WHY" ]]; then P_WHY="$1"; fi; }
p_file() {   # <path> [bgzf]: must exist and be non-empty (and end with a BGZF EOF block); joins the stage's mtime range
  local f="$1" t
  if [[ ! -s "$f" ]]; then p_bad "missing or empty: $f"; return 0; fi
  if [[ "${2:-}" == bgzf ]] && ! bgzf_ok "$f"; then p_bad "truncated (no BGZF end-of-file block): $f"; return 0; fi
  t="$(mtime "$f")"
  if [[ -z "$P_OLD" || "$t" -lt "$P_OLD" ]]; then P_OLD="$t"; fi
  if [[ -z "$P_NEW" || "$t" -gt "$P_NEW" ]]; then P_NEW="$t"; fi
}
p_index() {  # <data> <index>: index present and not older than its data file
  if [[ ! -s "$2" ]]; then p_bad "missing index: $2"; return 0; fi
  if [[ -s "$1" ]] && (( $(mtime "$2") + MTIME_TOL_S < $(mtime "$1") )); then p_bad "index older than its file: $2"; fi
}
p_nofail() { local f; for f in "$@"; do if [[ -s "$f" ]]; then p_bad "failure list present: $f ($(awk 'END{print NR}' "$f") entries)"; fi; done; }
p_done() { OK[$1]=$(( P_BAD == 0 ? 1 : 0 )); WHY[$1]="${P_WHY:-}"; OLDEST[$1]="$P_OLD"; NEWEST[$1]="$P_NEW"; }

check_00() {
  p_reset
  local f n
  for f in "${REF_FA}.fai" "$REF_DICT" "${REF_FA}.0123" "${REF_FA}.amb" "${REF_FA}.ann" "${REF_FA}.bwt.2bit.64" "${REF_FA}.pac"; do
    if [[ ! -s "$f" ]]; then p_bad "missing reference index: $f"; fi     # built once next to REF_FA and reused: not part of the mtime range
  done
  for f in "$MASK_BED" "$MASK_SORTED_BED" "$GENOME_TXT" "$GENOME_SORTED_TXT" "$NUCLEAR_BED" "$ANALYZABLE_BED" "$CONTIG_LIST"; do p_file "$f"; done
  p_done 00
}
check_01() {
  p_reset
  local s
  for s in "${S_NAME[@]}"; do
    p_file "$FASTP_OUT/$s/${s}_clean_R1.fq.gz"; p_file "$FASTP_OUT/$s/${s}_clean_R2.fq.gz"; p_file "$FASTP_OUT/$s/${s}.json"
  done
  p_nofail "$FASTP_OUT/FAILED_SAMPLES.txt"
  p_done 01
}
check_02() {
  p_reset
  local s
  for s in "${S_NAME[@]}"; do
    p_file "$MAP_OUT/${s}.md.bam" bgzf; p_index "$MAP_OUT/${s}.md.bam" "$MAP_OUT/${s}.md.bam.bai"
    p_file "$MAP_OUT/qc_metrics_map/flagstat/${s}.flagstat.txt"; p_file "$MAP_OUT/qc_metrics_map/stats/${s}.stats.txt"
  done
  p_nofail "$MAP_OUT/FAILED_SAMPLES.txt"
  p_done 02
}
check_03() {
  p_reset
  local s
  for s in "${S_NAME[@]}"; do
    p_file "$CALL_OUT/callable/${s}.cna3.bed.gz" bgzf; p_index "$CALL_OUT/callable/${s}.cna3.bed.gz" "$CALL_OUT/callable/${s}.cna3.bed.gz.tbi"
    p_file "$CALL_OUT/summary/${s}.callable.tsv"; p_file "$CALL_OUT/summary/${s}.depth_by_chrom.tsv"
  done
  p_nofail "$CALL_OUT/FAILED_SAMPLES.txt"
  p_done 03
}
check_04() {
  p_reset
  local s g name hdr
  for s in "${S_NAME[@]}"; do
    g="$GVCF_OUT/${s}.g.vcf.gz"
    p_file "$g" bgzf; p_index "$g" "${g}.tbi"
    if [[ -s "$g" ]]; then
      hdr="$(vcf_header "$g")"
      name="$(printf '%s\n' "$hdr" | awk -F'\t' '/^#CHROM/ {print $10}')"
      if [[ "$name" != "$s" ]]; then p_bad "gVCF sample name '${name}' != '${s}': $g"; fi
      if [[ "$hdr" != *"--pcr-indel-model ${HC_PCR_INDEL_MODEL}"* ]]; then p_bad "gVCF not made with --pcr-indel-model ${HC_PCR_INDEL_MODEL} (older 04): $g"; fi
    fi
  done
  p_nofail "$GVCF_OUT/FAILED_SAMPLES.txt"
  p_done 04
}
check_05() {
  p_reset
  local cs f got
  for cs in "${CALLSETS[@]}"; do
    f="$JOINT_OUT/$cs/$(callset_prefix "$cs").cna3.raw.vcf.gz"
    p_file "$f" bgzf; p_index "$f" "${f}.tbi"
    if [[ -s "$f" ]]; then
      got="$(vcf_header "$f" | awk -F'\t' '/^#CHROM/ {for (i = 10; i <= NF; i++) print $i}' | LC_ALL=C sort | tr '\n' ' ')"
      if [[ "$got" != "$(roster_of "$cs")" ]]; then p_bad "sample list of $f differs from the roster"; fi
    fi
  done
  p_nofail "$JOINT_OUT/FAILED_SAMPLES.txt" "$JOINT_OUT/FAILED_TASKS.txt"
  p_done 05
}
check_06() {
  p_reset
  local cs f od hdr got n
  local -a want_filters
  mapfile -t want_filters < <(filter_names)
  for cs in "${CALLSETS[@]}"; do
    od="$FILT_OUT/$cs"; f="$od/$(callset_prefix "$cs").cna3.snp.filtered.vcf.gz"
    p_file "$f" bgzf; p_index "$f" "${f}.tbi"
    p_file "$od/filter_summary.tsv"; p_file "$od/gt_mask_stats.tsv"
    if [[ -s "$f" ]]; then
      hdr="$(vcf_header "$f")"
      if [[ "$hdr" != *"##FORMAT=<ID=ADV,"* ]]; then p_bad "base 1 written by an older 06 (no FORMAT/ADV): $f"; fi
      for n in "${want_filters[@]}"; do
        if [[ "$hdr" != *"##FILTER=<ID=${n},"* ]]; then p_bad "base 1 written with other site filters (no ##FILTER ${n}): $f"; break; fi
      done
      got="$(printf '%s\n' "$hdr" | awk -F'\t' '/^#CHROM/ {for (i = 10; i <= NF; i++) print $i}' | LC_ALL=C sort | tr '\n' ' ')"
      if [[ "$got" != "$(roster_of "$cs")" ]]; then p_bad "sample list of $f differs from the roster"; fi
    fi
  done
  p_nofail "$FILT_OUT/FAILED_SAMPLES.txt"
  p_done 06
}
check_07() {
  p_reset
  local f n
  for f in "$QC_OUT/qc_gate.tsv" "$QC_OUT/qc_gate_summary.txt" "$QC_OUT/depth_by_chrom.tsv" "$QC_OUT/intermediate_af_sites_h99_case.tsv" "$QC_OUT/mask_union_summary.tsv" "$QC_OUT/mask_union.bed"; do p_file "$f"; done
  if [[ -s "$QC_OUT/qc_gate.tsv" ]]; then
    if ! head -n 1 "$QC_OUT/qc_gate.tsv" | tr '\t' '\n' | grep -qx 'intermediate_af_frac' || head -n 1 "$QC_OUT/qc_gate.tsv" | tr '\t' '\n' | grep -qx 'pass_mixture'; then
      p_bad "qc_gate.tsv written by an older 07 (column layout)"
    fi
  fi
  p_done 07
}

# Decide what runs
if [[ "$MODE" == new ]]; then
  if [[ -d "$RUN_ROOT" && -n "$(ls -A "$RUN_ROOT" 2>/dev/null)" ]]; then
    die "run directory already has content: $RUN_ROOT. Continue it with: bash submit_all.sh ${RUN_TAG} --resume (or --from <stage>), or choose a new RUN_TAG."
  fi
  for S in "${STAGES[@]}"; do RUN[$S]=1; REASON[$S]="new run"; OK[$S]=0; WHY[$S]="-"; done
else
  [[ -d "$RUN_ROOT" ]] || die "run directory not found: $RUN_ROOT"
  echo ">>> checking the products of run ${RUN_TAG} under $RUN_ROOT ..."
  for S in "${STAGES[@]}"; do "check_$S"; done
  for S in "${STAGES[@]}"; do
    FORCED[$S]=0
    if [[ "$S" == "$FROM" ]]; then FORCED[$S]=1; fi
    for U in ${UP[$S]}; do if (( FORCED[$U] )); then FORCED[$S]=1; fi; done
  done
  stale_reason() {   # <stage>: non-empty when a complete stage is older than a complete stage it reads
    local S="$1" U
    (( OK[$S] )) || return 0
    for U in ${UP[$S]}; do
      if (( OK[$U] )) && [[ -n "${NEWEST[$U]}" && -n "${OLDEST[$S]}" ]] && (( OLDEST[$S] + MTIME_TOL_S < NEWEST[$U] )); then echo "out of date: older than the products of $U"; return 0; fi
    done
  }
  decide() {   # <stage>: RUN / REASON from forcing, completeness, freshness and the stages it reads
    local S="$1" U st
    RUN[$S]=0; REASON[$S]=""
    if (( FORCED[$S] )); then RUN[$S]=1; REASON[$S]="forced (--from ${FROM})"; return 0; fi
    if (( ! OK[$S] )); then RUN[$S]=1; REASON[$S]="incomplete: ${WHY[$S]}"; return 0; fi
    st="$(stale_reason "$S")"
    if [[ -n "$st" ]]; then RUN[$S]=1; REASON[$S]="$st"; return 0; fi
    for U in ${UP[$S]}; do
      if (( RUN[$U] )); then RUN[$S]=1; REASON[$S]="stage $U reruns"; return 0; fi
    done
  }
  RUN[01]=0
  decide 00
  decide 02                                  # 01 has not been decided yet (RUN[01]=0): judge 02 first ...
  if (( FORCED[01] )); then RUN[01]=1; REASON[01]="forced (--from ${FROM})"
  elif (( RUN[02] && ! OK[01] )); then RUN[01]=1; REASON[01]="02 has to run and the cleaned reads are incomplete: ${WHY[01]}"
  else REASON[01]=$( (( OK[01] )) && echo "" || echo "not needed (02 is reused)" ); fi
  if (( RUN[01] )); then RUN[02]=1; REASON[02]="${REASON[02]:-stage 01 reruns}"; fi   # ... then make 02 follow a rerun of 01
  for S in 03 04 05 06 07; do decide "$S"; done
fi

# Plan
N_RUN=0
echo ">>> plan for RUN_TAG=${RUN_TAG} (${MODE}${FROM:+ ${FROM}})  outputs: $RUN_ROOT"
printf '    %-6s %-12s %-8s %s\n' stage products action reason
for S in "${STAGES[@]}"; do
  prod="$( [[ "$MODE" == new ]] && echo "-" || { (( OK[$S] )) && echo complete || echo incomplete; } )"
  act="$( (( RUN[$S] )) && echo RUN || echo reuse )"
  if (( RUN[$S] )); then N_RUN=$((N_RUN + 1)); fi
  printf '    %-6s %-12s %-8s %s\n' "$S" "$prod" "$act" "${REASON[$S]:-}"
done
if (( N_RUN == 0 )); then echo ">>> nothing to submit: every stage is complete and up to date."; exit 0; fi
if (( X1_ACTIVE )) && (( RUN[02] || RUN[03] || RUN[04] || RUN[05] )); then
  die "the X1 side check (g0_X1_case) is running and reads the 02-05 products this plan would rewrite; wait for it or scancel it first."
fi
if (( RUN[01] )); then
  n_missing=0
  for i in "${!S_NAME[@]}"; do
    for fq in "$RAW_ROOT/${S_R1[$i]}" "$RAW_ROOT/${S_R2[$i]}"; do
      if [[ ! -s "$fq" ]]; then echo "[ERROR] missing raw FASTQ: $fq" >&2; n_missing=$((n_missing + 1)); fi
    done
  done
  (( n_missing == 0 )) || die "${n_missing} raw FASTQ file(s) missing; nothing submitted."
  echo ">>> raw FASTQs OK: $(( 2 * ${#S_NAME[@]} )) present under $RAW_ROOT"
fi
if (( DRY_RUN )); then echo ">>> --dry-run: nothing submitted, run_tag.txt not written."; exit 0; fi

# Submit
echo "$RUN_TAG" > "${PIPELINE_DIR}/run_tag.txt"
echo ">>> RUN_TAG = $RUN_TAG  (written to run_tag.txt)"
EXPORTS="ALL,PIPELINE_DIR=${PIPELINE_DIR},RUN_TAG=${RUN_TAG}"
sb() { local j; j="$(sbatch --parsable "$@")"; echo "${j%%;*}"; }   # strip ";cluster" if present
for S in 00 01 02 04 03 05 06 07; do   # 04 before 03: 03 waits on 04 with afterany
  (( RUN[$S] )) || continue
  ok_ids=""; any_ids=""
  for U in ${SCHED[$S]}; do if (( RUN[$U] )); then ok_ids+=":${JID[$U]}"; fi; done
  if [[ "$S" == 03 ]] && (( RUN[04] )); then any_ids=":${JID[04]}"; fi
  dep=""
  if [[ -n "$ok_ids" ]]; then dep="afterok${ok_ids}"; fi
  if [[ -n "$any_ids" ]]; then dep="${dep:+${dep},}afterany${any_ids}"; fi
  args=( --export="${EXPORTS}" )
  if [[ -n "$dep" ]]; then args=( --dependency="$dep" "${args[@]}" ); fi
  case "$S" in
    02) args=( --array=0-$((MAP_SHARDS - 1)) -o 02_map_%a.out "${args[@]}" ) ;;
    04) args=( --array=0-$((HC_SHARDS - 1)) -o 04_hc_gvcf_%a.out "${args[@]}" ) ;;
  esac
  JID[$S]="$(sb "${args[@]}" "${SCRIPT[$S]}")"
  printf '>>> %-18s: %s%s\n' "${SCRIPT[$S]%.slurm}" "${JID[$S]}" "${dep:+   ($dep)}"
done

SUBMITTED=""; for S in "${STAGES[@]}"; do if (( RUN[$S] )); then SUBMITTED+=" ${JID[$S]}"; fi; done
cat <<EOF

>>> Monitor
    squeue -u \$USER
    grep PROGRESS ${PIPELINE_DIR}/<stage>.out        # progress bars

>>> If a stage fails, its dependants stay PENDING (DependencyNeverSatisfied); cancel them, fix, then resume:
    scancel${SUBMITTED}
    bash submit_all.sh ${RUN_TAG} --resume

>>> Per-stage failures
    cat ${RUN_ROOT}/0*/FAILED_SAMPLES.txt ${RUN_ROOT}/05_joint/FAILED_TASKS.txt 2>/dev/null

>>> G0 results (after 07)
    ${RUN_ROOT}/06_filter/main/st5_joint.cna3.snp.filtered.vcf.gz   (base 1; + h99/, case/)
    ${RUN_ROOT}/03_callable/callable/<sample>.cna3.bed.gz           (base 2)
    ${RUN_ROOT}/07_qc_gate/qc_gate.tsv  qc_gate_summary.txt  mask_union_summary.tsv

>>> X1 side check (independent, after 05; can run beside 06 -> 07):  sbatch X1_case_h99_check.slurm
EOF
