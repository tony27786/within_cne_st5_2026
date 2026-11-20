#!/usr/bin/env bash
# submit_all.sh -- submit the G1/G2 pipeline (new run, or rerun from a stage) as a SLURM dependency chain.
#
#   01_prep ─┬─ 02_l0 ──────────────────────────────┐
#            └─ 03_l1 ─┬─ 04_tree ──────────────────┼─ 07_archive
#                      ├─ 05_sens_l1 ···(threads)···┤
#                      └─ 06_display ···(threads)···┘
#
#   02, 03 afterok:01   04, 05, 06 afterok:03   07 afterok:02:04:05:06
#   05 and 06 wait inside the job until 04 has written 04_tree/iqtree_threads.txt, then build their trees beside 04's.
#
# Usage:
#   cd /path/to/hpc_g1_phylogeny
#   bash submit_all.sh                      # new run, RUN_TAG = today (yymmdd); refused if g1_<tag> already has content
#   bash submit_all.sh <run_tag>            # new run with a pinned RUN_TAG
#   bash submit_all.sh <run_tag> --from 04  # rerun 04 and everything downstream of it in an existing run; every upstream stage must be complete
#   add --dry-run to print the plan without submitting anything
set -euo pipefail
shopt -s nullglob

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_DIR"

usage() { sed -n '/^# Usage:/,/^#   add --dry-run/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }
die() { echo "[ERROR] $*" >&2; exit 1; }

FROM=""; DRY_RUN=0; TAG_ARG=""
while (( $# > 0 )); do
  case "$1" in
    --from)    FROM="${2:-}"; shift; [[ "$FROM" =~ ^0[1-7]$ ]] || die "--from needs a stage 01-07, got '${FROM}'" ;;
    --dry-run) DRY_RUN=1 ;;
    -h|--help) usage; exit 0 ;;
    -*)        usage >&2; die "unknown option $1" ;;
    *)         [[ -z "$TAG_ARG" ]] || die "more than one RUN_TAG given"; TAG_ARG="$1" ;;
  esac
  shift
done
if [[ -n "$FROM" && -z "$TAG_ARG" ]]; then die "--from needs the RUN_TAG of the run, e.g. bash submit_all.sh <run_tag> --from 04"; fi
[[ -z "$TAG_ARG" || "$TAG_ARG" =~ ^[0-9A-Za-z_-]+$ ]] || die "unsafe RUN_TAG '$TAG_ARG'"

RUN_TAG="${TAG_ARG:-$(date +%y%m%d)}"
export PIPELINE_DIR RUN_TAG
source "${PIPELINE_DIR}/00_config.sh"

# Pipeline files, manifests, G0 input
for f in 00_config.sh g1lib.py r_list.py xmask.py pseudo_check.py l0_distances.py rulers_address.py \
         l1_export.py pca_kin.py asr_branches.py state_follow.py sens_compare.py adoption_archive.py \
         01_prep.slurm 02_l0.slurm 03_l1.slurm 04_tree.slurm 05_sens_l1.slurm 06_display.slurm 07_archive.slurm; do
  [[ -f "$f" ]] || die "missing pipeline file: $PIPELINE_DIR/$f"
done
for f in "$SAMPLES_TSV" "$METADATA_TSV"; do
  [[ -f "$f" ]] || die "missing data table: $f"
done
load_samples
[[ -d "$G0_ROOT" ]] || die "G0 run not found: $G0_ROOT"
for f in "$(base1_vcf main)" "$(base1_vcf h99)" "$QC_GATE_TSV" "$ANALYZABLE_UNION_BED" "$REF_FA"; do
  [[ -s "$f" ]] || die "missing G0 input: $f"
done

# Nothing of this pipeline may still be queued or running
if command -v squeue >/dev/null 2>&1; then
  ACTIVE="$(squeue -h -u "${USER:-$(id -un)}" -o '%i %j %T' 2>/dev/null | awk '$2 ~ /^g1_0[1-7]_/' || true)"
  if [[ -n "$ACTIVE" ]]; then
    echo "[ERROR] G1 jobs are still queued or running; nothing submitted:" >&2
    printf '    %s\n' "$ACTIVE" >&2
    exit 1
  fi
fi

# Stage table
STAGES=(01 02 03 04 05 06 07)
declare -A SCRIPT=( [01]=01_prep.slurm [02]=02_l0.slurm [03]=03_l1.slurm [04]=04_tree.slurm [05]=05_sens_l1.slurm [06]=06_display.slurm [07]=07_archive.slurm )
declare -A DIR=( [01]="$PREP_OUT" [02]="$L0_OUT" [03]="$L1_OUT" [04]="$TREE_OUT" [05]="$SENS_OUT" [06]="$DISPLAY_OUT" [07]="$ARCHIVE_OUT" )
declare -A DEPS=( [01]="" [02]="01" [03]="01" [04]="03" [05]="03" [06]="03" [07]="02 04 05 06" )   # direct SLURM afterok edges
declare -A DATA=( [01]="" [02]="01" [03]="01" [04]="03" [05]="01 03 04" [06]="01 04" [07]="02 04 05 06" )   # what each stage reads (--from closure, completeness)
declare -A RUN=() JID=()

upstream_of() {   # transitive data dependencies of a stage
  local S="$1" U out=""
  for U in ${DATA[$S]}; do out+=" $U $(upstream_of "$U")"; done
  echo "$out"
}

if [[ -z "$FROM" ]]; then
  if [[ -d "$RUN_ROOT" && -n "$(ls -A "$RUN_ROOT" 2>/dev/null)" ]]; then
    die "run directory already has content: $RUN_ROOT. Rerun part of it with --from <stage>, or choose a new RUN_TAG."
  fi
  for S in "${STAGES[@]}"; do RUN[$S]=1; done
else
  [[ -d "$RUN_ROOT" ]] || die "run directory not found: $RUN_ROOT"
  for S in "${STAGES[@]}"; do
    RUN[$S]=0
    if [[ "$S" == "$FROM" ]] || [[ " $(upstream_of "$S") " == *" $FROM "* ]]; then RUN[$S]=1; fi
  done
  for S in "${STAGES[@]}"; do
    (( RUN[$S] )) || continue
    for U in $(upstream_of "$S"); do
      if (( ! RUN[$U] )) && [[ ! -e "${DIR[$U]}/$STAGE_DONE" ]]; then die "stage $U is not complete (${DIR[$U]}/$STAGE_DONE missing); rerun from it first"; fi
    done
  done
fi

echo ">>> plan for RUN_TAG=${RUN_TAG}${FROM:+ (--from ${FROM})}  outputs: $RUN_ROOT  (G0 input: $G0_ROOT)"
printf '    %-6s %-20s %s\n' stage script action
for S in "${STAGES[@]}"; do
  printf '    %-6s %-20s %s\n' "$S" "${SCRIPT[$S]}" "$( (( RUN[$S] )) && echo RUN || echo "reuse ($( [[ -e "${DIR[$S]}/$STAGE_DONE" ]] && echo complete || echo incomplete))" )"
done
if (( DRY_RUN )); then echo ">>> --dry-run: nothing submitted, run_tag.txt not written."; exit 0; fi

# Submit
echo "$RUN_TAG" > "${PIPELINE_DIR}/run_tag.txt"
mkdir -p "$RUN_ROOT"
if (( RUN[04] )); then rm -f "$TREE_OUT/iqtree_threads.txt"; fi   # 05 / 06 must wait for the thread count of this 04
EXPORTS="ALL,PIPELINE_DIR=${PIPELINE_DIR},RUN_TAG=${RUN_TAG}"
sb() { local j; j="$(sbatch --parsable "$@")"; echo "${j%%;*}"; }
for S in "${STAGES[@]}"; do
  (( RUN[$S] )) || continue
  ok_ids=""
  for U in ${DEPS[$S]}; do if (( RUN[$U] )); then ok_ids+=":${JID[$U]}"; fi; done
  args=( --export="${EXPORTS}" )
  if [[ -n "$ok_ids" ]]; then args=( --dependency="afterok${ok_ids}" "${args[@]}" ); fi
  JID[$S]="$(sb "${args[@]}" "${SCRIPT[$S]}")"
  printf '>>> %-14s: %s%s\n' "${SCRIPT[$S]%.slurm}" "${JID[$S]}" "${ok_ids:+   (afterok${ok_ids})}"
done

SUBMITTED=""; for S in "${STAGES[@]}"; do if (( RUN[$S] )); then SUBMITTED+=" ${JID[$S]}"; fi; done
cat <<EOF

>>> Monitor
    squeue -u \$USER
    grep PROGRESS ${PIPELINE_DIR}/01_prep.out

>>> If a stage fails, its dependants stay PENDING (DependencyNeverSatisfied); cancel them, fix, then rerun from it:
    scancel${SUBMITTED}
    bash submit_all.sh ${RUN_TAG} --from <stage>

>>> Results
    ${RUN_ROOT}/02_l0/main/        rulers A / B / C, anchors, SNP address   (02_l0/sensR/: R excluded)
    ${RUN_ROOT}/03_l1/             L1 alignment, PCA, K
    ${RUN_ROOT}/04_tree/           ML tree, ancestral reconstruction, topology parameters, supported clades
    ${RUN_ROOT}/07_archive/        adoption_candidates.tsv
EOF
