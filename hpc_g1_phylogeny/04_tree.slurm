#!/bin/bash
#SBATCH -J g1_04_tree
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH -o 04_tree.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
mkdir -p "$TREE_OUT"/{probe,ml,asr} "$RUNINFO_DIR"
rm -f "$TREE_OUT/$STAGE_DONE" "$TREE_OUT/iqtree_threads.txt"   # 05 / 06 wait for a fresh thread-count file
cd "$TREE_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_G1"
print_config_banner
require_tools iqtree3 python3 || exit 1
log_versions 04_tree iqtree3 python3
[[ -e "$L1_OUT/$STAGE_DONE" ]] || { echo "[ERROR] stage 03 not complete" >&2; exit 1; }
ALN="$L1_OUT/L1.aln.fasta"; COLS="$L1_OUT/L1.columns.tsv"; FCONST="$(tr -d '[:space:]' < "$L1_OUT/L1.fconst.txt")"
require_files "$ALN" "$COLS" || exit 1

# 1) fixed thread count
if [[ -n "$IQTREE_THREADS" ]]; then
  K="$IQTREE_THREADS"; SRC="override (exported IQTREE_THREADS)"
else
  OMP_NUM_THREADS=1 iqtree3 -s "$ALN" -fconst "$FCONST" -m "$PROBE_MODEL" -n 0 -T AUTO --threads-max "$THREADS_TOTAL" \
    --seed "$IQTREE_SEED" --prefix probe/probe --redo > probe/probe.stdout 2>&1 || { echo "[ERROR] thread probe failed (probe/probe.stdout)" >&2; exit 1; }
  K="$(grep -h 'BEST NUMBER OF THREADS' probe/probe.log probe/probe.stdout 2>/dev/null | head -n 1 | grep -oE '[0-9]+' | tail -n 1 || true)"
  [[ "$K" =~ ^[0-9]+$ ]] || { echo "[ERROR] no 'BEST NUMBER OF THREADS' line in probe/probe.log; export IQTREE_THREADS to set it" >&2; exit 1; }
  SRC="probe -T AUTO --threads-max ${THREADS_TOTAL} (-m ${PROBE_MODEL} -n 0)"
fi
printf 'threads\t%s\nsource\t%s\nseed\t%s\niqtree\t%s\ncpu\t%s\nalignment_md5\t%s\n' "$K" "$SRC" "$IQTREE_SEED" \
  "$(iqtree3 --version 2>&1 | head -n 1)" "$(awk -F': ' '/model name/ {print $2; exit}' /proc/cpuinfo 2>/dev/null || echo NA)" "$(md5_of "$ALN")" > iqtree_threads.txt.tmp
mv -f iqtree_threads.txt.tmp iqtree_threads.txt   # atomic: 05 and 06 start their trees as soon as this file exists
echo "[$(date '+%F %T')] >>> IQ-TREE threads fixed at ${K} (${SRC})"

# 2) ML tree
iqtree3 -s "$ALN" -fconst "$FCONST" -m "$IQTREE_MODEL" -B "$IQTREE_UFBOOT" --alrt "$IQTREE_ALRT" -T "$K" --seed "$IQTREE_SEED" \
  --prefix ml/L1 --redo > ml/L1.stdout 2>&1 || { echo "[ERROR] ML tree failed (ml/L1.log)" >&2; exit 1; }
BEST="$(awk -F': ' '/Best-fit model according to BIC/ {print $2; exit}' ml/L1.iqtree | tr -d '[:space:]')"
[[ -n "$BEST" ]] || { echo "[ERROR] best-fit model not found in ml/L1.iqtree" >&2; exit 1; }
echo "[$(date '+%F %T')] >>> ML tree done; best-fit model (BIC): ${BEST}"

# 3) ancestral reconstruction on the ML topology
python3 "$PIPELINE_DIR/asr_branches.py" strip --tree ml/L1.treefile --out asr/ml_topology.nwk || exit 1
rm -f asr/L1.state asr/L1.aln_sites.state
iqtree3 -s "$ALN" -fconst "$FCONST" -t asr/ml_topology.nwk --tree-fix -m "$BEST" --ancestral -T "$K" --seed "$IQTREE_SEED" \
  --prefix asr/L1 --redo > asr/L1.stdout 2>&1 &
IQ_PID=$!
python3 "$PIPELINE_DIR/state_follow.py" --path asr/L1.state --pid "$IQ_PID" --columns "$COLS" --out asr/L1.aln_sites.state &
FOLLOW_PID=$!
if wait "$IQ_PID"; then rc_iq=0; else rc_iq=$?; fi
if wait "$FOLLOW_PID"; then rc_follow=0; else rc_follow=$?; fi
(( rc_iq == 0 )) || { echo "[ERROR] ancestral reconstruction failed (exit ${rc_iq}; asr/L1.log)" >&2; exit 1; }
(( rc_follow == 0 )) || { echo "[ERROR] reading / reducing asr/L1.state failed (exit ${rc_follow})" >&2; exit 1; }
rm -f asr/L1.state
echo "[$(date '+%F %T')] >>> .state reduced to the alignment columns: asr/L1.aln_sites.state"

# 4) per-branch counts and topology
python3 "$PIPELINE_DIR/asr_branches.py" count --state asr/L1.aln_sites.state --asr-tree asr/L1.treefile --ml-tree ml/L1.treefile \
  --aln "$ALN" --columns "$COLS" --qc-gate "$QC_GATE_TSV" --lowpost "$ASR_LOWPOST" \
  --alrt-min "$SUPPORT_ALRT_MIN" --ufboot-min "$SUPPORT_UFBOOT_MIN" --out-dir "$TREE_OUT" || exit 1
printf 'best_model\t%s\n' "$BEST" >> iqtree_threads.txt

echo "=============================================================="
echo "[$(date '+%F %T')] >>> tree: $TREE_OUT/ml/L1.treefile (SH-aLRT/UFBoot), $TREE_OUT/asr/L1.aln_sites.state"
sed 's/^/    /' topology_parameters.tsv
echo "=============================================================="
mark_stage_done "$TREE_OUT"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
