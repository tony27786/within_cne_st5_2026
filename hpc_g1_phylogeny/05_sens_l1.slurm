#!/bin/bash
#SBATCH -J g1_05_sens_l1
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH -o 05_sens_l1.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

mkdir -p "$SENS_OUT"/{ml,pca} "$RUNINFO_DIR"
rm -f "$SENS_OUT/$STAGE_DONE"
cd "$SENS_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_G1"
print_config_banner
require_tools iqtree3 python3 || exit 1
log_versions 05_sens_l1 iqtree3 python3
for d in "$PREP_OUT" "$L1_OUT"; do [[ -e "$d/$STAGE_DONE" ]] || { echo "[ERROR] upstream stage not complete: $d" >&2; exit 1; }; done

OVERLAP="$(awk -F'\t' '$1 == "R_in_L1_any_class" {print $2}' "$L1_OUT/L1.summary.tsv")"
[[ "$OVERLAP" =~ ^[0-9]+$ ]] || { echo "[ERROR] R overlap not found in $L1_OUT/L1.summary.tsv" >&2; exit 1; }
if (( OVERLAP == 0 )); then
  printf 'metric\tvalue\nstatus\tnot_needed\nR_in_L1_any_class\t0\nnote\texcluding R does not change the L1 input; tree, PCA and K are not rerun\n' > sens_status.tsv
  echo "[$(date '+%F %T')] >>> R overlaps no L1 column: L1 part of the sensitivity analysis not needed"
  mark_stage_done "$SENS_OUT"
  echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
  exit 0
fi
echo "[$(date '+%F %T')] >>> R overlaps ${OVERLAP} L1 columns: rebuilding L1 without R"

python3 "$PIPELINE_DIR/l1_export.py" --taxa "$L1_OUT/taxa.tsv" --region "$ANALYZABLE_UNION_BED" --genome "$GENOME_TXT" \
  --max-n-frac "$MAX_N_FRAC" --exclude-bed "$PREP_OUT/r/R.bed" --r-bed "$PREP_OUT/r/R.bed" --out-prefix "$SENS_OUT/L1R" || exit 1
python3 "$PIPELINE_DIR/pca_kin.py" --aln L1R.aln.fasta --columns L1R.columns.tsv --usable-stats "$PREP_OUT/sets/main/stats.tsv" \
  --qc-gate "$QC_GATE_TSV" --analyzable-bp "$(bed_bp "$ANALYZABLE_UNION_BED")" --out-dir "$SENS_OUT/pca" || exit 1

wait_for_iqtree_threads || exit 1
K="$(awk -F'\t' '$1 == "threads" {print $2}' "$TREE_OUT/iqtree_threads.txt")"
[[ "$K" =~ ^[0-9]+$ ]] || { echo "[ERROR] fixed thread count not found in $TREE_OUT/iqtree_threads.txt" >&2; exit 1; }
FCONST="$(tr -d '[:space:]' < L1R.fconst.txt)"
iqtree3 -s L1R.aln.fasta -fconst "$FCONST" -m "$IQTREE_MODEL" -B "$IQTREE_UFBOOT" --alrt "$IQTREE_ALRT" -T "$K" --seed "$IQTREE_SEED" \
  --prefix ml/L1R --redo > ml/L1R.stdout 2>&1 || { echo "[ERROR] ML tree (R excluded) failed" >&2; exit 1; }
printf 'metric\tvalue\nstatus\trerun\nR_in_L1_any_class\t%s\nthreads\t%s\nseed\t%s\nalignment_md5\t%s\n' \
  "$OVERLAP" "$K" "$IQTREE_SEED" "$(md5_of L1R.aln.fasta)" > sens_status.tsv
echo "[$(date '+%F %T')] >>> L1 without R: tree $SENS_OUT/ml/L1R.treefile, PCA/K $SENS_OUT/pca/ (compared with the main L1 in 07)"

mark_stage_done "$SENS_OUT"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
