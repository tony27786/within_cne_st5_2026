#!/bin/bash
#SBATCH -J g1_07_archive
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH -o 07_archive.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

mkdir -p "$ARCHIVE_OUT" "$RUNINFO_DIR"
rm -f "$ARCHIVE_OUT/$STAGE_DONE"
cd "$ARCHIVE_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_G1"
print_config_banner
require_tools python3 || exit 1
log_versions 07_archive python3
for d in "$PREP_OUT" "$L0_OUT" "$L1_OUT" "$TREE_OUT" "$SENS_OUT" "$DISPLAY_OUT"; do
  [[ -e "$d/$STAGE_DONE" ]] || { echo "[ERROR] upstream stage not complete: $d" >&2; exit 1; }
done

NEED=(
  "$PIPELINE_DIR/00_config.sh" "$SAMPLES_TSV" "$METADATA_TSV"
  "$PREP_OUT/r/R.tsv"
  "$L0_OUT/main/pairwise_L0.tsv" "$L0_OUT/main/n_median.tsv" "$L0_OUT/main/anchors.tsv"
  "$L0_OUT/main/ruler_A_pairs.tsv" "$L0_OUT/main/ruler_B_pairs.tsv" "$L0_OUT/main/ruler_C_summary.tsv"
  "$L0_OUT/main/snp_address.tsv" "$L0_OUT/main/snp_address_levels.tsv" "$L0_OUT/main/snp_address_clusters.tsv"
  "$L0_OUT/main/snp_address_dnorm.tsv" "$L0_OUT/sensR/pairwise_L0.tsv" "$L0_OUT/sensR/snp_address.tsv"
  "$L1_OUT/L1.aln.fasta" "$L1_OUT/L1.fconst.txt" "$L1_OUT/L1.columns.tsv" "$L1_OUT/L1.summary.tsv"
  "$L1_OUT/pca_scores.tsv" "$L1_OUT/pca_variance.tsv" "$L1_OUT/K.tsv" "$L1_OUT/pc_vs_usable_rate.tsv"
  "$TREE_OUT/iqtree_threads.txt" "$TREE_OUT/ml/L1.treefile" "$TREE_OUT/ml/L1.iqtree" "$TREE_OUT/asr/L1.treefile" "$TREE_OUT/asr/L1.aln_sites.state"
  "$TREE_OUT/branches.tsv" "$TREE_OUT/topology_parameters.tsv" "$TREE_OUT/supported_clades.tsv"
  "$SENS_OUT/sens_status.tsv" "$DISPLAY_OUT/tree/display.treefile"
)
if [[ "$(awk -F'\t' '$1 == "status" {print $2}' "$SENS_OUT/sens_status.tsv")" == rerun ]]; then
  python3 "$PIPELINE_DIR/sens_compare.py" --main-tree "$TREE_OUT/ml/L1.treefile" --sens-tree "$SENS_OUT/ml/L1R.treefile" \
    --main-pca "$L1_OUT/pca_scores.tsv" --sens-pca "$SENS_OUT/pca/pca_scores.tsv" --main-kin "$L1_OUT/K.tsv" --sens-kin "$SENS_OUT/pca/K.tsv" \
    --alrt-min "$SUPPORT_ALRT_MIN" --ufboot-min "$SUPPORT_UFBOOT_MIN" --out "$ARCHIVE_OUT/sens_R_compare.tsv" || exit 1
  NEED+=( "$ARCHIVE_OUT/sens_R_compare.tsv" "$SENS_OUT/ml/L1R.treefile" "$SENS_OUT/L1R.summary.tsv" )
fi
require_files "${NEED[@]}" || exit 1

python3 "$PIPELINE_DIR/adoption_archive.py" --clusters "$L0_OUT/main/snp_address_clusters.tsv" --supported "$TREE_OUT/supported_clades.tsv" \
  --sample-table "$PREP_OUT/samples_callset.tsv" --metadata "$METADATA_TSV" --usable-stats "$PREP_OUT/sets/main/stats.tsv" \
  --analyzable-bp "$(bed_bp "$ANALYZABLE_UNION_BED")" --ladder "$LADDER" --out-dir "$ARCHIVE_OUT" || exit 1

echo "=============================================================="
echo "[$(date '+%F %T')] >>> adoption table: $ARCHIVE_OUT/adoption_candidates.tsv (fill the user_decision column)"
echo "=============================================================="
mark_stage_done "$ARCHIVE_OUT"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
