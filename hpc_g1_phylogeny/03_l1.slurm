#!/bin/bash
#SBATCH -J g1_03_l1
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH -o 03_l1.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

mkdir -p "$L1_OUT" "$RUNINFO_DIR"
rm -f "$L1_OUT/$STAGE_DONE"
cd "$L1_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_G1"
print_config_banner
require_tools python3 || exit 1
log_versions 03_l1 python3
[[ -e "$PREP_OUT/$STAGE_DONE" ]] || { echo "[ERROR] stage 01 not complete" >&2; exit 1; }
load_samples || exit 1

: > taxa.tsv
for i in "${!S_NAME[@]}"; do
  if [[ "${S_ROLE[$i]}" == "$L1_ROLE" ]]; then printf '%s\t%s\n' "${S_NAME[$i]}" "$PREP_OUT/sets/main/fasta/${S_NAME[$i]}.fa" >> taxa.tsv; fi
done

python3 "$PIPELINE_DIR/l1_export.py" --taxa taxa.tsv --region "$ANALYZABLE_UNION_BED" --genome "$GENOME_TXT" \
  --max-n-frac "$MAX_N_FRAC" --r-bed "$PREP_OUT/r/R.bed" --out-prefix "$L1_OUT/L1" || exit 1
python3 "$PIPELINE_DIR/pca_kin.py" --aln L1.aln.fasta --columns L1.columns.tsv --usable-stats "$PREP_OUT/sets/main/stats.tsv" \
  --qc-gate "$QC_GATE_TSV" --analyzable-bp "$(bed_bp "$ANALYZABLE_UNION_BED")" --out-dir "$L1_OUT" || exit 1

echo "=============================================================="
echo "[$(date '+%F %T')] >>> L1: $L1_OUT/L1.aln.fasta  -fconst $(cat L1.fconst.txt)"
sed 's/^/    /' L1.summary.tsv
echo "    R overlap with L1 (any class): $(awk -F'\t' '$1 == "R_in_L1_any_class" {print $2}' L1.summary.tsv)  -> 05_sens_l1 reruns L1 only when > 0"
sed 's/^/    /' pc_vs_usable_rate.tsv
echo "=============================================================="
mark_stage_done "$L1_OUT"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
