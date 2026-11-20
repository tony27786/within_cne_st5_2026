#!/bin/bash
#SBATCH -J g1_06_display
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH -o 06_display.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

mkdir -p "$DISPLAY_OUT/tree" "$RUNINFO_DIR"
rm -f "$DISPLAY_OUT/$STAGE_DONE"
cd "$DISPLAY_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_G1"
print_config_banner
require_tools iqtree3 python3 || exit 1
log_versions 06_display iqtree3 python3
[[ -e "$PREP_OUT/$STAGE_DONE" ]] || { echo "[ERROR] upstream stage not complete: $PREP_OUT" >&2; exit 1; }
load_samples || exit 1

: > taxa.tsv; H99=()
for i in "${!S_NAME[@]}"; do
  case "${S_ROLE[$i]}" in
    representative)      printf '%s\t%s\n' "${S_NAME[$i]}" "$PREP_OUT/sets/main/fasta/${S_NAME[$i]}.fa" >> taxa.tsv ;;
    technical_replicate) printf '%s\t%s\n' "${S_NAME[$i]}" "$PREP_OUT/sets/h99/fasta/${S_NAME[$i]}.fa" >> taxa.tsv; H99+=( "${S_NAME[$i]}" ) ;;
  esac
done

python3 "$PIPELINE_DIR/l1_export.py" --taxa taxa.tsv --region "$ANALYZABLE_UNION_BED" --genome "$GENOME_TXT" \
  --max-n-frac "$MAX_N_FRAC" --out-prefix "$DISPLAY_OUT/display" || exit 1
wait_for_iqtree_threads || exit 1
K="$(awk -F'\t' '$1 == "threads" {print $2}' "$TREE_OUT/iqtree_threads.txt")"
[[ "$K" =~ ^[0-9]+$ ]] || { echo "[ERROR] fixed thread count not found in $TREE_OUT/iqtree_threads.txt" >&2; exit 1; }
FCONST="$(tr -d '[:space:]' < display.fconst.txt)"
OUTGROUP="$(IFS=,; echo "${H99[*]}")"
iqtree3 -s display.aln.fasta -fconst "$FCONST" -m "$IQTREE_MODEL" -B "$IQTREE_UFBOOT" --alrt "$IQTREE_ALRT" -T "$K" --seed "$IQTREE_SEED" \
  -o "$OUTGROUP" --prefix tree/display --redo > tree/display.stdout 2>&1 || { echo "[ERROR] display tree failed" >&2; exit 1; }
printf 'threads\t%s\nseed\t%s\noutgroup\t%s\nalignment_md5\t%s\n' "$K" "$IQTREE_SEED" "$OUTGROUP" "$(md5_of display.aln.fasta)" > display_run.tsv

echo "[$(date '+%F %T')] >>> display tree: $DISPLAY_OUT/tree/display.treefile (rooted on ${OUTGROUP}; draw with a broken axis)"
mark_stage_done "$DISPLAY_OUT"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
