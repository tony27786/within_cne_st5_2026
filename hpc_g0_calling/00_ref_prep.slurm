#!/bin/bash
#SBATCH -J g0_00_ref
#SBATCH -p <your_partition>   # your cluster partition
#SBATCH -N 1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH -o 00_ref_prep.out

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

mkdir -p "$REF_DIR" "$INTERVAL_DIR" "$RUNINFO_DIR"
cd "$REF_DIR"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MAIN"
print_config_banner

require_tools samtools bwa-mem2 gatk bedtools md5sum || exit 1
log_versions 00_ref_prep samtools bwa-mem2 gatk bedtools
[[ -f "$REF_FA" ]] || { echo "[ERROR] REF_FA not found: $REF_FA" >&2; exit 1; }

# 1) Indexes next to REF_FA (built only if missing)
if [[ ! -s "${REF_FA}.fai" ]]; then
  echo "[$(date '+%F %T')] >>> samtools faidx"
  samtools faidx "$REF_FA"
else
  echo "[$(date '+%F %T')] >>> .fai exists, skip"
fi
bwa_ok=1
for ext in 0123 amb ann bwt.2bit.64 pac; do
  if [[ ! -s "${REF_FA}.${ext}" ]]; then bwa_ok=0; fi
done
if (( bwa_ok == 0 )); then
  echo "[$(date '+%F %T')] >>> bwa-mem2 index"
  bwa-mem2 index "$REF_FA"
else
  echo "[$(date '+%F %T')] >>> bwa-mem2 index exists, skip"
fi
if [[ ! -s "$REF_DICT" ]]; then
  echo "[$(date '+%F %T')] >>> gatk CreateSequenceDictionary -> $REF_DICT"
  gatk CreateSequenceDictionary -R "$REF_FA" -O "$REF_DICT"
else
  echo "[$(date '+%F %T')] >>> .dict exists, skip"
fi

# 2) Static mask (BED, 0-based half-open, .fai order), written verbatim
# Mask: 14 centromeres (flanking-gene intergenic regions), 25 Tcn elements (+/-500 bp), 20 kb subtelomeres, mitochondrion.
MASK_TMP="${MASK_BED}.tmp"
awk 'BEGIN{OFS="\t"} NF{print $1, $2, $3}' > "$MASK_TMP" <<'EOF'
NC_026745.1	0	20000
NC_026745.1	970168	1006713
NC_026745.1	1913873	1915374
NC_026745.1	2271499	2291499
NC_026746.1	0	20000
NC_026746.1	148994	150484
NC_026746.1	266563	268054
NC_026746.1	835383	892713
NC_026746.1	1601675	1621675
NC_026747.1	0	20000
NC_026747.1	1340040	1341518
NC_026747.1	1372001	1409632
NC_026747.1	1555141	1575141
NC_026748.1	0	20000
NC_026748.1	708803	750078
NC_026748.1	1064805	1084805
NC_026749.1	0	20000
NC_026749.1	214169	215634
NC_026749.1	703118	705065
NC_026749.1	1559982	1580220
NC_026749.1	1794975	1814975
NC_026750.1	0	20000
NC_026750.1	450734	453231
NC_026750.1	782055	819722
NC_026750.1	1402463	1422463
NC_026751.1	0	20000
NC_026751.1	525713	570802
NC_026751.1	1379503	1399503
NC_026752.1	0	20000
NC_026752.1	458008	505241
NC_026752.1	1165014	1166473
NC_026752.1	1378693	1398693
NC_026753.1	0	20000
NC_026753.1	803599	840435
NC_026753.1	1166808	1186808
NC_026754.1	0	20000
NC_026754.1	199433	230149
NC_026754.1	1039964	1059964
NC_026755.1	0	20000
NC_026755.1	866081	929452
NC_026755.1	1386605	1388743
NC_026755.1	1389510	1392312
NC_026755.1	1393007	1394641
NC_026755.1	1541994	1561994
NC_026756.1	0	20000
NC_026756.1	139632	170629
NC_026756.1	231791	233269
NC_026756.1	754062	774062
NC_026757.1	0	20000
NC_026757.1	124714	126217
NC_026757.1	128015	129561
NC_026757.1	579771	628587
NC_026757.1	736744	756744
NC_026758.1	0	20000
NC_026758.1	443128	477986
NC_026758.1	922867	942867
NC_018792.1	0	24919
EOF

mask_n=$(awk 'END{print NR}' "$MASK_TMP")
mask_bp=$(awk '{s += $3 - $2} END{print s + 0}' "$MASK_TMP")
mask_md5=$(md5sum "$MASK_TMP" | awk '{print $1}')
mv -f "$MASK_TMP" "$MASK_BED"
echo "[$(date '+%F %T')] >>> static mask OK: ${mask_n} intervals / ${mask_bp} bp -> $MASK_BED"

# 3) Derived BEDs (lexicographic order for bedtools)
awk 'BEGIN{OFS="\t"} {print $1, $2}' "${REF_FA}.fai" > "$GENOME_TXT"
LC_ALL=C sort -k1,1 "$GENOME_TXT" > "$GENOME_SORTED_TXT"
awk -v m="$MITO_CONTIG" 'BEGIN{OFS="\t"} $1 != m {print $1, 0, $2}' "$GENOME_TXT" | LC_ALL=C sort -k1,1 -k2,2n > "$NUCLEAR_BED"
LC_ALL=C sort -k1,1 -k2,2n "$MASK_BED" > "$MASK_SORTED_BED"
bedtools subtract -a "$NUCLEAR_BED" -b "$MASK_SORTED_BED" | LC_ALL=C sort -k1,1 -k2,2n > "$ANALYZABLE_BED"

analyzable_bp=$(awk '{s += $3 - $2} END{print s + 0}' "$ANALYZABLE_BED")

# 4) Per-contig GATK intervals (.fai order)
rm -f "$INTERVAL_DIR"/*.intervals
awk '{print $1}' "$GENOME_TXT" > "$CONTIG_LIST"
while IFS=$'\t' read -r contig clen; do
  printf '%s:1-%s\n' "$contig" "$clen" > "$INTERVAL_DIR/${contig}.intervals"
done < "$GENOME_TXT"
INTERVAL_FILES=( "$INTERVAL_DIR"/*.intervals )
echo "[$(date '+%F %T')] >>> intervals: ${#INTERVAL_FILES[@]} contigs -> $INTERVAL_DIR"

echo "=============================================================="
echo "[$(date '+%F %T')] >>> 00_ref_prep complete"
echo "    static mask : $MASK_BED (${mask_n} intervals, ${mask_bp} bp, md5 ${mask_md5})"
echo "    analyzable  : $ANALYZABLE_BED (${analyzable_bp} bp)"
echo "    intervals   : $CONTIG_LIST"
echo "=============================================================="
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
