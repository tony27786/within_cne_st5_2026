# *Cryptococcus neoformans* ST5 Genomic Analysis

This repository contains the analysis scripts and de-identified input tables to reproduce the genomic, phylogenetic, and genome–phenotype association analyses in *"Genomic diversity within and between patients among Cryptococcus neoformans sequence type 5 clinical isolates."*

## Repository Structure

```text
.
├── data/               # De-identified clinical tables and sample manifests
├── hpc_g0_calling/     # SLURM pipeline: read QC, mapping, variant calling, filtering, QC gate
├── hpc_g1_phylogeny/   # SLURM pipeline: pseudo-genomes, SNP distances, ML phylogeny, PCA/kinship
├── association/        # Genome–phenotype mixed-model and sensitivity analyses (Python, R)
├── tables/             # Scripts generating summary TSV tables (R)
└── figures/            # Scripts generating manuscript and supplementary figures (R)
```

## Requirements

| Stage | Environment | Key Tools / Packages |
|---|---|---|
| Variant calling (`hpc_g0_calling/`) | SLURM, Conda | `fastp` 1.3.7, `bwa-mem2` 2.3, `samtools` 1.24, `gatk` 4.6.2.0, `bcftools` 1.24, `bedtools` 2.31.1, `mosdepth` 0.3.14, Python 3.12 |
| Phylogeny & distances (`hpc_g1_phylogeny/`) | SLURM, Conda | `bcftools` 1.24, `bedtools` 2.31.1, `samtools` 1.24, `snp-dists` 1.2.0, `iqtree3` 3.1.4, Python 3.12 (`numpy`, `scipy`, `networkx`, `dendropy`) |
| Association (`association/`) | Local (R, Python) | Python 3 (stdlib), R 4.6.0 (`gaston`, `RLRsim`, `CompQuadForm`, `car`, `boot`) |
| Tables & figures (`tables/`, `figures/`) | Local (R) | R 4.6.0 (`ggplot2`, `patchwork`, `ggtree`, `ape`, `ggnewscale`, `ggrepel`, `scales`) |

## Usage

### 1. HPC pipelines (SLURM)

Before running, update the cluster partition (`#SBATCH -p <your_partition>`) and path placeholders (`CONDA_LAUNCHER`, `CONDA_ENV_*`, `CLIN_ROOT`, `REF_FA`, `RAW_ROOT`) in each pipeline's `00_config.sh`.

```bash
# Stage 1: reads -> filtered SNP VCF, callable masks, QC gate
cd hpc_g0_calling
bash submit_all.sh <run_tag>
# Optional side check (after stage 05): sbatch X1_case_h99_check.slurm

# Stage 2: pseudo-genomes, pairwise SNP distances, ML trees, PCA/kinship
cd ../hpc_g1_phylogeny
export G0_RUN_TAG=<g0_run_tag>
bash submit_all.sh <run_tag>
```

### 2. Link HPC outputs for local analysis

From the repository root, symlink (or copy) the completed HPC output directories under `results/`:

```bash
mkdir -p results
ln -s /path/to/g0_<g0_run_tag> results/g0
ln -s /path/to/g1_<g1_run_tag> results/g1
```

If reproducing `tables/table2.R`, generate the Cn13 allele-depth file under `results/misscheck/`:

```bash
mkdir -p results/misscheck
bcftools query -s Cn14,Cn95,Cn211,Cn13 -f '%CHROM\t%POS[\t%SAMPLE|%GT|%AD]\n' \
  results/g0/06_filter/main/st5_joint.cna3.snp.filtered.vcf.gz | \
  awk -F'\t' 'BEGIN{OFS="\t"} {for(i=3;i<=6;i++){split($i,a,"|"); gt[a[1]]=a[2]; ad[a[1]]=a[3]}
    for(k=1;k<=2;k++){c=(k==1?"Cn95":"Cn211"); if(gt["Cn14"]!="." && gt[c]!="." && gt["Cn14"]!=gt[c])
      print c,$1,$2,gt["Cn14"],gt[c],gt["Cn13"],ad["Cn13"]}}' > results/misscheck/cn13_mix_check.tsv
```

### 3. Local statistical analysis, tables, and figures

Run from the repository root:

```bash
# Genome-phenotype association
python3 association/build_design.py
Rscript association/host_prep.R        # optional diagnostics
Rscript association/a1_run.R
Rscript association/a1_treat_sens.R

# Summary tables (TSV -> tables/out/)
for s in table1 table2 table_a1; do
  Rscript tables/${s}.R
done

# Figures (PDF & PNG -> figures/out/)
for s in fig2 fig3 figS_trees; do
  Rscript figures/${s}.R
done
```

## Data Availability

- Raw sequencing reads are deposited under `<repository and accession number>`.
- De-identified clinical and isolate tables (with personal identifiers removed and admission numbers replaced by anonymous codes `ADM001`–`ADM194`) are provided in `data/`.

## **Acknowledgments**

This workflow relies on open-source bioinformatics tools (including `fastp`, `bwa-mem2`, `samtools`, `GATK`, `bcftools`, `bedtools`, `mosdepth`, `snp-dists`, and `IQ-TREE`) and the R/Python scientific ecosystems. We sincerely thank their authors and maintainers.

## Citation

`<Authors, Title, Journal, Year, DOI>`

## License

`<License>`
