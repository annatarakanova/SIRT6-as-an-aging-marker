# Differential expression analysis of SIRT6 datasets

This directory contains script for differential expression (DE) analysis of bulk RNA-seq datasets related to SIRT6.

The analysis is designed to work with the **SIRT6 database** stored on the lab server and performs DESeq2-based differential expression analysis in a reproducible way.

---

## Script overview

### `DE_analysis_script.R`

This script performs differential expression analysis of SIRT6-related bulk RNA-seq datasets for a specified organism.

Key features:

- Reads count matrices stored in **Parquet** format
- Uses sample metadata from the SIRT6 database
- Performs genotype-resolved comparisons:
  - WT vs each non-WT genotype (e.g. KO, Het, OE, mutant)
- Stratifies experiments by **treatment**, **cell type**, **tissue** or **condition** when applicable
- Applies conservative low-expression filtering
- Uses **DESeq2** for differential expression analysis
- Dynamically incorporates **biological covariates** (e.g. sex, age, strain) when available and informative
- Writes results in **Parquet** format for downstream meta-analysis

---

## Analysis logic

```text
For each organism, the script performs the following steps:
	1.	Iterate over all experiments (GSEs)
	2.	Within each experiment:
	•	Split samples into homogeneous strata based on treatment, cell type, tissue or condition
	3.	Within each stratum:
	•	Identify all genotypes present
	•	For each non-WT genotype:
		- Perform a WT vs genotype differential expression analysis
		- Apply conservative expression filtering
	•	Fit a DESeq2 model with genotype as the main effect
	•	Include additional covariates (sex, age, strain) only if informative
```
Each contrast is analyzed independently and saved as a separate result file.

---

## Requirements

### Software

- R (version ≥ 4.1 recommended)

### R packages

The following R packages are required:

- `arrow`
- `dplyr`
- `tibble`
- `DESeq2`

---

## Input data

The script expects the SIRT6 database to have the following structure:

```text
SIRT6_db/
├── expression/
│   └── <organism>/
│       └── <GSE_ID>.parquet
└── metadata/
    └── <organism>/
        ├── samples.parquet
        ├── samples_to_experiment.parquet
```

---

## Usage 

Run the script from the command line using Rscript.

### Example command

```bash
Rscript DE_analysis_script.R \
  --organism homo_sapiens \
  --expr_path /path/to/directory/with/count/matrices/for/human \
  --meta_path /path/to/directory/with/human/metadata \
  --out_dir /path/to/output/directory/for/homo_sapiens
```

---

## Output

The script creates the following directory structure in the specified output directory:

```text
out_dir/
├── results/
│   ├── <GSE_ID>/
│   │   ├── <GSE_ID>_<stratum>_<genotype>_vs_WT_deseq2.parquet
│   │   ├── <GSE_ID>_<stratum>_<genotype>_vs_WT_deseq2.parquet
│   │   └── ...
│   ├── <GSE_ID>/
│   │   └── ...
│   └── ...
└── logs/
    └── errors.log
```
### Output description

- results/

Each file is a data frame in parquet format that contains DESeq2 results (log2FC, adjusted p-value)for a single contrast (WT vs one genotype) within a specific experiment and stratum.

- logs/errors.log

Contains a record of experiments or strata that failed during processing,
along with error messages.

