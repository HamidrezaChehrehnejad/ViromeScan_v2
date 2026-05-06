# ViromeScan v2.2

**A multi-host, probabilistic pipeline for virome detection and profiling from shotgun metagenomic data.**

Developed by [Hamidreza Chehrehnejad](https://github.com/HamidrezaChehrehnejad) — MSc Molecular Biology and Genetics, University of Pavia (2025/2026).
Supervised by Prof. Davide Sassera (University of Pavia) and Prof. Simone Rampelli (University of Bologna).

## Overview

ViromeScan v2.2 addresses the two core limitations of conventional virome profiling tools:

1. **Multi-mapping false positives:** resolved by replacing Bowtie2 best-hit assignment with **Salmon Expectation-Maximization (EM)** probabilistic quantification.
2. **Human-centric databases:** resolved by integrating host-specific depletion and viral indices for **human, cattle, pig, and grapevine** under a One Health framework.

## Pipeline

Raw reads (FASTQ)
│
├── Step 1: Broad viral pre-screen
│ Bowtie2 (--sensitive-local --no-unal)
│ Retains only viral candidate reads
│
├── Step 2: Host genome depletion
│ Bowtie2 (--un-conc)
│ Human: GRCh38
│ Cattle: Bos taurus
│ Pig: Sus scrofa
│ Grapevine: Vitis vinifera
│
├── Step 3: Microbiome depletion
│ Bowtie2 (--un-conc)
│ Human: HRGMv2 hybrid index
│ Other hosts: host-specific MAG databases
│
├── Step 4: Fungal depletion (human only)
│ Bowtie2 (--un-conc)
│ CGF catalogue + HumanMycobiomeScan index
│
├── Step 5: Probabilistic viral quantification
│ Salmon EM (--validateMappings,
│ --rangeFactorizationBins 4)
│ Output: salmon_quant/quant.sf
│
└── Step 6: Detection filtering and report generation
NumReads >= minreads AND TPM >= mintpm
Human blacklist applied (EVEs removed)
Output: detected_viruses.txt, REPORT.txt
cat << 'README' > ~/ViromeScan_v2/README.md
# ViromeScan v2.2

**A multi-host, probabilistic pipeline for virome detection and profiling from shotgun metagenomic data.**

Developed by [Hamidreza Chehrehnejad](https://github.com/HamidrezaChehrehnejad) — MSc Molecular Biology and Genetics, University of Pavia (2025/2026).
Supervised by Prof. Davide Sassera (University of Pavia) and Prof. Simone Rampelli (University of Bologna).

## Overview

ViromeScan v2.2 addresses the two core limitations of conventional virome profiling tools:

1. **Multi-mapping false positives:** resolved by replacing Bowtie2 best-hit assignment with **Salmon Expectation-Maximization (EM)** probabilistic quantification.
2. **Human-centric databases:** resolved by integrating host-specific depletion and viral indices for **human, cattle, pig, and grapevine** under a One Health framework.

## Pipeline

Raw reads (FASTQ)
│
├── Step 1: Broad viral pre-screen
│ Bowtie2 (--sensitive-local --no-unal)
│ Retains only viral candidate reads
│
├── Step 2: Host genome depletion
│ Bowtie2 (--un-conc)
│ Human: GRCh38
│ Cattle: Bos taurus
│ Pig: Sus scrofa
│ Grapevine: Vitis vinifera
│
├── Step 3: Microbiome depletion
│ Bowtie2 (--un-conc)
│ Human: HRGMv2 hybrid index
│ Other hosts: host-specific MAG databases
│
├── Step 4: Fungal depletion (human only)
│ Bowtie2 (--un-conc)
│ CGF catalogue + HumanMycobiomeScan index
│
├── Step 5: Probabilistic viral quantification
│ Salmon EM (--validateMappings,
│ --rangeFactorizationBins 4)
│ Output: salmon_quant/quant.sf
│
└── Step 6: Detection filtering and report generation
NumReads >= minreads AND TPM >= mintpm
Human blacklist applied (EVEs removed)
Output: detected_viruses.txt, REPORT.txt


## Dependencies

| Tool      | Version tested |
|-----------|----------------|
| Bowtie2   | >= 2.5         |
| samtools  | >= 1.17        |
| Salmon    | >= 1.10        |
| R         | >= 4.3         |
| tximport  | >= 1.28        |
| vegan     | >= 2.6         |
| tidyverse | >= 2.0         |

## Installation

```bash
git clone https://github.com/HamidrezaChehrehnejad/ViromeScan_v2.git
cd ViromeScan_v2
chmod +x viromescan_v2.2.sh
```

> **Note:** Databases (Bowtie2 indices and Salmon index) must be built separately and placed under `database/bowtie2/` and `database/salmon_index/`. See the thesis for database construction details.

## Usage

```bash
Paired-end, human host

./viromescan_v2.2.sh
-1 sample_R1.fastq.gz
-2 sample_R2.fastq.gz
--host human
-o results/sample_name
-p 16
Single-end, cattle host

./viromescan_v2.2.sh
-1 sample.fastq.gz
--host cattle
-o results/cattle_sample
-p 8



### Parameters

| Flag              | Description                                 | Default  |
|-------------------|---------------------------------------------|----------|
| `-1, --input1`    | R1 (or single-end) FASTQ (.gz accepted)     | required |
| `-2, --input2`    | R2 FASTQ (paired-end only)                  | optional |
| `--host`          | `human`, `cattle`, `pig`, or `grapevine`    | required |
| `-o, --output`    | Output directory                            | required |
| `-p, --threads`   | Number of threads                           | 8        |
| `--minreads`      | Min Salmon NumReads for detection           | 50       |
| `--mintpm`        | Min TPM for detection                       | 1.0      |
| `--salmon-index`  | Override default Salmon index path          | auto     |

## Output Files

| File                          | Description                                        |
|-------------------------------|----------------------------------------------------|
| `detected_viruses.txt`        | Filtered detections (NumReads >= minreads, TPM >= mintpm) |
| `detected_viruses_all.txt`    | All taxa with > 0 assigned reads                   |
| `salmon_quant/quant.sf`       | Full Salmon quantification output                  |
| `REPORT.txt`                  | Summary report with stats and top detections       |

## Downstream R Analysis

The `analysis/hadza_italian_virome/` folder contains the complete R pipeline used for the Hadza vs. Italian virome comparison (PRJNA278393, n = 20).

analysis/
└── hadza_italian_virome/
│
├── R/
│ ├── 01_tximport_abundance_matrix.R # Build TPM matrix from quant.sf files
│ ├── 02_differential_abundance.R # Wilcoxon rank-sum test + BH correction
│ ├── 03_virome_plots_species_names.R # Composition barplot + significant boxplots
│ ├── 04_beta_diversity_braycurtis.R # Bray-Curtis dissimilarity + PCoA
│ └── 05_run_all.R # Driver: sources scripts 01 through 04
│
└── figures/
├── Fig_Barplot_ViromeComposition_SpeciesNames.pdf
├── Fig_Top3_Boxplots_SpeciesNames.pdf
└── Fig_BetaDiversity_BrayCurtis_PCoA.pdf

To reproduce the full analysis, set the `PROJ` path in each script and run:

```r
source("analysis/hadza_italian_virome/R/05_run_all.R")
```

## Citation

If you use ViromeScan v2.2, please cite:

> Chehrehnejad H. (2026). *Development and Validation of a New Bioinformatics Pipeline for Detection and Profiling of the Virome from Shotgun Metagenomic Data.* MSc Thesis, University of Pavia.

## License

MIT License. See `LICENSE` for details.
