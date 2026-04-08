# ViromeScan v2.2
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-2.2-blue.svg)]()
[![Conda](https://img.shields.io/badge/conda-compatible-green.svg)](https://docs.conda.io/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20Ubuntu-lightgrey.svg)]()
A multi-host bioinformatics pipeline for metagenomic viral detection and
quantification using EM-based read assignment. Developed as part of a
Master's Thesis in Molecular Biology and Genetics at Università degli studi di Pavia (Dipartamento di Biologia e Biotechnologia)/Alma mater studiorum - Università di Bologna (Dipartamento di Farmacia e Biotechnologia).
## How It Works
The pipeline runs six sequential steps on paired-end or single-end Illumina data:

Note: for human samples, a curated blacklist of endogenous viral elements (EVEs)
and integrated sequences (e.g. NC_022518.1, AC000005.1) is automatically applied
at Step 6 to suppress false positives.
## Prerequisites
### System Requirements
| Resource   | Minimum      | Recommended  |
|------------|--------------|--------------|
| RAM        | 16 GB        | 32 GB+       |
| CPU Cores  | 4            | 8+           |
| Disk Space | 60 GB       | 200 GB+      |
| OS         | Ubuntu 20.04+| Ubuntu 22.04 |
### Dependencies (managed via Conda)
| Tool     | Version | Role |
|----------|---------|------|
| Bowtie2  | 2.4.5+  | All filtering steps and PhiX removal |
| Salmon   | 1.9.0+  | EM-based viral quantification |
| Samtools | 1.15+   | SAM/BAM to FASTQ conversion |
| Seqkit   | 2.3.0+  | FASTQ manipulation |
| Python   | 3.9+    | Output parsing |
## Installation
bash
git clone https://github.com/YourUsername/ViromeScan2.git
cd ViromeScan2
conda env create -f envs/viromescan_v2.yaml
conda activate viromescan_v2

## Usage
bash
bash bin/viromescan_v2.2.sh \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    --host human \
    -o results/sample_out/ \
    -p 8

Single-end mode is also supported by omitting -2.
### Parameters
| Flag              | Description                                         | Default |
|-------------------|-----------------------------------------------------|---------|
| -1, --input1    | R1 FASTQ file (.gz accepted)                        | required |
| -2, --input2    | R2 FASTQ file, paired-end only                      | optional |
| --host          | Host: human, cattle, pig, grapevine         | required |
| -o, --output    | Output directory                                    | required |
| -p, --threads   | CPU threads                                         | 8 |
| --minreads      | Minimum Salmon NumReads threshold for detection     | 50 |
| --mintpm        | Minimum TPM threshold for detection                 | 1.0 |
| --salmon-index  | Override default Salmon index path                  | auto |
Note: human uses a curated Salmon index; all other hosts use the full viral database.
## Output

Intermediate step files (step1–step4 FASTQs) are automatically removed
after the run to save disk space.
## Validation
Validated against 5 synthetic mock communities across human, cattle and pig
hosts spanning 3 taxonomic levels. Key results:
- Precision: >92% across all replicates
- Recall: >88% at species level (50-read threshold)
- False Positive Rate: <5% with combined NumReads + TPM filtering
See /validation for mock community configurations and full benchmark scripts.
## Case Study
Applied to 10 Hadza hunter-gatherer and 10 urban Italian gut metagenomic
samples (NCBI BioProject: PRJNA278393). Key findings:
- Illumina PhiX174 spike-in artifact (up to 96% of viral reads in some
  samples) identified and removed via upstream Bowtie2 filtering.
- Only 4 of 105 detected viral species shared across both cohorts,
  reflecting extreme population-level virome specificity.
- Reference database bias caused ~13% viral read loss in Italian samples
  vs ~7.5% in Hadza due to prophage cross-mapping.
See /case_study for R analysis scripts and figures.
## How to Cite
> Chehrehnejad, H. (2026). ViromeScan v2.2: An optimized multi-host pipeline
> for metagenomic viral detection and quantification. Master's Thesis,
> Università degli studi di Pavia/ Alma mater studiorum - Università di Bologna.
## License
MIT License. See LICENSE file for details.
## Author
Hamidreza Chehrehnejad
Master's Thesis — Università degli studi di Pavia/ Alma mater studiorum - Università di Bologna
