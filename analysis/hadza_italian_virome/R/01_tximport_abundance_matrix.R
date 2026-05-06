# 01_tximport_abundance_matrix.R — Build TPM matrix from Salmon quant.sf files
library(tximport)

PROJ <- "/home/hamidreza/Virome_Project"

meta  <- read.csv(file.path(PROJ, "data/metadata/metadata.csv"))
rownames(meta) <- meta$sample_id

files <- file.path(PROJ, "results/viromescan", meta$sample_id, "salmon_quant/quant.sf")
names(files) <- meta$sample_id

txi <- tximport(files, type = "salmon", txOut = TRUE)
tpm <- txi$abundance

dir.create(file.path(PROJ, "results"), recursive = TRUE, showWarnings = FALSE)
write.csv(tpm, file.path(PROJ, "results/abundance_matrix_raw.csv"), quote = FALSE)

cat("[01] Abundance matrix saved.\n")
