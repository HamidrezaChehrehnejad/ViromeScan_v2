source("R/00_config.R")
tpm_list <- lapply(all_samples, function(s) {
    f  <- file.path(BASE_PATH, s, "salmon_quant", "quant.sf")
    df <- read.delim(f, header=TRUE, stringsAsFactors=FALSE)
    tpm <- df$TPM; names(tpm) <- df$Name; tpm
})
names(tpm_list) <- all_samples
all_viruses <- unique(unlist(lapply(tpm_list, names)))
tpm_matrix  <- matrix(0, nrow=length(all_viruses), ncol=length(all_samples),
                      dimnames=list(all_viruses, all_samples))
for (s in all_samples) tpm_matrix[names(tpm_list[[s]]), s] <- tpm_list[[s]]
tpm_matrix  <- as.data.frame(tpm_matrix)
tpm_matrix  <- tpm_matrix[rowSums(tpm_matrix) > 0, ]
tpm_matrix  <- tpm_matrix[!rownames(tpm_matrix) %in% "NC_001422.1", ]
rel_matrix  <- sweep(tpm_matrix, 2, colSums(tpm_matrix), "/")
cat("Matrix ready:", nrow(rel_matrix), "viruses x", ncol(rel_matrix), "samples
")
