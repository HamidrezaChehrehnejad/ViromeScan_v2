library(tidyverse)
PROJ <- "/home/hamidreza/Virome_Project"
tpm  <- read.csv(file.path(PROJ, "results/abundance_matrix_raw.csv"), row.names = 1)
meta <- read.csv(file.path(PROJ, "data/metadata/metadata.csv"))
rownames(meta) <- meta$sample_id
hadza   <- meta$sample_id[meta$group == "Hadza"]
italian <- meta$sample_id[meta$group == "Italian"]
results <- do.call(rbind, lapply(rownames(tpm), function(virus) {
  h <- as.numeric(tpm[virus, hadza]); i <- as.numeric(tpm[virus, italian])
  if (sum(c(h, i) > 0) < 4) return(NULL)
  test <- wilcox.test(h, i, exact = FALSE)
  data.frame(Accession = virus, p_value = test$p.value, median_Hadza = median(h), median_Italian = median(i))
}))
results$p_adj <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$p_value), ]
sig <- results[results$p_adj < 0.05, ]
write.csv(results, file.path(PROJ, "results/differential_abundance.csv"), row.names = FALSE)
write.csv(sig, file.path(PROJ, "results/significant_viruses.csv"), row.names = FALSE)
