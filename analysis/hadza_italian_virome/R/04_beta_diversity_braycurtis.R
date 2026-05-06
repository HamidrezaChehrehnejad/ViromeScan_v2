library(tidyverse); library(vegan)
PROJ <- "/home/hamidreza/Virome_Project"
tpm  <- read.csv(file.path(PROJ, "results/abundance_matrix_raw.csv"), row.names = 1)
meta <- read.csv(file.path(PROJ, "data/metadata/metadata.csv"))
rownames(meta) <- meta$sample_id
tpm <- tpm[, meta$sample_id]
bray_dist <- vegdist(t(tpm), method = "bray")
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(Sample = rownames(pcoa$points), PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2], group = meta[rownames(pcoa$points), "group"])
p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, colour = group)) + geom_point(size = 3, alpha = 0.9) + stat_ellipse(aes(fill = group), geom = "polygon", alpha = 0.15, colour = NA, show.legend = FALSE) +
  scale_colour_manual(values = c("Hadza" = "#2196F3", "Italian" = "#FF9800")) + theme_bw() + theme(legend.position = "bottom")
ggsave(file.path(PROJ, "results/figures/Fig_BetaDiversity_BrayCurtis_PCoA.pdf"), p, width = 6, height = 5)
