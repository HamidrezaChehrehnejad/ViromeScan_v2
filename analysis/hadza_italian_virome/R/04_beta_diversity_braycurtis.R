# 04_beta_diversity_braycurtis.R — Bray-Curtis dissimilarity + PCoA
library(tidyverse)
library(vegan)

PROJ <- "/home/hamidreza/Virome_Project"

tpm  <- read.csv(file.path(PROJ, "results/abundance_matrix_raw.csv"), row.names = 1)
meta <- read.csv(file.path(PROJ, "data/metadata/metadata.csv"))
rownames(meta) <- meta$sample_id

# Align column order to metadata
tpm <- tpm[, meta$sample_id]

# Bray-Curtis distance (samples in rows)
bray_dist <- vegdist(t(tpm), method = "bray")

# PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
var_expl <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100

pcoa_df <- data.frame(
  Sample = rownames(pcoa$points),
  PC1    = pcoa$points[, 1],
  PC2    = pcoa$points[, 2],
  group  = meta[rownames(pcoa$points), "group"]
)

p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, colour = group)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(fill = group), geom = "polygon",
               alpha = 0.15, colour = NA, show.legend = FALSE) +
  scale_colour_manual(values = c("Hadza" = "#2196F3", "Italian" = "#FF9800")) +
  scale_fill_manual(values   = c("Hadza" = "#2196F3", "Italian" = "#FF9800")) +
  labs(
    x      = paste0("PC1 (", round(var_expl[1], 1), "%)"),
    y      = paste0("PC2 (", round(var_expl[2], 1), "%)"),
    colour = NULL,
    title  = "Bray\u2013Curtis PCoA of Virome Composition"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid      = element_blank()
  )

dir.create(file.path(PROJ, "results/figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(PROJ, "results/figures/Fig_BetaDiversity_BrayCurtis_PCoA.pdf"),
       p, width = 6, height = 5)

cat("[04] Bray-Curtis PCoA saved.\n")
