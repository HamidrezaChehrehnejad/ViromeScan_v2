# 00_config.R — Shared configuration (sourced by individual scripts if needed)
# NOTE: Each script also sets PROJ independently for safe standalone execution.
library(RColorBrewer)
library(vegan)

PROJ <- "/home/hamidreza/Virome_Project"
OUTPUT_DIR <- file.path(PROJ, "results/figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

hadza_samples   <- paste0("H",  c(1,2,3,4,5,6,7,8,9,10))
italian_samples <- paste0("IT", c(1,2,3,4,5,6,7,8,11,13))
all_samples     <- c(hadza_samples, italian_samples)

group <- as.factor(c(rep("Hadza",   length(hadza_samples)),
                     rep("Italian", length(italian_samples))))
names(group) <- all_samples

GROUP_COLORS <- c("Hadza" = "#2196F3", "Italian" = "#FF9800")
