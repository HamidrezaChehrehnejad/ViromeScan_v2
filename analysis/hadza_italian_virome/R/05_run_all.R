# 05_run_all.R — Driver script: runs the full Hadza vs Italian virome analysis
PROJ <- "/home/hamidreza/Virome_Project"
setwd(file.path(PROJ, "analysis/hadza_italian_virome"))

source("R/01_tximport_abundance_matrix.R")
source("R/02_differential_abundance.R")
source("R/03_virome_plots_species_names.R")
source("R/04_beta_diversity_braycurtis.R")

cat("\n=== Analysis complete. All outputs saved in results/figures/ ===\n")
