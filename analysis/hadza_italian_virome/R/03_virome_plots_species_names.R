# 03_virome_plots_species_names.R — Composition barplot + significant virus boxplots
library(tidyverse)

PROJ <- "/home/hamidreza/Virome_Project"

tpm  <- read.csv(file.path(PROJ, "results/abundance_matrix_raw.csv"), row.names = 1)
meta <- read.csv(file.path(PROJ, "data/metadata/metadata.csv"))
rownames(meta) <- meta$sample_id
sig  <- read.csv(file.path(PROJ, "results/significant_viruses.csv"))

# Load accession -> species name map (TSV: accession \t full title)
# File: data/metadata/accession_names.tsv  (two columns, no header)
# If missing, accession IDs are used as fallback labels.
name_map_path <- file.path(PROJ, "data/metadata/accession_names.tsv")
if (file.exists(name_map_path)) {
  name_map <- read.delim(name_map_path, header = FALSE,
                         col.names = c("Accession", "FullTitle")) %>%
    mutate(Species = str_remove(FullTitle, ",.*") %>%
                     str_remove("\\s*(complete|partial).*") %>%
                     str_trim())
  acc2species <- setNames(name_map$Species, name_map$Accession)
} else {
  warning("accession_names.tsv not found — using accession IDs as labels.")
  acc2species <- setNames(rownames(tpm), rownames(tpm))
}

dir.create(file.path(PROJ, "results/figures"), showWarnings = FALSE, recursive = TRUE)

# --- Boxplots of significant viruses ---
if (nrow(sig) > 0) {
  long_box <- tpm[sig$Accession, , drop = FALSE] %>%
    rownames_to_column("Accession") %>%
    pivot_longer(-Accession, names_to = "Sample", values_to = "TPM") %>%
    mutate(
      Group   = meta[Sample, "group"],
      logTPM  = log2(TPM + 1),
      Species = ifelse(is.na(acc2species[Accession]), Accession, acc2species[Accession])
    )

  p1 <- ggplot(long_box, aes(x = Group, y = logTPM, fill = Group)) +
    geom_boxplot(outlier.shape = 21, width = 0.6) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
    facet_wrap(~ Species, scales = "free_y") +
    scale_fill_manual(values = c("Hadza" = "#2196F3", "Italian" = "#FF9800")) +
    labs(y = "log2(TPM + 1)", x = NULL, fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(PROJ, "results/figures/Fig_Top3_Boxplots_SpeciesNames.pdf"),
         p1, width = 8, height = 4)
  cat("[03] Boxplots saved.\n")
} else {
  cat("[03] No significant viruses — boxplot skipped.\n")
}

# --- Stacked barplot of relative abundance (top 15 + Other) ---
top15    <- names(sort(rowMeans(tpm), decreasing = TRUE)[1:15])
tpm_plot <- rbind(tpm[top15, ], Other = colSums(tpm) - colSums(tpm[top15, ]))
tpm_rel  <- sweep(tpm_plot, 2, colSums(tpm_plot), "/") * 100

rownames(tpm_rel) <- ifelse(
  rownames(tpm_rel) == "Other", "Other",
  ifelse(is.na(acc2species[rownames(tpm_rel)]),
         rownames(tpm_rel),
         acc2species[rownames(tpm_rel)])
)

long_bar <- tpm_rel %>%
  rownames_to_column("Virus") %>%
  pivot_longer(-Virus, names_to = "Sample", values_to = "RelAb") %>%
  mutate(
    Group = meta[Sample, "group"],
    Virus = factor(Virus, levels = rownames(tpm_rel))
  )

p2 <- ggplot(long_bar, aes(x = Sample, y = RelAb, fill = Virus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  labs(y = "Relative Abundance (%)", x = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text     = element_text(size = 7)
  )

ggsave(file.path(PROJ, "results/figures/Fig_Barplot_ViromeComposition_SpeciesNames.pdf"),
       p2, width = 14, height = 7)

cat("[03] Barplot saved.\n")
