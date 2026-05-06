library(RColorBrewer)
library(vegan)
BASE_PATH  <- "/home/hamidreza/Virome_Project/results/viromescan"
OUTPUT_DIR <- "output"
dir.create(OUTPUT_DIR, showWarnings = FALSE)
hadza_samples   <- paste0("H",  c(1,2,3,4,5,6,7,8,9,10))
italian_samples <- paste0("IT", c(1,2,3,4,5,6,7,8,11,13))
all_samples     <- c(hadza_samples, italian_samples)
group <- as.factor(c(rep("Hadza", length(hadza_samples)),
                     rep("Italian", length(italian_samples))))
names(group) <- all_samples
GROUP_COLORS <- c("Hadza"="forestgreen","Italian"="steelblue")
