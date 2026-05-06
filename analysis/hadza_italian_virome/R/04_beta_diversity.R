source("R/00_config.R"); source("R/01_build_matrix.R")
x.beta <- vegdist(t(rel_matrix), "bray")
vare.x <- capscale(as.dist(x.beta) ~ 1)
coords <- as.data.frame(scores(vare.x, display="sites", choices=c(1,2)))
colore <- ifelse(group == "Hadza", "forestgreen", "steelblue")
eig     <- vare.x$CA$eig
var_exp <- round(eig / sum(eig) * 100, 1)
pdf(file.path(OUTPUT_DIR, "PCoA_BrayCurtis.pdf"), width=8, height=7)
plot(vare.x, type="n", main="Beta Diversity - Bray-Curtis PCoA", choices=c(1,2),
     xlab=paste0("MDS1 (", var_exp[1], "%)"),
     ylab=paste0("MDS2 (", var_exp[2], "%)"))
ordiellipse(vare.x, choices=c(1,2), group, kind="se", conf=0.95,
            show.groups="Hadza",   col="forestgreen", lwd=2)
ordiellipse(vare.x, choices=c(1,2), group, kind="se", conf=0.95,
            show.groups="Italian", col="steelblue",   lwd=2)
ordispider(vare.x, choices=c(1,2), group, col="forestgreen",
           show.groups="Hadza",   label=FALSE, lwd=1.5)
ordispider(vare.x, choices=c(1,2), group, col="steelblue",
           show.groups="Italian", label=FALSE, lwd=1.5)
with(coords, points(MDS1, MDS2, pch=21, bg=colore, cex=1.5))
legend("topright", legend=c("Hadza","Italian"), fill=c("forestgreen","steelblue"))
dev.off()
set.seed(42)
perm_result <- adonis2(x.beta ~ group, permutations=9999)
print(perm_result)
capture.output(print(perm_result), file=file.path(OUTPUT_DIR, "PERMANOVA_result.txt"))
cat("Saved: PCoA_BrayCurtis.pdf and PERMANOVA_result.txt
")
