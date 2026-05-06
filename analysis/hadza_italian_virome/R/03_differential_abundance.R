source("R/00_config.R"); source("R/01_build_matrix.R")
kw <- vector("numeric", nrow(rel_matrix))
for (i in 1:nrow(rel_matrix))
    kw[i] <- kruskal.test(as.numeric(rel_matrix[i,]) ~ group)$p.value
names(kw) <- rownames(rel_matrix)
cat("Significant viruses (p<0.05):", sum(kw < 0.05), "
")
sig_table <- data.frame(Accession=names(kw[kw<0.05]), p.value=round(kw[kw<0.05],5))
sig_table <- sig_table[order(sig_table$p.value), ]
write.csv(sig_table, file.path(OUTPUT_DIR, "significant_viruses.csv"), row.names=FALSE)
print(sig_table)
pdf(file.path(OUTPUT_DIR, "boxplots_all_viruses.pdf"), width=10, height=12)
par(mfrow=c(4,2), mar=c(4,12,2,2))
for (i in 1:nrow(rel_matrix)) {
    boxplot(as.numeric(rel_matrix[i, group=="Hadza"]),
            as.numeric(rel_matrix[i, group=="Italian"]),
            las=1, col=c("forestgreen","steelblue"),
            names=c("Hadza","Italian"), horizontal=TRUE, lwd=0.8,
            xlab=if (kw[i]<0.05) paste("p =", signif(kw[i],2)) else "",
            main=list(rownames(rel_matrix)[i], cex=0.9, font=4))
}
dev.off()
cat("Saved: boxplots_all_viruses.pdf
")
