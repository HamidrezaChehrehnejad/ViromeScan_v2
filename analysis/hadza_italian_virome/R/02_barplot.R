source("R/00_config.R"); source("R/01_build_matrix.R")
top_n <- 30
top_species <- names(sort(rowSums(rel_matrix), decreasing=TRUE))[1:top_n]
plot_data   <- rel_matrix[top_species, ]
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector    <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(42); colori <- sample(col_vector, top_n)
pdf(file.path(OUTPUT_DIR, "barplot_virome.pdf"), width=12, height=7)
par(mar=c(9, 5, 4, 8))
barplot(as.matrix(plot_data), col=colori, lwd=0.6, border=FALSE, space=0,
        ylab="Relative Abundance", las=2,
        main=list("Gut Virome - Hadza vs Italian", cex=1.3, font=4))
abline(v=10, lwd=2, lty=2, col="black")
mtext("Hadza",   side=3, at=5,  font=2, cex=1.1)
mtext("Italian", side=3, at=15, font=2, cex=1.1)
plot(as.matrix(plot_data), type="n", axes=FALSE, xlab="", ylab="")
legend("center", legend=rownames(plot_data), fill=colori,
       horiz=FALSE, text.font=2, cex=0.75, ncol=2)
dev.off()
cat("Saved: barplot_virome.pdf
")
