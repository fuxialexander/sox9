#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
x <- read.delim(args[1], row.names = "Geneid")
x[, 1] <- NULL
x <- x[, c(8, 1, 3, 4)]
library(edgeR)
y <- DGEList(counts = x, group = c(1, 1, 2, 2))

keep <- filterByExpr(y)
y <- y[keep, keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
# pdf("mds.pdf")
# plotMDS(y, col = 1:4)
# dev.off()

cormat <- cor(x)
library(corrplot)
pdf("sample_correlation.pdf")
corrplot(cormat, method = "square", order = "hclust", hclust.method = "ward.D2")
dev.off()

library(ggplot2)
library(matrixStats)
rowv <- rowVars(as.matrix(log2(x)))
pdf("gene_variance_distribution.pdf")
qplot(rowv) + geom_density()
dev.off()

design <- model.matrix(~ factor(c(1, 1, 2, 2)))
y <- estimateDisp(y, design)


fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)

topTags(lrt, n = 100)

# pdf("bcv.pdf")
# plotBCV(y)
# dev.off()

pdf("logfc_logcpm.pdf")
plotMD(lrt)
abline(h = c(-1, 1), col = "blue")
dev.off()


library("org.Mm.eg.db")
eid <- mapIds(org.Mm.eg.db, rownames(topTags(lrt, n = 100)), "ENTREZID", "SYMBOL")
universe <- mapIds(org.Mm.eg.db, rownames(lrt$table), "ENTREZID", "SYMBOL")

go <- goana(eid, universe = universe, species = "Mm")

write.csv(x, "all_read_counts.csv", quote=FALSE)
write.csv(go, "go.csv", quote = FALSE)
write.csv(go[go$N < 200 & go$P.DE < (0.05 / nrow(go)), ], "go_N200_Bonferroni_sig.csv", quote = FALSE)
write.csv(topTags(lrt, n = 100), "top100_deg.csv", quote = FALSE)