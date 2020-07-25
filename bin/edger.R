#!/usr/bin/env Rscript
library(pheatmap)
library(edgeR)
library(ggplot2)
library(matrixStats)
library("org.Mm.eg.db")
library(corrplot)
library(RColorBrewer)
args <- commandArgs(trailingOnly = TRUE)
x <- read.delim(args[1], row.names = "Geneid")
# x <- read.delim("nf_output/featureCounts/merged_gene_counts.txt", row.names="Geneid")
x[, 1] <- NULL
a <- x

# x <- x[, c(1, 4, 5, 3)]

gene_length <- read.table("/mnt/d/ucsc/mm10/gene_length.txt", header = T)

# y <- DGEList(counts = x, group = c(1, 1, 2, 2), genes = gene_length)
# keep <- filterByExpr(y)
# y <- y[keep, keep.lib.sizes = FALSE]
# y <- calcNormFactors(y)

expgroup <- as.factor(c(
    "WT_POS",
    "WT_NEG",
    "CD_POS",
    "WT_POS",
    "CD_POS",
    "CD_NEG",
    "WT_NEG",
    "CD_NEG"
))
design <- model.matrix(~ 0 + expgroup)
ya <- DGEList(counts = a, group = expgroup, genes = gene_length)
keep <- filterByExpr(ya)
ya <- ya[keep, keep.lib.sizes = FALSE]
ya <- calcNormFactors(ya)

# mds
p <- plotMDS(ya)
dev.off()
mds <- as.data.frame(cbind(p$x, p$y))
colnames(mds) <- c("dim1", "dim2")

mds$Group <- expgroup

pdf("mds.pdf", width = 4, height = 3)
ggplot(mds, aes(x = dim1, y = dim2, color = Group)) +
    geom_point() +
    xlab("Leading logFC dim 1") +
    ylab("Leading logFC dim 2") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()


# sample correlation
cormat <- cor(x)
pdf("sample_correlation.pdf")
corrplot(cormat, method = "square", order = "hclust", hclust.method = "ward.D2")
dev.off()

# deg
ya <- estimateDisp(ya, design)
fit <- glmQLFit(ya, design, robust = TRUE)


con_pos <- makeContrasts(expgroupWT_POS - expgroupCD_POS, levels = design)
con_neg <- makeContrasts(expgroupWT_NEG - expgroupCD_NEG, levels = design)
qlf_pos <- glmQLFTest(fit, contrast = con_pos)
qlf_neg <- glmQLFTest(fit, contrast = con_neg)

qlf_pos_sig <- topTags(qlf_pos, n = 20000, p.value = 0.05)
qlf_neg_sig <- topTags(qlf_neg, n = 20000, p.value = 0.05)



# heatmap
logrpkm <- rpkm(ya, log = TRUE)

plot_and_go <- function(logrpkm, filename, cutree) {
    pdf(filename, width = 10, height = 2.5)
    ph = pheatmap(t(logrpkm), colorRampPalette(rev(brewer.pal(
        n = 7, name = "RdBu")))(100), border_color = NA, breaks = seq(-1, 1, 0.02), scale = "column", show_colnames = F, treeheight_col = 0, clustering_method = "ward.D2", cutree_cols = cutree)
    dev.off()

    cutt = cutree(ph$tree_col, k = cutree)

    return(lapply(seq(cutree), function(x) names(cutt[cutt==x])))
}





pdf("heatmap_pos_comparison.pdf", width = 10, height = 2.5)
ph_pos = pheatmap(t(logrpkm[rownames(qlf_pos_sig),]), colorRampPalette(rev(brewer.pal(
    n = 7, name =
        "RdBu"
)))(100), border_color = NA, breaks = seq(-1, 1, 0.02), scale = "column", show_colnames = F, treeheight_col = 0, clustering_method = "ward.D2", cutree_cols = 5)
dev.off()

pdf("heatmap_neg_comparison.pdf", width = 10, height = 2.5)
ph_neg = pheatmap(t(logrpkm[rownames(qlf_neg_sig),]), colorRampPalette(rev(brewer.pal(
    n = 7, name =
        "RdBu"
)))(100), border_color = NA, breaks = seq(-1, 1, 0.02), scale = "column", show_colnames = F, treeheight_col = 0, clustering_method = "ward.D2", cutree_cols = 7)
dev.off()




pdf("logfc_logcpm.pdf")
plotMD(lrt)
abline(h = c(-1, 1), col = "blue")
dev.off()

eid <- mapIds(org.Mm.eg.db, rownames(topTags(lrt, n = 100)), "ENTREZID", "SYMBOL")

go <- goana(eid, universe = universe, species = "Mm")

all_gene <- merge(lrt$table, a, by = "row.names", all.x = TRUE)
colnames(all_gene)[1] <- "gene"

write.csv(all_gene, "all_genes.csv", quote = FALSE, row.names = FALSE)

deg <- merge(topTags(lrt, n = 100), a, by = "row.names", all.x = TRUE)
colnames(deg)[1] <- "gene"
write.csv(deg, "top100_deg.csv", quote = FALSE, row.names = FALSE)

write.csv(go, "go.csv", quote = FALSE)
write.csv(go[go$N < 200 & go$P.DE < (0.05 / nrow(go)), ], "go_N200_Bonferroni_sig.csv", quote = FALSE)