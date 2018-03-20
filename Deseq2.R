# Inputs ------------------------------------------------------------------
WD <- ""
condition1 <- ""
condition2 <- ""
CHR <- "all"   #specify individual chromosome, deafult is "all"
filter_data <- FALSE
abund_filter <- 1
cv_filter <- 100

# Load Libraries ----------------------------------------------------------
library(BiocParallel)
library(tximport)
library(readr)
library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
register(MulticoreParam(10))

# TxImport ----------------------------------------------------------------
setwd(WD)
comparison <- paste(condition1, "_vs_", condition2, sep="")
samples <- read.table(file.path(WD, "samples.txt"), header = TRUE)
samples <- samples[ which(samples$condition == condition1 | samples$condition == condition2), ]
files <- file.path(WD, "quants", samples$location, "quant.sf")
names(files) <- samples$sample
all(file.exists(files))
tx2gene <- read.csv(file.path(WD, "tx2gene.csv"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
txi.tx <- tximport(files, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)
write.csv(txi.tx, paste(comparison, ".txi.tx.csv", sep=""))

# DESeq2 ------------------------------------------------------------------
sampleTable <- data.frame(condition = samples$condition)
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds, contrast=c("condition", condition1, condition2), parallel = TRUE)
gene_synonym <- unique(tx2gene[, -1])
z <- data.frame(res)
z$gene_symbol <- gene_synonym$gene_symbol[match(rownames(res), gene_synonym$gene_id)]
z$gene_biotype <- gene_synonym$gene_biotype[match(rownames(res), gene_synonym$gene_id)]
z$mgi_id <- gene_synonym$mgi_id[match(rownames(res), gene_synonym$gene_id)]
z$chr <- gene_synonym$chr[match(rownames(res), gene_synonym$gene_id)]
z$start <- gene_synonym$start[match(rownames(res), gene_synonym$gene_id)]
z$end <- gene_synonym$end[match(rownames(res), gene_synonym$gene_id)]
z$description <- gene_synonym$description[match(rownames(res), gene_synonym$gene_id)]
if(CHR!="all") {
  z<- z [ which(z$chr == CHR), ]
  }
write.csv(merge(z,txi, by=0), paste(comparison, ".res.csv", sep=""))

# Expression Profile Plot -------------------------------------------------
p <- ggplot(z, aes(x = baseMean + 0.01, y = log2FoldChange,
                   col = paste(padj < 0.05, is.na(padj)),
                   shape = paste(padj < 0.05, is.na(padj))))
p + geom_point() +
  scale_x_log10() +
  scale_y_continuous(limits = c(-2, 2)) +
  scale_color_manual(labels=c("p>0.05", "Below threshold", "p<0.05"),
                     values = c("lightblue", "salmon", "green"),
                     guide_legend(title="")) +
  scale_shape_manual(guide="none", values = c(0, 1, 2)) +
  labs(title= paste(paste(condition1, "vs.", condition2, "\n"),
                    "Expression Profile", sep=""),
                    x="Expression") +
  guides(colour = guide_legend(override.aes = list(shape = c(0,1,2))))

# Volcano Plot ------------------------------------------------------------
p <- ggplot(z, aes(x = log2FoldChange, y = 1 - padj, col = abs(log2FoldChange)))
p + geom_point(size = 1) + scale_y_log10()

# Clustered Sample Distance Plot ------------------------------------------
rld <- rlog(dds, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste0(samples$sample, rld$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# PCA Plot ------------------------------------------------------------
pca_plot_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pca_plot_data, "percentVar"))
ggplot(pca_plot_data, aes(PC1, PC2, color=condition)) +
  geom_point(size=2) + geom_text_repel(aes(label=name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_color_manual(values=c("green", "red", "purple"))

# Data Filtering for Heatmap ----------------------------------------------
if(filter_data==TRUE) {
  filtered_samples <- txi
  filtered_samples$control_abund_avg <-apply(filtered_samples$abundance[, c(1:3)], 1, FUN = mean)
  filtered_samples$control_abund_sd <-apply(filtered_samples$abundance[, c(1:3)], 1, FUN = sd)
  filtered_samples$control_abund_cv <- (filtered_samples$control_abund_sd / filtered_samples$control_abund_avg) * 100
  filtered_samples2 <- as.data.frame(which(filtered_samples$control_abund_avg >=abund_filter & filtered_samples$control_abund_cv <=cv_filter))
  z <- z[which(rownames(z) %in% rownames(filtered_samples2)), ]
}

# Clustered Top 50 Heatmap ------------------------------------------------
top.count <- 100
top.genes <- rownames(z)[order(z$padj)][1:top.count]
abund <- txi$abundance
abund <- abund[which(rownames(abund) %in% top.genes), ]
abund.scale <- t(scale(t(abund), center = TRUE, scale = TRUE))
range(abund.scale)
heatmap.labels <- z$gene_symbol[which(rownames(z) %in% top.genes)]
no.symbol.ind <- which(is.na(heatmap.labels))
heatmap.labels[no.symbol.ind] <- z$gene_biotype[which(rownames(z) %in% top.genes)][no.symbol.ind]
samplecol = character(length = length(samples$condition))
for (i in which(samples$condition==condition1)) {
  samplecol[i] = "forestgreen"
}
for (i in which(samples$condition==condition2)) {
  samplecol[i] = "red"
}
pdf(file = paste(comparison,".heatmap.pdf", sep=""))
heatmap.2(abund.scale, trace = "none", Colv = TRUE, margins = c(7, 10),
          breaks = seq(-2, 2, by = 4/11), col = brewer.pal(11, "BrBG"),
          labRow = heatmap.labels, cexRow = 0.8, cexCol = 1, ColSideColors = samplecol
)
dev.off()
