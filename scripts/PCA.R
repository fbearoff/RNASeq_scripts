# Inputs ------------------------------------------------------------------
WD <- ""
project <- "" # subdir of WD
data_dir <- "" #Salmon output dir

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
samples <- read.table(file.path(WD, project, "samples.txt"), colClasses = 'character', header = TRUE)
files <- file.path(WD, data_dir, "quants", samples$sample.id, "quant.sf")
names(files) <- samples$name
all(file.exists(files))
tx2gene <- read.csv(file.path(WD, "tx2gene.csv"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# DESeq2 ------------------------------------------------------------------
sampleTable <- data.frame(condition = samples$condition)
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
dds <- DESeq(dds, parallel = TRUE)
gene_synonym <- unique(tx2gene[, -1])

# Clustered Sample Distance Plot ------------------------------------------
rld <- rlog(dds, blind=FALSE)
# PCA Plot ------------------------------------------------------------
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
pcaData$run <- as.factor(samples$run) #other criteria important for sample separation
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file.path(WD, project, "PCA.pdf"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=run)) + #shape is other criteria
  geom_point(size=2) + geom_text_repel(aes(label=name)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()
dev.off()
