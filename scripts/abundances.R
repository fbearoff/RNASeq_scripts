# Inputs ------------------------------------------------------------------
WD <- ""
project <- "" #becomes subdir of WD
dataDir <- "" #Salmon output dir

# Load Libraries ----------------------------------------------------------
library(BiocParallel)
library(tximport)
library(readr)
register(MulticoreParam(10))

# TxImport ----------------------------------------------------------------
setwd(WD)
dir.create(project)
samples <- read.table(file.path(WD,  "samples.txt"), colClasses = 'character', header = TRUE)
files <- file.path(WD, dataDir, "quants", samples$sample.id, "quant.sf")
names(files) <- samples$name
all(file.exists(files))
tx2gene <- read.csv(file.path(WD, "tx2gene.csv"))
gene_synonym <- unique(tx2gene[, -1])
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

txi_df <- data.frame(txi)

txi_df$gene_symbol <- gene_synonym$gene_symbol[match(rownames(txi_df), gene_synonym$gene_id)]
txi_df$gene_biotype <- gene_synonym$gene_biotype[match(rownames(txi_df), gene_synonym$gene_id)]
txi_df$mgi_id <- gene_synonym$mgi_id[match(rownames(txi_df), gene_synonym$gene_id)]
txi_df$chr <- gene_synonym$chr[match(rownames(txi_df), gene_synonym$gene_id)]
txi_df$start <- gene_synonym$start[match(rownames(txi_df), gene_synonym$gene_id)]
txi_df$end <- gene_synonym$end[match(rownames(txi_df), gene_synonym$gene_id)]
txi_df$description <- gene_synonym$description[match(rownames(txi_df), gene_synonym$gene_id)]

txi.tx <- tximport(files, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)
write.csv(txi.tx, file.path(WD, project, paste("abundances", ".txi.tx.csv", sep="")))
write.csv(txi_df, file.path(WD, project, paste("abundances", ".txi.csv", sep="")))



