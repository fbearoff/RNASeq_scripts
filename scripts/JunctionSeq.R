# Inputs ------------------------------------------------------------------
QC_dir <- "/home/frank/analyzed/RNASeq1_call_variants/QC"
WD <- "/home/frank/R_projects/RNASeq1"
gene_names <- "/home/frank/GRCm38_construct/gene_names.csv"
condition1 <- "B6"
condition2 <- "103"
cores = 24
# Load Libraries ----------------------------------------------------------
library(JunctionSeq)
library(BiocParallel)
# Load samples ------------------------------------------------------------
setwd(WD)
gene_names <- read.csv(gene_names, header = TRUE)
comparison <- paste(condition1, "_vs_", condition2, sep = "")
samples <- read.table(file.path(WD, "samples.txt"),
                      header = TRUE)
samples <-
  samples[which(samples$condition == condition1 |
                  samples$condition == condition2),]
GFF <- file.path(QC_dir, "withNovel.forJunctionSeq.gff.gz")
countFiles <-
  file.path(
    QC_dir,
    samples$sample.ID,
    "QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz"
  )
# JunctionSeq -------------------------------------------------------------
jscs <- runJunctionSeqAnalyses(
  sample.files = countFiles,
  sample.names = samples$sample.ID,
  condition = factor(samples$condition),
  flat.gff.file = GFF,
  nCores = cores,
  analysis.type = "junctionsAndExons",
  gene.names = gene_names
)

dir.create(path =  paste0(WD, "/JunctionSeq/", comparison, "/results"),
           recursive = TRUE)

writeCompleteResults(jscs,
                     outfile.prefix = paste0(WD,"/JunctionSeq/", comparison, "/results/"),
                     save.jscs = TRUE)

buildAllPlots(
  jscs = jscs,
  outfile.prefix = paste0(WD,"/JunctionSeq/", comparison, "/"),
  use.plotting.device = "png",
  FDR.threshold = 0.01,
  colorRed.FDR.threshold = 0.05
)