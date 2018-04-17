# Inputs ------------------------------------------------------------------
QC_dir <- "/home/frank/analyzed/ALS_RNASeq2_call_variants/QC"
WD <- "/scratch/frank/R_projects/ALS_RNASeq2"
condition1 <- "jax_BB60"
condition2 <- "SNO24"
cores = 32
# Load Libraries ----------------------------------------------------------
library(JunctionSeq)
library(BiocParallel)
# Load samples ------------------------------------------------------------
setwd(WD)
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


jscs <- runJunctionSeqAnalyses(
  sample.files = countFiles,
  sample.names = samples$sample.ID,
  condition = factor(samples$condition),
  flat.gff.file = GFF,
  nCores = cores,
  analysis.type = "junctionsAndExons"
)

dir.create(path =  paste0(WD, "/JunctionSeq/", comparison),
           recursive = TRUE)

buildAllPlots(
  jscs = jscs,
  outfile.prefix = paste(WD,"/JunctionSeq/", comparison, "/plots/", sep = ""),
  use.plotting.device = "png",
  FDR.threshold = 0.01
)
