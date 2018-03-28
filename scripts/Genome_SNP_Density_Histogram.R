# Inputs ------------------------------------------------------------------
#input is .lqual file from vcftools (--site-quality option)
snp_quality <- ""
genotype <- ""

# main --------------------------------------------------------------------
library(ggplot2)
snps <-
  read.table(snp_quality,
    sep = "\t",
    header = TRUE
  )
#sort chromosomes and remove patches
goodchrorder <- c(1:22, "X", "Y", "MT")
snps$CHROM <- factor(snps$CHROM, levels = goodchrorder)
snps <- subset(snps, snps$CHROM != "NA")

snpdensity <- ggplot(snps) +
  geom_histogram(aes(x = POS), binwidth = 1e6) +
  facet_wrap( ~ CHROM, ncol = 2, strip.position = "left") +
  ggtitle(paste("Density of SNPs in", genotype)) +
  ylim(0, 1500) +
  xlab("Position") +
  ylab("SNP Density") +
  theme_bw()

snpdensity
pdf(file = paste(genotype, "_SNP_Density", sep = ""))
snpdensity
dev.off()