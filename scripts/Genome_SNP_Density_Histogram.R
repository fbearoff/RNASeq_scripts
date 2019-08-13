# Inputs ------------------------------------------------------------------
#input is .lqual file from bcftools (bcftools query -i'%FILTER="PASS"' -f'%CHROM %POS %QUAL\n')
WD <- ""
genotype <- "" #prefix of .lqual file
project <- "" #subdir of WD

# main --------------------------------------------------------------------
library(ggplot2)
snps <-
  read.table(file.path(WD, "SNPs", paste0(genotype, ".filtered.lqual")),
    header = FALSE
  )
#insert header
colnames(snps) <- c("CHROM", "POS", "QUAL")

#sort chromosomes and remove patches
goodchrorder <- c(1:22, "X", "Y", "MT")
snps$CHROM <- factor(snps$CHROM, levels = goodchrorder)
snps <- subset(snps, snps$CHROM != "NA")
#snps <- subset(snps, snps$CHROM == 10) #just one CHROM?

snpdensity <- ggplot(snps) +
  geom_histogram(aes(x = POS), binwidth = 1e6) +
  facet_wrap( ~ CHROM, ncol = 2, strip.position = "left") +
  ggtitle(paste("Density of SNPs in", genotype)) +
  ylim(0, 1500) +
  xlab("Position") +
  ylab("SNP Density") +
  theme_bw()

dir.create(file.path(WD, project, "SNPs"))

pdf(paste0(file.path(WD, project, "SNPs/"), genotype, "_SNP_Density.pdf"))
snpdensity
dev.off()
