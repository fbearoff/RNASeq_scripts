library(biomaRt)

snp_list <- readLines("snps.txt")

snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "mmusculus_snp")

snp_sequence <-
  getBM(
    attributes = c("refsnp_id", "snp", "chr_name", "chrom_start", "chrom_end", "allele"),
    filters = c("snp_filter", "upstream_flank", "downstream_flank"),
    checkFilters = FALSE,
    values = list(snp_list, 200, 200),
    mart = snpmart,
    bmHeader = TRUE
  )
snp_sequence
