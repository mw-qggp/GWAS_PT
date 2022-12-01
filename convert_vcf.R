library(GWASpoly)
setwd("/1data/mariusW/PotatoTools_GWAS/")

VCF2dosage("Reheader_filtered_tetraploids.vcf", "geno_tetraploids.csv", "GT", ploidy = 4, samples=NULL, min.DP=10, max.missing = 0.05, min.minor=2)
