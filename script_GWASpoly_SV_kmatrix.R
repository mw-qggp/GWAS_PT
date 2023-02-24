library(GWASpoly)
library(ggplot2)
setwd("/1data/mariusW/PotatoTools_GWAS/")

#2)
genofile <- file("inputfiles/matrix_manta_header_dosage_100genotypes_GWAS.csv")
phenofile <- file("inputfiles/AEM_table_reseq_renamed.csv")

#3) #26 trait in AEM file; genotype and phenotype sample names have to be equal
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile, format="numeric", n.traits=26, delim=",")

#4) Population structure
#data.loco <- set.K(data,LOCO=TRUE,n.core=10)

data.original <- set.K(data,LOCO=FALSE,n.core=10)

#To use it in GWAS script for each chromosome
save(data.original, file = "kmatrix_SV_dosage.R")
