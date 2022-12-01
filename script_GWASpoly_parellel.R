library(GWASpoly)
library(ggplot2)
setwd("/1data/mariusW/PotatoTools_GWAS/")

args = commandArgs(trailingOnly=TRUE)

chrom = args[1]
trait = args[2]

genofile <- file(paste("inputfiles/geno_tetraploids_",chrom,".csv",sep=""))
phenofile <- file("inputfiles/AEM_table_reseq_renamed.csv")

#3) #26 trait in AEM file; genotype and phenotype sample names have to be equal
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile, format="numeric", n.traits=26, delim=",")

data.original <- set.K(data,LOCO=FALSE,n.core=1)

data.original.scan <- GWASpoly(data.original,models=c("additive","1-dom","2-dom"), traits=trait, n.core=1)

data2 <- set.threshold(data.original.scan,method="Bonferroni",level=0.05)

#get QTL
qtl <- get.QTL(data=data2,traits=trait,models=c("additive","1-dom","2-dom"),bp.window=5e6)
write.csv(qtl,file = paste("sig_QTLs_",trait,"_",chrom,".csv",sep=""), sep = "\t")