library(GWASpoly)
library(ggplot2)
setwd("/1data/mariusW/PotatoTools_GWAS/")

#VCF2dosage(VCF.file, dosage.file, geno.code, ploidy, samples=NULL, min.DP=1, max.missing, min.minor=5) - example/documentation
#vcf = inputfile, dosage = outputfile, gencode = GT, ploidy = 4 etc.
#1)
#VCF2dosage("Reheader_filtered_tetraploids.vcf", "geno_tetraploids.csv", "GT", ploidy = 4, samples=NULL, min.DP=10, max.missing = 0.05, min.minor=2)
?set.params
?set.K
?GWASpoly
?read.GWASpoly
?write.GWASpoly

#2)
genofile <- file("inputfiles/geno_tetraploids_Chr02.csv") #chr02 for testing script (11.01.2023) #K matrix all chromosomes
phenofile <- file("inputfiles/AEM_table_reseq_renamed.csv")

#3) #26 trait in AEM file; genotype and phenotype sample names have to be equal
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile, format="numeric", n.traits=26, delim=",")

#4) Population structure
#data.loco <- set.K(data,LOCO=TRUE,n.core=10)

data.original <- set.K(data,LOCO=FALSE,n.core=10)

#To use it in GWAS script for each chromosome
save(data.original, file = "kmatrix.R")

#testing script 11.01.2023
load("kmatrix.R")

#PCA matrix?

#5) Marker curation
#N <- 100 #Population size
#params <- set.params(geno.freq = 1 - 5/N, fixed = "env", fixed.type = "factor", n.PC = 2) #Parameters: PC1 + PC2 as fixed effect (SNP matrix)

#testing script 11.01.2023
params <- set.params(n.PC = 2)

#6)
#without params=params; no LOCO

data.original.scan <- GWASpoly(data.original,models=c("additive", "general"), traits=c("AUF"), n.core=10,params=params) #additive and general

#7)

#8) Filter - threshold
data2 <- set.threshold(data.original.scan,method="Bonferroni",level=0.05)
#data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)

#9) Manhattan plot
#p <- manhattan.plot(data2,traits="AUF", models = "additive")
#p + theme(axis.text.x = element_text(angle=90,vjust=0.5))

p <- manhattan.plot(data2,traits="AUF",chrom="Chr01")
ggsave("test.png", p, bg = "transparent") #Found somewhere

#p <- LD.plot(data2, max.loci=1000)

#get QTL
qtl <- get.QTL(data=data2,traits="AUF",models=c("additive","general"),bp.window=5e6)
write.csv(qtl,file = "sig_qtls_AUF.csv", sep = "\t")

#write output GWAS
write.GWASpoly(data=data2, trait="AUF", filename="test_GWAS_AUF_effects.csv", what = "effects", delim = "\t")
write.GWASpoly(data=data2, trait="AUF", filename="test_GWAS_AUF_scores.csv", what = "scores", delim = "\t")

?write.GWASpoly

#???
#knitr::kable(qtl)
#fit.ans <- fit.QTL(data=data2,trait="AUF",qtl=qtl[,c("Marker","Model")],fixed=data.frame(Effect="env",Type="factor"))
#knitr::kable(fit.ans,digits=3)
