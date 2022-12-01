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

#2)
genofile <- file("inputfiles/geno_tetraploids.csv")
phenofile <- file("inputfiles/AEM_Resequencing2.csv")

#3) #20 trait in AEM file; genotype and phenotype samle names have to be equal
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile, format="numeric", n.traits=20, delim=",")

#4) Population structure
data.loco <- set.K(data,LOCO=TRUE,n.core=8)
data.original <- set.K(data,LOCO=FALSE,n.core=8)

#5) Marker curation
N <- 100 #Population size
params <- set.params(geno.freq = 1 - 5/N, fixed = "env", fixed.type = "factor") #Ask BST

#6)
data.loco.scan <- GWASpoly(data=data.loco,models=c("additive"), traits=c("AUF"), n.core=8) #without params=params
data.original.scan <- GWASpoly(data.original,models=c("additive"), traits=c("AUF"), n.core=8) #additive or 1-dom

#7)

#8) Filter - threshold
data2 <- set.threshold(data.loco.scan,method="Bonferroni",level=0.05)
#data2 <- set.threshold(data.loco.scan,method="M.eff",level=0.05)

?manhattan.plot

#9) Manhattan plot
p <- manhattan.plot(data2,traits="AUF", models = "additive")
p + theme(axis.text.x = element_text(angle=90,vjust=0.5))

#manhattan.plot(data2,traits="AUF",chrom="Chr01")
#ggsave("DWSte.png", p, bg = "transparent") Found somewhere

#get QTL
qtl <- get.QTL(data=data2,traits="AUF",models="additive",bp.window=5e6)

#???
#knitr::kable(qtl)
#fit.ans <- fit.QTL(data=data2,trait="vine.maturity",qtl=qtl[,c("Marker","Model")],fixed=data.frame(Effect="env",Type="factor"))
#knitr::kable(fit.ans,digits=3)