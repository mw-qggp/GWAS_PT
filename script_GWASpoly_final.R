library(GWASpoly)
library(ggplot2)
setwd("/1data/mariusW/PotatoTools_GWAS/")

args = commandArgs(trailingOnly=TRUE)

trait = args[1]

#testing script 11.01.2023
load("kmatrix.R")

#testing script 11.01.2023
params <- set.params(n.PC = 2)

data.original.scan <- GWASpoly(data.original,models=c("additive", "general"), traits=trait, n.core=10,params=params) #additive and general

#8) Filter - threshold
data2 <- set.threshold(data.original.scan,method="Bonferroni",level=0.05)

#get QTL
qtl <- get.QTL(data=data2,traits=trait,models=c("additive","general"),bp.window=5e6)
write.csv(qtl,file = paste("QTLs_sig_",trait,".csv",sep = ""), sep = "\t")

#write output GWAS
write.GWASpoly(data=data2, trait=trait, filename= paste("GWAS_",trait,"_effects.csv", sep = ""), what = "effects", delim = "\t")
write.GWASpoly(data=data2, trait=trait, filename= paste("GWAS_",trait,"_scores.csv",sep = ""), what = "scores", delim = "\t")
