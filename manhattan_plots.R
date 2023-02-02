library(qqman)

setwd("/1data/mariusW/PotatoTools_GWAS")

args = commandArgs(trailingOnly=TRUE)
trait = args[1]

sig_qtls <- read.csv(paste("GWAS_",trait,"_scores_general_subset.csv", sep = ""), sep = "\t")

in_manhattan <- data.frame(sig_qtls$Marker,sig_qtls$om,sig_qtls$Position,sig_qtls$general)
colnames(in_manhattan) <- c("SNP","CHR","BP","P")

setEPS()
postscript(paste("manhattan_",trait,"general.eps", sep = ""))
manhattan(in_manhattan, logp = FALSE, genomewideline = 8.5, ylim = c(0,19), main = paste(trait," - general", sep = ""))
dev.off()

#awk + this script in for loop for traits ...
#additive
#awk '$6 != "NA" && $1 !~ /Unknown/ && $6 > 2' GWAS_org_file >GWAS_RHI_scores_additive_subset.csv
#sed  -i 's/Chr0//g' GWAS_RHI_scores_additive_subset.csv
#sed  -i 's/Chr//g' GWAS_RHI_scores_additive_subset.csv

#general:
#s.o.
#something wrong with general model - how to fix it?
#manhattan plots for traits additive model tomorrow. SNP in genes