#!/bin/bash

#additive model

for trait in AUF AUG EIN ENT ERT_P FLE FORML FORMQ KNG RFE RHI SCH SCHF SCO STAE PPO SFL2 KON GES GRA CHI8 UEpop Npop ZW HH WR
do
	echo "${trait}"
	awk '$6 != "NA" && $1 !~ /Unknown/' GWAS_files_SV_dosage/GWAS_${trait}_scores.csv | sed 's/Chr0//g' >GWAS_files_SV_dosage/GWAS_${trait}_scores_additive_subset.csv
	sed -i 's/Chr//g' GWAS_files_SV_dosage/GWAS_${trait}_scores_additive_subset.csv

	Rscript manhattan_plots.R ${trait}
	
done


#dosage:
#add: 5.86    gen: 5.83