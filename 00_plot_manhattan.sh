#!/bin/bash

#additive model

for trait in AUF AUG EIN ENT ERT_P FLE FORML FORMQ KNG RFE RHI SCH SCHF SCO STAE PPO SFL2 KON GES GRA CHI8 UEpop Npop ZW HH WR
do
	echo "${trait}"
	awk '$7 > 2 && $7 != "NA" && $1 !~ /Unknown/' GWAS_${trait}_scores.csv | sed 's/Chr0//g' >GWAS_${trait}_scores_general_subset.csv
	sed -i 's/Chr//g' GWAS_${trait}_scores_general_subset.csv

	Rscript manhattan_plots.R ${trait}
	
done
