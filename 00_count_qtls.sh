#!/bin/bash

for trait in AUF AUG EIN ENT ERT_P FLE FORML FORMQ KNG RFE RHI SCH SCHF SCO STAE PPO SFL2 KON GES GRA CHI8 UEpop Npop ZW HH WR
do
	#awk -F , '$3 ~ "additive"' QTLs_sig_${trait}.csv | wc -l
	#awk -F , '$3 ~ "general"' QTLs_sig_${trait}.csv | wc -l
	
	echo "${trait}"
	awk -F , '$3 ~ "additive"' QTLs_sig_${trait}.csv | cut -d"," -f10 | sort -n | sed -n '$s/^//p'
	awk -F , '$3 ~ "general"' QTLs_sig_${trait}.csv | cut -d"," -f10 | sort -n | sed -n '$s/^//p'
done
