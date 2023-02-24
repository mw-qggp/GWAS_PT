#!/bin/bash

for trait in AUF AUG EIN ENT ERT_P FLE FORML FORMQ KNG RFE RHI SCH SCHF SCO STAE PPO SFL2 KON GES GRA CHI8 UEpop Npop ZW HH WR
do
        Rscript script_GWASpoly_finalSV.R ${trait}

        awk '$6 >= 5.86 && $6 != "NA" {print $2"\t"$3"\t"$3"\t"$6}' GWAS_files_SV_dosage/GWAS_${trait}_scores.csv >inputfiles_sigSV_dosage/input_SVGenes_${trait}_additive_sig.csv

        for chrom in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chromosome_Unknown
        do
                python3 00_SNPs_genes_overlap.py inputfiles_sigSV_dosage/input_SVGenes_${trait}_additive_sig.csv ${chrom} SV_dosage ${trait}
        done
done

