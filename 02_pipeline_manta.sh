#!/bin/bash

module load Python/3.6.5

#mkdir manta_dosage/

>manta_dosage/genotyped_deletions_genome.csv

for sample in Adretta Agria Albtras Allians Altus Ambition Atlntic Atzimba belena Bintje BNA_1 BNA_2 BNA_3 BNA_4 BNA5 Cara Celtane CGN17881 CGN_18114 Charlte Cherie Colomba Dark Desiree Donata EaryRse Edison Eurogrande Europrim Felsina Flava Fontane Gala Gladitor GLKS3 Granola H98D12 Harpun hermes Innovat Jelly jukijiro Karelia karolin Kathadin KENN KingRust Kolibri Krone Kuba Kuras LadyRose Laura Leyla Lilly Marabel MarsPipr Natalia Nevsky Nicola Odysseus Olympus ONA Otolia P3 P40 PentDell Pirol Princess p_russet Quadriga Quarta record Regina Rode_Est Rooster Rosara Rudolph russet_burbank s13_017 S14214 S14317_1 S14_9122 S5061 S5089 S5202 S5220a S5342 SA222_2 SA225_2 saskia Semlo seresta shc909 Shepody Skawa snowden solist Spunta Talent tyoshi Udacha Velox verdi Vitablla vr808 Yangana zorba
do
	for chrom in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12
	do
		awk -v chr=${chrom} -v name=${sample} '$1 == chr && $3 !~ /BND/ {split($8,a,";");split(a[1],b,"="); print $1"\t"$2"\t"b[2]"\t"name"\t"$13}' manta_${sample}/results/variants/diploidSV_INV_dosage2.vcf >manta_dosage/${sample}_${chrom}_dosage.bed
	done
done

for chrom in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12
do
        cat manta_dosage/*_${chrom}_dosage.bed | sort -k2,2n -k3,3n >manta_dosage/genos_${chrom}_dosage.bed

        python3 mergeInbreds_newApproach_dels_dosage.py manta_dosage/genos_${chrom}_dosage.bed 500 ${chrom} manta_dosage/genotyped_deletions_${chrom}.csv

	awk '$1 !~ /Chrom/ {print}' manta_dosage/genotyped_deletions_${chrom}.csv >>manta_dosage/genotyped_deletions_genome_final.csv
done

sort -k1,1 -k2,2n -k3,3n manta_dosage/genotyped_deletions_genome_final.csv >manta_dosage/genotyped_deletions_genome_final_sorted.csv