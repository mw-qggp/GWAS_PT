#!/bin/bash
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -l walltime=71:59:00
#PBS -A "denap"

cd $PBS_O_WORKDIR
#for sample in tyoshi SA225_2 S5220a P3 Olympus Nevsky Gladitor


for sample in Allians Adretta Agria Albtras Allians Altus Ambition Atlntic Atzimba belena Bintje BNA_1 BNA_2 BNA_3 BNA_4 BNA5 Cara Celtane CGN17881 CGN_18114 Charlte Cherie Colomba Dark Desiree Donata EaryRse Edison Eurogrande Europrim Felsina Flava Fontane Gala Gladitor GLKS3 Granola H98D12 Harpun hermes Innovat Jelly jukijiro Karelia karolin Kathadin KENN KingRust Kolibri Krone Kuba Kuras LadyRose Laura Leyla Lilly Marabel MarsPipr Natalia Nevsky Nicola Odysseus Olympus ONA Otolia P3 P40 PentDell Pirol Princess p_russet Quadriga Quarta record Regina Rode_Est Rooster Rosara Rudolph russet_burbank s13_017 S14214 S14317_1 S14_9122 S5061 S5089 S5202 S5220a S5342 SA222_2 SA225_2 saskia Semlo seresta shc909 Shepody Skawa snowden solist Spunta Talent tyoshi Udacha Velox verdi Vitablla vr808 Yangana zorba
do
	awk '$7 == "PASS"' manta_${sample}/results/variants/diploidSV_INV.vcf >manta_${sample}/results/variants/diploidSV_INV_PASS.vcf
	awk '$1 !~ /#/ && $9 ~ /SR/ {split($10,a,":");split(a[6],b,","); split(a[5], c, ",");  ref = (b[1]+c[1])/(b[1]+b[2]+c[1]+c[2]); alt = (b[2]+c[2])/(b[1]+b[2]+c[1]+c[2]); print $0"\t"ref"\t"alt} $9 !~ /SR/ && $9 ~ /PR/ {split($10,a,":"); split(a[5],b,","); ref = b[1]/(b[1]+b[2]); alt = b[2]/(b[1]+b[2]); print $0"\t"ref"\t"alt} $1 ~ /#/ {print}' manta_${sample}/results/variants/diploidSV_INV_PASS.vcf >manta_${sample}/results/variants/diploidSV_INV_dosage1.vcf

	awk '$12 == 0 && $1 !~ /#/ {print $0"\t"0} $12 > 0 && $12 <= 0.25 && $1 !~ /#/ {print $0"\t"1} $12 > 0.25 && $12 <= 0.5 && $1 !~ /#/ {print $0"\t"2} $12 > 0.5 && $12 <= 0.75 && $1 !~ /#/ {print $0"\t"3} $12 > 0.75 && $1 !~ /#/ {print $0"\t"4} $1 ~ /#/ {print}' manta_${sample}/results/variants/diploidSV_INV_dosage1.vcf >manta_${sample}/results/variants/diploidSV_INV_dosage2.vcf

	rm manta_${sample}/results/variants/diploidSV_INV_dosage1.vcf
done
