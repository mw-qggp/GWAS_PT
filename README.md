### GWAS_PT

##### 28.11.2022 - Convert SNPs/Indels in vcf file to inputfile for GWASpoly

```
Rscript convert_vcf.R
```

##### 01.12.2022 - 02.12.2022

Subset marker inputfile by chromosome and let GWASpoly script run for each trait:

```
for chrom in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12
do
        awk -F "," -v chr="${chrom}" '$1 ~ /Marker/ || $2 == chr' inputfiles/geno_tetraploids.csv >inputfiles/geno_tetraploids_${chrom}.csv

        for trait in "AUF" "AUG" "EIN" "ENT" "ERT_P" "FLE" "FORML" "FORMQ" "KNG" "RFE" "RHI" "SCH" "SCHF" "SCO" "STAE" "PPO" "SFL2" "KON" "GES" "GRA" "CHI8" "UEpop" "Npop" "ZW" "HH" "WR"
        do
                Rscript script_GWASpoly_parellel.R ${chrom} ${trait}
        done
done
```

Concatenate all significant markers of all traits sorted by trait, chromosome, and position:

```
cat sig_QTLs_* | sort -t "," -k2,2 -k6,6 -k7,7n | uniq | tac | awk 'NR==1 {line =$0; next} 1; END{print line}' | tac | sed '2d;s/,/\t/g' >QTLs_10traits_testrun.csv
```

##### 02.02.2023 - check significant SNPs and annotation of genes

```
cat Agria_PASA_AGAT.LC.gff3 Agria_PASA_AGAT.HC.gff3 | awk '$1 !~ /#/ && $3 == "gene"' | sort -k1,1 -k4,4n -k5,5n >Agria_PASA_AGAT.all.genesRow.gff3

cat Agria_PASA_AGAT.LC.gff3 Agria_PASA_AGAT.HC.gff3 | awk '$1 !~ /#/ && $3 == "intron"' | sort -k1,1 -k4,4n -k5,5n >Agria_PASA_AGAT.all.intronsRow.gff3
```

Next step: bash run2.sh:

```
#!/bin/bash

for trait in AUF AUG EIN ENT ERT_P FLE FORML FORMQ KNG RFE RHI SCH SCHF SCO STAE PPO SFL2 KON GES GRA CHI8 UEpop Npop ZW HH WR
do
        Rscript script_GWASpoly_final.R ${trait}

        awk '$6 >= 8.5 && $6 != "NA" {print $2"\t"$3"\t"$3"\t"$6}' GWAS_files/GWAS_${trait}_scores.csv >inputfiles_sigSNPs/input_SNPsGenes_${trait}_additive_sig.csv

	for chrom in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chromosome_Unknown
        do
		python3 00_SNPs_genes_overlap.py inputfiles_sigSNPs/input_SNPsGenes_${trait}_additive_sig.csv ${chrom} SNPs
        done
done
```

##### SV, SV dosage (manta), manhattan plots etc. - 21.02.2023 and following days

```
for trait in AUF AUG EIN ENT ERT_P FLE FORML FORMQ KNG RFE RHI SCH SCHF SCO STAE PPO SFL2 KON GES GRA CHI8 UEpop Npop ZW HH WR
do
        Rscript script_GWASpoly_finalSV.R ${trait}

        awk '$6 >= 5.86 && $6 != "NA" {print $2"\t"$3"\t"$3"\t"$6}' GWAS_files_SV_dosage/GWAS_${trait}_scores.csv >inputfiles_sigSV_dosage/input_SVGenes_${trait}_additive_sig.csv

        for chrom in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chromosome_Unknown
        do
                python3 00_SNPs_genes_overlap.py inputfiles_sigSV_dosage/input_SVGenes_${trait}_additive_sig.csv ${chrom} SV_dosage ${trait}
        done
done


for trait in AUF AUG EIN ENT ERT_P FLE FORML FORMQ KNG RFE RHI SCH SCHF SCO STAE PPO SFL2 KON GES GRA CHI8 UEpop Npop ZW HH WR
do
        echo "${trait}"
        awk '$6 != "NA" && $1 !~ /Unknown/' GWAS_files_SV_dosage/GWAS_${trait}_scores.csv | sed 's/Chr0//g' >GWAS_files_SV_dosage/GWAS_${trait}_scores_additive_subset.csv
        sed -i 's/Chr//g' GWAS_files_SV_dosage/GWAS_${trait}_scores_additive_subset.csv

        Rscript manhattan_plots.R ${trait}
        
done
```

##### check orthologs

```
#!/bin/bash

for file in $(ls SV_dosage/)
do
        python3 01_SNPs_genes_overlap.py SV_dosage/${file} >SV_dosage/${file}_orthologs.txt
done
```


##### Preparation of SV data on HPC

Merge VCFs (for 0/1 matrix, without dosage):

```
/gpfs/project/projects/qggp/src/SURVIVOR/bin/SURVIVOR merge list_manta_vcfs.txt 500 1 1 1 0 0 manta_merged.vcf
```

Create 0/1 matrix:

```
awk '$1 ~ /CHROM/' manta_merged.vcf >header_manta_merged.txt
awk '{split($8,a,";");split(a[2],b,"=")} $1 !~ /#/ {print b[2]}' manta_merged.vcf | sed 's/0/0\t/g;s/1/1\t/g' | sed -z 's/\t\n/\n/g' >matrix_0_1_manta.csv
awk '{split($8,a,";");split(a[2],b,"=")} $1 !~ /#/ {print b[2]}' manta_merged.vcf | sed 's/0/0\t/g;s/1/1\t/g' | sed -z 's/\t\n/\n/g' >matrix_0_1_manta.csv
awk '{split($8,a,";");split(a[2],b,"=")} $1 !~ /#/ {print $1"\t"$2}' manta_merged.vcf >chrom_pos_manta_merged.txt

nano header_manta_merged.txt #add chrom and pos at header
paste chrom_pos_manta_merged.txt matrix_0_1_manta.csv >matrix_0_1_manta_pos.csv
cat header_manta_merged.txt matrix_0_1_manta_pos.csv >matrix_0_1_manta_header.csv

rm matrix_0_1_manta_pos.csv matrix_0_1_manta.csv chrom_pos_manta_merged.txt

awk '{print $1"_"$2"\t"$1"\t"$2"\tA\tT\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31"\t"$32"\t"$33"\t"$34"\t"$35"\t"$36"\t"$38"\t"$40"\t"$41"\t"$42"\t"$43"\t"$44"\t"$45"\t"$46"\t"$47"\t"$48"\t"$49"\t"$50"\t"$51"\t"$52"\t"$53"\t"$54"\t"$55"\t"$56"\t"$57"\t"$58"\t"$59"\t"$60"\t"$61"\t"$62"\t"$63"\t"$64"\t"$65"\t"$66"\t"$69"\t"$70"\t"$71"\t"$72"\t"$73"\t"$74"\t"$75"\t"$76"\t"$77"\t"$78"\t"$79"\t"$80"\t"$81"\t"$82"\t"$83"\t"$84"\t"$85"\t"$86"\t"$87"\t"$88"\t"$89"\t"$90"\t"$93"\t"$94"\t"$95"\t"$96"\t"$97"\t"$98"\t"$99"\t"$100"\t"$101"\t"$102"\t"$103"\t"$104"\t"$105"\t"$106"\t"$107"\t"$108"\t"$109"\t"$110}' matrix_0_1_manta_header.csv | sed 's/\t/,/g' >matrix_0_1_manta_header_100genotypes_GWAS.csv
```

Dosage manta VCF files:

```
for sample in Allians Adretta Agria Albtras Allians Altus Ambition Atlntic Atzimba belena Bintje BNA_1 BNA_2 BNA_3 BNA_4 BNA5 Cara Celtane CGN17881 CGN_18114 Charlte Cherie Colomba Dark Desiree Donata EaryRse Edison Eurogrande Europrim Felsina Flava Fontane Gala Gladitor GLKS3 Granola H98D12 Harpun hermes Innovat Jelly jukijiro Karelia karolin Kathadin KENN KingRust Kolibri Krone Kuba Kuras LadyRose Laura Leyla Lilly Marabel MarsPipr Natalia Nevsky Nicola Odysseus Olympus ONA Otolia P3 P40 PentDell Pirol Princess p_russet Quadriga Quarta record Regina Rode_Est Rooster Rosara Rudolph russet_burbank s13_017 S14214 S14317_1 S14_9122 S5061 S5089 S5202 S5220a S5342 SA222_2 SA225_2 saskia Semlo seresta shc909 Shepody Skawa snowden solist Spunta Talent tyoshi Udacha Velox verdi Vitablla vr808 Yangana zorba
do
        awk '$7 == "PASS"' manta_${sample}/results/variants/diploidSV_INV.vcf >manta_${sample}/results/variants/diploidSV_INV_PASS.vcf
        awk '$1 !~ /#/ && $9 ~ /SR/ {split($10,a,":");split(a[6],b,","); split(a[5], c, ",");  ref = (b[1]+c[1])/(b[1]+b[2]+c[1]+c[2]); alt = (b[2]+c[2])/(b[1]+b[2]+c[1]+c[2]); print $0"\t"ref"\t"alt} $9 !~ /SR/ && $9 ~ /PR/ {split($10,a,":"); split(a[5],b,","); ref = b[1]/(b[1]+b[2]); alt = b[2]/(b[1]+b[2]); print $0"\t"ref"\t"alt} $1 ~ /#/ {print}' manta_${sample}/results/variants/diploidSV_INV_PASS.vcf >manta_${sample}/results/variants/diploidSV_INV_dosage1.vcf

        awk '$12 == 0 && $1 !~ /#/ {print $0"\t"0} $12 > 0 && $12 <= 0.25 && $1 !~ /#/ {print $0"\t"1} $12 > 0.25 && $12 <= 0.5 && $1 !~ /#/ {print $0"\t"2} $12 > 0.5 && $12 <= 0.75 && $1 !~ /#/ {print $0"\t"3} $12 > 0.75 && $1 !~ /#/ {print $0"\t"4} $1 ~ /#/ {print}' manta_${sample}/results/variants/diploidSV_INV_dosage1.vcf >manta_${sample}/results/variants/diploidSV_INV_dosage2.vcf

        rm manta_${sample}/results/variants/diploidSV_INV_dosage1.vcf
done
```

Merge vcf files:
```
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
```

After this, same awk command as above for 0/1 matrix ...