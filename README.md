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