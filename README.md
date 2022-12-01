### GWAS_PT

Subset marker inputfile by chromosome (test):

```awk '$1 ~ /Marker/ || $1 ~ /Chr01/' inputfiles/geno_tetraploids.csv >inputfiles/geno_tetraploids_Chr01.csv```
