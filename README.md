### GWAS_PT

##### 01.12.2022
Subset marker inputfile by chromosome (test):

```
#!/bin/bash

for chrom in Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12
do
        awk -F "," -v chr="${chrom}" '$1 ~ /Marker/ || $2 == chr' inputfiles/geno_tetraploids.csv >inputfiles/geno_tetraploids_${chrom}.csv
done
```

`bash run.sh`
