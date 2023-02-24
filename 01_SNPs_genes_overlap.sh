#!/bin/bash

for file in $(ls SV_dosage/)
do
	python3 01_SNPs_genes_overlap.py SV_dosage/${file} >SV_dosage/${file}_orthologs.txt
done
