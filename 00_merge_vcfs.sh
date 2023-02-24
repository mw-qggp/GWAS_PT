#!/bin/bash
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -l walltime=71:00:00
#PBS -A "denap"

#module load

cd $PBS_O_WORKDIR

/gpfs/project/projects/qggp/src/SURVIVOR/bin/SURVIVOR merge list_manta_vcfs.txt 500 1 1 1 0 0 manta_merged.vcf
