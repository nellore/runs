#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=15G,h_vmem=20G,h_fsize=30G
#$ -N test10-compare

echo '**** Job starts ****'
date


## Make logs dir
mkdir -p logs

# Format the pheno table
module load R/3.2.x
Rscript compare_results.R

## Move log files into the logs directory
mv test10-compare.* logs/

echo '**** Job ends ****'
date
