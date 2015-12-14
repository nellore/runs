#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=100G,h_vmem=120G,h_fsize=30G
#$ -N region-bed

echo '**** Job starts ****'
date


## Make logs dir
mkdir -p logs

# Find region cuts
module load R/3.2.x
Rscript find_region_bed.R

## Move log files into the logs directory
mv region-bed.* logs/

echo '**** Job ends ****'
date
