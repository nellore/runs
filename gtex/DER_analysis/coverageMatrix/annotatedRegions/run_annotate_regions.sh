#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=4G,h_vmem=5G,h_fsize=10G
#$ -N gtex-annoRegs

echo '**** Job starts ****'
date

## Make logs dir
mkdir -p logs

## Annotate regions
module load R/3.2.x
Rscript annotate_regions.R

## Move log files into the logs directory
mv gtex-annoRegs.* logs/

echo '**** Job ends ****'
date
