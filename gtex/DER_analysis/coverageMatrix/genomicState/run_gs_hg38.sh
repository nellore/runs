#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=50G,h_vmem=60G,h_fsize=100G
#$ -N gs-hg38

echo '**** Job starts ****'
date


## Make logs dir
mkdir -p logs

## Create genomic state for hg38
module load R/3.2.x
Rscript hg38-genomicState.R

## Move log files into the logs directory
mv gs-hg38.* logs/

echo '**** Job ends ****'
date
