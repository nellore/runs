#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=15G,h_vmem=20G,h_fsize=30G
#$ -N genes-hg38-ucsc

echo '**** Job starts ****'
date


## Make logs dir
mkdir -p logs

# Create BED file
module load R/3.3
Rscript create_bed.R

## Move log files into the logs directory
mv genes-hg38-ucsc.* logs/

echo '**** Job ends ****'
date
