#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=2G,h_vmem=3G
#$ -pe local 25
#$ -N ers_gtex-check-tsv
#$ -hold_jid sra_ers_gtex_bwtool

echo '**** Job starts ****'
date

## Make logs dir
mkdir -p logs

## Annotate regions
module load R/3.3
Rscript check_tsv.R

## Move log files into the logs directory
mv ers_gtex-check-tsv.* logs/

echo '**** Job ends ****'
date
