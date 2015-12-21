#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=20G,h_vmem=30G,h_fsize=10G
#$ -N gtex-bp-annoRegs
#$ -hold_jid gtex-annoRegs

echo '**** Job starts ****'
date

## Make logs dir
mkdir -p logs

## Annotate regions
module load R/3.2.x
Rscript annotated_bp.R

## Move log files into the logs directory
mv gtex-bp-annoRegs.* logs/

echo '**** Job ends ****'
date
