#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=15G,h_vmem=20G,h_fsize=200G
#$ -N gtex-railMat
#$ -pe local 24

echo '**** Job starts ****'
date

# Make logs directory
mkdir -p logs

# Run railMatrix()
module load R/devel
Rscript run_railMatrix.R

## Move log files into the logs directory
mv gtex-railMat.* logs/

echo '**** Job ends ****'
date
