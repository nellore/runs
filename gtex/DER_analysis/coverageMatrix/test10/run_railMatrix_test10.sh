#!/bin/sh

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix
WDIR=${MAINDIR}/test10

# Define variables
#CORES=4

mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs


#for chr in chrY chrX chr21 chr22 chr20 chr1 chr2 chr3
cut -f1 /dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes | while read chr 
    do 
    sname="test10-${chr}"
    ## Create script
    echo "Creating script for chromosome ${chr}"
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=15G,h_vmem=25G,h_fsize=100G
#$ -N ${sname}
#$ -hold_jid pheno-format

echo '**** Job starts ****'
date

# Run railMatrix()
module load R/3.2.x
Rscript run_railMatrix_test10.R -c "${chr}"

## Move log files into the logs directory
mv ${sname}.* logs/

echo '**** Job ends ****'
date

EOF

    call="qsub ${WDIR}/.${sname}.sh"
    echo $call
    $call
done


sname="gtex-railMat-merge"
## Create script
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=150G,h_vmem=400G,h_fsize=100G
#$ -N ${sname}
#$ -hold_jid test10-chr*

echo '**** Job starts ****'
date

# Merge railMatrix() results
module load R/3.2.x
Rscript merge_railMatrix_test10.R

## Move log files into the logs directory
mv ${sname}.* logs/

echo '**** Job ends ****'
date

EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call
