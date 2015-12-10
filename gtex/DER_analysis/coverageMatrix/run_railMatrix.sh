#!/bin/sh

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis
WDIR=${MAINDIR}/coverageMatrix

# Define variables
CORES=40

mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs


for chrnum in Y X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
    do 
    sname="gtex-railMat-chr${chrnum}"
    ## Create script
    echo "Creating script for chromosome ${chrnum}"
    cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=5G,h_vmem=10G,h_fsize=100G
#$ -N ${sname}
#$ -pe local ${CORES}

echo '**** Job starts ****'
date

# Run railMatrix()
module load R/3.2.x
Rscript run_railMatrix.R -c "${chrnum}"

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
#$ -hold_jid gtex-railMat-chr*

echo '**** Job starts ****'
date

# Merge railMatrix() results
module load R/3.2.x
Rscript merge_railMatrix.R

## Move log files into the logs directory
mv ${sname}.* logs/

echo '**** Job ends ****'
date

EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call
