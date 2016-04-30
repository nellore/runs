#!/bin/sh

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix
WDIR=${MAINDIR}/ers_gtex


mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

sname="sra_ers_gtex_bwtool-merge"
## Create script
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=10G,h_vmem=15G,h_fsize=100G
#$ -N ${sname}
#$ -pe local 25
#$ -hold_jid sra_ers_gtex_bwtool

echo '**** Job starts ****'
date

# Merge bwtool results
module load R/3.3
Rscript merge_bwtool.R

## Move log files into the logs directory
mv ${sname}.* logs/

echo '**** Job ends ****'
date

EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call
