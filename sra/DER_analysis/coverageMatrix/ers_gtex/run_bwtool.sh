#!/bin/sh

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/sra/DER_analysis/coverageMatrix
WDIR=${MAINDIR}/ers_gtex


mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell files
# For testing use: "head -n 10 |" before the "while read bwtoolcmd", without quotes
sh /dcl01/leek/data/gtex_work/runs/gtex/generate_sums.sh /dcl01/leek/data/bwtool/bwtool-1.0/bwtool /dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/region_bed/regions-cut0.5.bed /dcl01/leek/data/sra/v2 /dcl01/leek/data/gtex_work/runs/sra/v2/auc.tsv /dcl01/leek/data/sra_work/mean_cov_ers_gtex | head -n 10 > ${WDIR}/bwtool_cmds.txt


# Count how many commands there are
LINES=$(wc -l ${WDIR}/bwtool_cmds.txt | cut -f1 -d " ")

echo "Creating script for quantifying the SRA data with ERs defined by the GTEx data"
sname="sra_ers_gtex_bwtool"

cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m a
#$ -l mem_free=1G,h_vmem=2G,h_fsize=100G
#$ -N ${sname}
#$ -t 1:${LINES}
#$ -hold_jid region-bed

## Get the bwtool command
bwtoolcmd=\$(awk "NR==\${SGE_TASK_ID}" ${WDIR}/bwtool_cmds.txt)

## Extract the sample and print it
bwfile=\$(echo "\${bwtoolcmd}" | cut -f5 -d " ")
bwsample=\$(basename \${bwfile} .bw)

echo "**** Job starts sample \${bwsample} ****"
date

## Run bwtool
echo "\${bwtoolcmd}"
\${bwtoolcmd}

echo "**** Job ends ****"
date

## Move log files
mv ${WDIR}/${sname}.*.\${SGE_TASK_ID} ${WDIR}/logs/

date
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call
