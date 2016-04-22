#!/bin/sh

# Usage
# sh run_bwtool.sh sra
# sh run_bwtool gtex

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/recount2
WDIR=${MAINDIR}/bwtool

# Define variables
PROJECT=$1

# Create log dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell file
# For testing use: "head -n 10 |" before > ${WDIR}/bwtool_cmds_${PROJECT}.txt

if [[ "${PROJECT}" == "sra" ]]
then
    sh /dcl01/leek/data/gtex_work/runs/recount2/generate_sums.sh /dcl01/leek/data/bwtool/bwtool-1.0/bwtool /dcl01/leek/data/gtex_work/runs/recount2/genes/ucsc-knowngene-hg38.bed /dcl01/leek/data/sra/v2 /dcl01/leek/data/recount2/coverage > ${WDIR}/bwtool_cmds_${PROJECT}.txt
elif [[ "${PROJECT}" == "gtex" ]]
then
    sh /dcl01/leek/data/gtex_work/runs/recount2/generate_sums.sh /dcl01/leek/data/bwtool/bwtool-1.0/bwtool /dcl01/leek/data/gtex_work/runs/recount2/genes/ucsc-knowngene-hg38.bed /dcl01/leek/data/gtex /dcl01/leek/data/recount2/coverage > ${WDIR}/bwtool_cmds_${PROJECT}.txt
else
    echo "Specify a valid project: gtex, sra"
fi

# Count how many commands there are
LINES=$(wc -l ${WDIR}/bwtool_cmds_${PROJECT}.txt | cut -f1 -d " ")


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.bwtool"
    
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m a
#$ -l mem_free=1G,h_vmem=2G,h_fsize=100G
#$ -N ${sname}
#$ -t 1:${LINES}
#$ -hold_jid genes-hg38-ucsc

## Get the bwtool command
bwtoolcmd=\$(awk "NR==\${SGE_TASK_ID}" ${WDIR}/bwtool_cmds_${PROJECT}.txt)

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
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call


