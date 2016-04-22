#!/bin/sh

# Usage
# sh run_rse.sh sra
# sh run_rse.sh gtex

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/recount2
WDIR=${MAINDIR}/rse

# Define variables
PROJECT=$1

# Create log dir
mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell file

if [[ "${PROJECT}" == "sra" ]]
then
    echo "$PROJECT"
elif [[ "${PROJECT}" == "gtex" ]]
then
    echo "$PROJECT"
else
    echo "Specify a valid project: gtex, sra"
fi

# Count how many commands there are
# For testing use: LINES=10
#LINES=$(wc -l ${WDIR}/project_ids_${PROJECT}.txt | cut -f1 -d " ")
LINES=10
METADATA="${MAINDIR}/metadata/metadata_${PROJECT}.Rdata"


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.rse"
    
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m a
#$ -l mem_free=10G,h_vmem=15G,h_fsize=20G
#$ -N ${sname}
#$ -t 1:${LINES}

PROJECTNAME=\$(awk "NR==\${SGE_TASK_ID}" ${WDIR}/bwtool_cmds_${PROJECT}.txt)

echo "**** Job starts project \${PROJECTNAME} ****"
date

## Run the R script that creates the count and rse files for exons and genes
module load R/3.3
Rscript create_rse.R -p "${PROJECT}" -m "${METADATA}" -i "\${PROJECTNAME}"

echo "**** Job ends ****"
date

## Move log files
mv ${WDIR}/${sname}.*.\${SGE_TASK_ID} ${WDIR}/logs/
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call


