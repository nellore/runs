#!/bin/sh

# Usage
# sh run_fileinfo.sh sra
# sh run_fileinfo.sh gtex

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/recount2
WDIR=${MAINDIR}/fileinfo

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
#LINES=10
LINES=$(wc -l ${MAINDIR}/metadata/project_ids_${PROJECT}.txt | cut -f1 -d " ")
METADATA="${MAINDIR}/metadata/metadata_${PROJECT}.Rdata"


echo "Creating script for project ${PROJECT}"
sname="${PROJECT}.fileinfo"
    
cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m a
#$ -l mem_free=1G,h_vmem=2G,h_fsize=10G
#$ -N ${sname}
#$ -t 1:${LINES}
#$ -hold_jid ${PROJECT}.mean

PROJECTNAME=\$(awk "NR==\${SGE_TASK_ID}" ${MAINDIR}/metadata/project_ids_${PROJECT}.txt)

echo "**** Job starts project \${PROJECTNAME} ****"
date

## Run the R script that computes md5sum and file sizes as well as the final
## list of files to upload
module load R/3.3
Rscript get_fileinfo.R -p "${PROJECT}" -m "${METADATA}" -i "\${PROJECTNAME}"

echo "**** Job ends ****"
date

## Move log files
mv ${WDIR}/${sname}.*.\${SGE_TASK_ID} ${WDIR}/logs/
EOF

call="qsub ${WDIR}/.${sname}.sh"
echo $call
$call


