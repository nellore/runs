#!/bin/sh

# Directories
MAINDIR=/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix
WDIR=${MAINDIR}/bwtool


mkdir -p ${WDIR}
mkdir -p ${WDIR}/logs

# Construct shell files
sh /dcl01/leek/data/gtex_work/runs/gtex/generate_means.sh /dcl01/leek/data/bwtool/bwtool-1.0/bwtool /dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/region_bed/regions-cut0.5.bed /dcl01/leek/data/gtex /dcl01/leek/data/gtex_work/runs/gtex/SraRunInfo.csv /dcl01/leek/data/gtex_work/mean_cov | while read bwtoolcmd
    do
    bwfile=$(echo "${bwtoolcmd}" | cut -f5 -d " ")
    bwsample=$(basename ${bwfile} .bw)
    echo "Creating script for ${bwsample}"
    sname="${bwsample}.bwtool"
    
	cat > ${WDIR}/.${sname}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -m e
#$ -l mem_free=1G,h_vmem=2G,h_fsize=100G
#$ -N ${sname}
#$ -hold_jid region-bed
echo "**** Job starts ****"
date

## Run bwtool
${bwtoolcmd}

mv ${WDIR}/${sname}.* ${WDIR}/logs/

echo "**** Job ends ****"
date
EOF

	call="qsub ${WDIR}/.${sname}.sh"
	echo $call
	$call
done
