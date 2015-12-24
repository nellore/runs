#!/usr/bin/env bash
# Computes AUC for each GTEx bigwig with wiggletools
# $1: path to wiggletools binary; we used v1.1, available at https://github.com/Ensembl/WiggleTools/releases/download/v1.1/wiggletools_x86_64_linux
# $2: path to output directory from download.sh containing batch_k subdirectories for k the batch number
# $3: output directory; final output dumped to auc.tsv with format 1st column: SRR number; second column: AUC
# $4: number of processes to run simultaneously
WIGGLETOOLS=$1
INPUT=$2
OUT=$3
PROCESSES=$4

function nrwait() {
    local nrwait_my_arg
    if [[ -z $1 ]] ; then
	nrwait_my_arg=2
    else
	nrwait_my_arg=$1
    fi
    
    while [[ $(jobs -p | wc -l) -ge $nrwait_my_arg ]] ; do
	sleep 0.33;
    done
  }

mkdir -p $OUT/tempaucs
cd $OUT/tempaucs
for bigwig in $(for i in $INPUT/*; do find $i/coverage_bigwigs -name \*.bw | grep -Ev 'mean|median|unique'; done)
do
	$WIGGLETOOLS AUC $bigwig >$(basename $bigwig | cut -d'_' -f1) &
	nrwait $PROCESSES
done
wait

for auc in *.auc
do
	echo -e $auc "\t" $(cat $auc) >>../auc.tsv
done

rm -rf $OUT/tempaucs
