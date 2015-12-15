#!/usr/bin/env bash
# Writes commands that should be executed to generate coverage sum files by sample
# $1: path to bwtool 1.0; download it from https://github.com/CRG-Barcelona/bwtool/releases/tag/1.0
# $2: path to BED file with regions
# $3: path to GTEx output directory; assumes directory structure of output directory in download.sh
# $4: path to SraRunInfo.csv from GTEx repository
# $5: path to directory in which to dump summed coverages
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BWTOOL=$1
BED=$2
export GTEX=$3
RUNINFO=$4
DUMP=$5
for bw in $(find $GTEX -name "SRR*.bw" | grep -v unique)
do
	echo "sh ${DIR}/sum.sh $BWTOOL $BED $bw $RUNINFO $DUMP"
done