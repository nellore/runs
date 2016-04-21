#!/usr/bin/env bash
# Writes commands that should be executed to generate coverage sum files by sample
# $1: path to bwtool 1.0; download it from https://github.com/CRG-Barcelona/bwtool/releases/tag/1.0
# $2: path to BED file with regions
# $3: path to GTEx/SRA output directory; assumes directory structure of output directory in download.sh
# $4: path to directory in which to dump summed coverages
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BWTOOL=$1
BED=$2
export PROJECT=$3
DUMP=$4
for bw in $(find $PROJECT -name "*.bw" | grep -v -E "\\.[ACGNT]\\.bw" | grep -v -E "unique|mean|median" )
do
	echo "sh ${DIR}/sum.sh $BWTOOL $BED $bw $DUMP"
done