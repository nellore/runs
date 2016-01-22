#!/usr/bin/env bash
# Helper script that uses bwtool to generate coverage sums in regions from bed file
# $1: path to bwtool 1.0; download it from https://github.com/CRG-Barcelona/bwtool/releases/tag/1.0
# $2: path to BED file with regions
# $3: path to bigWig file in GTEx output directory; assumes directory structure of output directory in download.sh
# $4: path to AUC file
# $5: path to directory in which to dump coverage sums
BWTOOL=$1
BED=$2
export BW=$3
export AUC=$4
DUMP=$5
mkdir -p $DUMP
export SAMPLENAME=$(echo $(basename $BW) | rev | awk '{print substr($0, 4)}' | rev)
SRR=$(echo $SAMPLENAME | cut -d'_' -f1)
export COUNTS=$(dirname $(dirname $BW))/cross_sample_results/counts.tsv.gz
CURRENTAUC=$(grep $SRR $AUC | cut -f2)
$BWTOOL summary $BED $BW /dev/stdout -fill=0 -with-sum | cut -f1-3,10 | awk -v currauc=$CURRENTAUC -v CONVFMT=%.17g '{print $1 "\t" $2 "\t" $3 "\t" $4*4000000000/currauc}' >$DUMP/$SRR.sum.tsv