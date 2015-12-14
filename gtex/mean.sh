#!/usr/bin/env bash
# $1: path to bwtool 1.0; download it from https://github.com/CRG-Barcelona/bwtool/releases/tag/1.0
# $2: path to BED file with regions
# $3: path to bigWig file in GTEx output directory; assumes directory structure of output directory in download.sh
# $4: path to SraRunInfo.csv from GTEx repository
# $5: path to directory in which to dump mean coverage
BWTOOL=$1
BED=$2
export BW=$3
RUNINFO=$4
DUMP=$5
mkdir -p $DUMP
export SAMPLENAME=$(echo $(basename $BW) | rev | awk '{print substr($0, 4)}')
SRR=$(echo $SAMPLENAME | cut -d'_' -f1)
COUNTS=$(dirname $(dirname $BW))/cross_sample_results/counts.tsv.gz
READLENGTH=$(grep $SRR $RUNINFO | cut -d',' -f9)
READCOUNT=$(gzip -cd $COUNTS | grep $SAMPLENAME | rev | cut -f1 | rev | cut -d',' -f1)
$BWTOOL summary $BED $BW /dev/stdout | cut -f1-3,8 | awk -v rcount=$READCOUNT -v rlength=$READLENGTH '{print $1 "\t" $2 "\t" $3 "\t" $4/$rcount*40000000/$rlength*2}' >$DUMP/$SRR.mean.tsv