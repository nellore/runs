#!/usr/bin/env bash
# Creates BED file with intropolis junctions
# $1: location of intropolis.v1.hg19.tsv.gz
# $2: destination directory
# $3: location of bedToBigBed executable (obtained from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed)
INTROPOLIS=$1
DEST=$2
BEDTOBIGBED=$3
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p $DEST
gzip -cd $INTROPOLIS | awk '{print $1 "\t" ($2 - 1) "\t" $3 "\tjunction_" NR "\t1000\t" $4}' | gzip >$DEST/intropolis.v1.hg19.bed.gz
cd $DEST
gzip -cd intropolis.v1.hg19.bed.gz >intropolis.v1.hg19.bed
$BEDTOBIGBED intropolis.v1.hg19.bed $DIR/hg19.chrom.sizes intropolis.v1.hg19.bb
