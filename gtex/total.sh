#!/usr/bin/env bash
# Computes total number of reads aligned from GTEx output data
# $1: GTEx output directory from download.sh
GTEX=$1
for i in $GTEX/*; do gzip -cd $i/cross_sample_results/counts.tsv.gz | tail -n +2 | rev | cut -f1 | rev | cut -d "," -f1 | paste -s -d+ | bc; done | paste -s -d+ | bc
