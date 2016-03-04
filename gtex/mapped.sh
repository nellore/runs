#!/usr/bin/env bash
# Outputs proportion of mapped reads in each sample across SRA
# $1: directory containing all the batch_{k} scripts
# Results of running this script are in mapped_prop.tsv
cd $1
paste <(for i in *; do gzip -cd $i/cross_sample_results/counts.tsv.gz | tail -n +2 | rev | cut -f2 | rev | cut -d',' -f1; done) <(for i in *; do gzip -cd $i/cross_sample_results/counts.tsv.gz | tail -n +2 | rev | cut -f1 | rev | cut -d',' -f1; done) | tr '\t' '/' | bc -l 
