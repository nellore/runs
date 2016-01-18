#!/usr/bin/env bash
# Compressed junctions file (intropolis.v1.tsv.gz) should be first command-line parameter
# According to http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;g=ENSG00000232810;r=6:31543344-31546113;t=ENST00000449264
# the TNF gene is on "Chromosome 6: 31,543,344-31,546,113 forward strand."
JUNC=$1
# Select junctions only in TNF
gzip -cd $JUNC | grep -w chr6 | grep -w "-" | awk '$2 >= 31543344 && $2 <= 31546113 && $3 >= 31543344 && $3 <= 31546113' | gzip >tnf_junctions.tsv.gz
