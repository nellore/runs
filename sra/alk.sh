#!/usr/bin/env bash
# Compressed junctions file (intropolis.v1.tsv.gz) should be first command-line parameter
# According to http://grch37.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000171094;r=2:29415640-30144432;t=ENST00000389048
# the ALK gene is on "Chromosome 2: 29,415,640-30,144,432 reverse strand."
JUNC=$1
# Select junctions only in ALK
gzip -cd $JUNC | grep -w chr2 | grep -w "-" | awk '$2 >= 29415640 && $2 <= 30144432 && $3 >= 29415640 && $3 <= 30144432' | gzip >alk_junctions.tsv.gz
# Now according to the Nature paper http://www.nature.com/nature/journal/v526/n7573/full/nature15258.html,
# many cancers have an alternative transcription initiation site after intron 19; in other words,
# for the so-called alternative ALK^{ATI} transcript, no junction with coordinates >= 29,446,395
# should be expressed. Divide these up into "start" (unexpressed in ATI; exons 1-19) expression and
# "end" (expressed in ATI; exons 20-29) expression.
gzip -cd alk_junctions.tsv.gz | awk '$2 >= 29446395' | gzip >alk_start_junctions.tsv.gz
gzip -cd alk_junctions.tsv.gz | awk '$2 < 29446395' | gzip >alk_end_junctions.tsv.gz