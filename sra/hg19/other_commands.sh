#!/usr/bin/env bash
# Miscellaneous commands used to analyze introns. Requires mawk 1.3.4 2015-05-03
cat all_SRA_introns.tsv | ./mawk-1.3.4-20150503/mawk -F ',' '{print $0 "\t" (NF - 1)/2 + 1}' | ./mawk-1.3.4-20150503/mawk '{print $1 "\t" $2 "\t" $3 "\t" $9}' >introns_and_sample_counts.tsv
cat all_SRA_introns.tsv | ./mawk-1.3.4-20150503/mawk '{print $1 "\t" $2 "\t" $3 "\t," $8}' | ./mawk-1.3.4-20150503/mawk -F ',' '{sum=0;for(i=2;i<=NF;i++){sum += $i} print $1 sum}' >introns_and_coverages.tsv
sort -k4,4nr -o introns_and_sample_counts.tsv introns_and_sample_counts.tsv
sort -k4,4nr -o introns_and_sample_counts.tsv introns_and_coverages.tsv
awk '{print "{ \"id\": " $1 ", \"project\": \"" $2 "\", \"sample\": \"" $3 "\", \"experiment\": \"" $4 "\", \"run\": \"" $5 "\" }"}' index_to_SRA_accession.tsv >index_to_SRA_accession.json
for i in $(cat index_to_SRA_accession.tsv | cut -f5); do curl "http://www.ncbi.nlm.nih.gov/gds?LinkName=sra_gds&from_uid=$(curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=$i" | sed -n 's|[^<]*<Id>\([^<]*\)</Id>[^<]*|\1|gp')" | grep Series | awk -F 'acc\=GSM' '{print "GSM" $2}' | grep -vFx GSM | awk -F '"' -v r=$i '{print r "\t" $1}'; done >SRRToGSM.tsv