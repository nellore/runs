#!/usr/bin/env bash
# Creates a table with sample metadata to facilitate manual construction of a tidier table.
# Requires curl, php, and index_to_sra_accession.tsv
for samp in $(cut -f3 index_to_SRA_accession.tsv | sort | uniq); do curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=biosample&id=$(curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term=$samp" | sed -n 's|[^<]*<Id>\([^<]*\)</Id>[^<]*|\1|gp')" | grep '<SampleData>' | php -r 'while(($line=fgets(STDIN)) !== FALSE) echo html_entity_decode($line, ENT_QUOTES|ENT_HTML401);' | awk -v mysamp=$samp '{print mysamp "\t" $0}'; done >biosample_data_to_parse.tsv
cat biosample_data_to_parse.tsv | python xmlparse.py >parsed_biosample_data.tsv
cat parsed_biosample_data.tsv | python tag.py >biosample_tags.tsv
#cat parsed_biosample_data.tsv | python classify_cell_line.py
#cat <(paste -d'\t' <(cat parsed_biosample_data.tsv | tr '[:upper:]' '[:lower:]' | grep 'line:' | awk -F 'line:' '{print $2}' | awk -F ';' '{print $1}') <(grep -i 'line:' parsed_biosample_data.tsv)) <(grep -iv 'line:' parsed_biosample_data.tsv | awk '{print "NA\t" $0}') | awk '$1=="--" {print "NA\t" substr($0,4)} $1!="--" {print $0}' >parsed_biosample_data_cell_line_after_colon.tsv