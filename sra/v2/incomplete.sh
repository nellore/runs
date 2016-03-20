#!/usr/bin/env bash
# Computes stats on incompletely downloaded samples. First command-line parameter should be SRA download dir with batch_* subdirectories.
SRADIR=$1
pypy incomplete.py --sra-dir $SRADIR | (read -r; printf "%s\n" "$REPLY"; sort -k7,7nr) >incomplete.tsv
join -1 1 -2 1 -a 1 <(cat incomplete.tsv | tail -n +2 | grep -v NA | cut -f1 | less | sort | uniq -c | sed 's/^ *//g' | awk '{print $2 "\t" $1}') <(cut -f2 intropolis.idmap.v2.hg38.tsv | sort | uniq -c | sed 's/^ *//g' | awk '{print $2 "\t" $1}') | awk '{print $1 "\t" $2 "\t" $3 "\t" $2/$3}' | sort -k2,2nr >incomplete_projects.tsv
paste <(head -n 15 incomplete_projects.tsv | cut -f1,2) <(for i in $(head -n 15 incomplete_projects.tsv | cut -f1); do grep -w $i <(grep -v NA incomplete.tsv) | cut -f7 | awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++) { print sum[i]/NR "\t" sqrt((sumsq[i]-sum[i]^2/NR)/(NR-1))} }'; done) >top_incomplete_projects.tsv
cat incomplete.tsv | awk '$7 < 0.5' | cut -f4 >noreads.txt
