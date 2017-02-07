#!/usr/bin/env bash
# Computes project stats in bash from intropolis.idmap.hg38.tsv
# Requires datamash: https://www.gnu.org/software/datamash/examples/ -- we used v1.1.1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
cut -f2 intropolis.idmap.v2.hg38.tsv | sort | uniq -c | rev | cut -d' ' -f2 | rev | sort -n | datamash -H mean 1 q1 1 median 1 q3 1 iqr 1 sstdev 1 jarque 1 >intropolis.project_stats.v2.hg38.txt
