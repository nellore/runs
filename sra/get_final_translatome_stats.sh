#!/usr/bin/env bash
# We ran sh get_translatome_stats_for_table.sh >translatome_stats_for_table.tsv
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 also in intropolis:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 also in intropolis and in annotation:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$15 == 1' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 also in intropolis and annotated exonSkip:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$13 == 1 && $14 == 1 && $15 == 0' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 also in intropolis and annotated altStartEnd:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$13 == 0 && $14 == 1 || $14 == 0 && $13 == 1' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 also in intropolis and novel:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$13 == 0 && $14 == 0 && $15 == 0' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 in >= 1000 samples from intropolis:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$7 >= 1000' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 in >= 1000 samples from intropolis and in annotation:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$7 >= 1000' | awk '$15 == 1' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 in >= 1000 samples from intropolis and annotated exonSkip:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$7 >= 1000' | awk '$13 == 1 && $14 == 1 && $15 == 0' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 in >= 1000 samples from intropolis and annotated altStartEnd:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$7 >= 1000' | awk '$13 == 0 && $14 == 1 || $14 == 0 && $13 == 1' | wc -l
echo "Number of translatome junctions successfully lifted over from mm10 to hg19 in >= 1000 samples from intropolis and novel:"
gzip -cd translatome_mm10_to_hg19_junctions.tsv.gz | awk '$7 >= 1000' | awk '$13 == 0 && $14 == 0 && $15 == 0' | wc -l
