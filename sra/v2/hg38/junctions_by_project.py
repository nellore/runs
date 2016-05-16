#!/usr/bin/env python
"""
junctions_by_project.py

Takes as input junction files from GTEx (first_pass_gtex_junctions.tsv.gz
output by gtex/combine_gtex.py) and public SRA (intropolis.v2.hg38.tsv.gz
output by sra/v2/hg38/combine_sra.py) in the following
tab-separated format:

1. chromosome
2. start position (1-based, inclusive)
3. end position (1-based, inclusive)
4. strand (+ or -)
5. start motif (e.g., GT)
6. end motif (e.g., AG)
7. comma-separated list of sample indexes
8. comma-separated list of coverages

as well as samples.tsv output by gtex/combine_gtex.py (mapping GTEx samples to
sample IDs in first_pass_gtex_junctions.tsv.gz)
and intropolis.idmap.v2.hg38.tsv output by sra/v2/hg38/combine_sra.py

and for each SRA project accession writes:
A) A BED.gz file assigning IDs to junctions. Columns are:
  1. chromosome
  2. start position (1-based, inclusive)
  3. end position (1-based, inclusive)
  4. ID (integer) -- this is the same across all projects
  5. the number 1000 (this is the score, and here just says use a dark color
                      in the UCSC genome browser)
  6. strand (+ or -)
B) A TSV.gz file mapping junction IDs to sample IDs and their corresponding
    coverages. Columns are:
  1. junction ID
  2. comma-separated list of IDs of samples in which junction is found
  3. comma-separated list of corresponding numbers of reads found to map across
      junction
"""
