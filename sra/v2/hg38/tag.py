#!/usr/bin/env python
"""
tag.py

Defines clusters of tags from NCBI Biosample metadata and writes a binary
variable for each SRA sample accession number: 1 means at least one tag from
the cluster is present in the sample's metadata, and 0 means no tags from
the cluster are present. Operates on output of xmlparse.py . Tags are provided
in header line.
"""
import sys
import re

if __name__ == '__main__':
    # Define tag clusters; each variable is a different cluster of keywords
    cell_line = set(['line:', 'cellline', 'cell line', 'a549', '\\bt24_',
                        'geuv', 'hela', 'hesc', 'hepg2', 'hep g2'])
    small_rna = set(['mirna', 'microrna', 'small rna', '\\bsrna'])
    single_cell = set(['single-cell', 'single cell'])
    fetal = set(['fetal', 'fetus'])
    stem_cell = set(['hesc', 'stem cell' 'stem-cell', 'ipsc', 'pluripotent'])
    primary = set(['patient', 'subject', 'primary tissue', 'donor'])
    cancer = set(['tumor', 'cancer', 'oma\\b'])
    polyA = set(['polya', 'poly-a'])
    totalrna = set(['total rna', 'ribozero', 'ribo-zero', 'ribominus', 'ribo-minus'])
    regexps = ['|'.join(cluster) for cluster
    			in [cell_line, small_rna, single_cell, fetal, stem_cell,
                    primary, cancer, polyA, totalrna]]
    original_header = sys.stdin.readline().strip().split('\t')
    header = ['cell line', 'small rna', 'single cell', 'fetal',
                'stem cell', 'primary', 'cancer', 'total RNA', 'polyA'
                ] + original_header
    print '\t'.join(header)
    for line in sys.stdin:
    	lower_line = line.lower()
        for regexp in regexps:
        	if re.search(regexp, lower_line):
        		sys.stdout.write('1\t')
        	else:
        		sys.stdout.write('0\t')
        sys.stdout.write(line)