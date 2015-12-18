#!/usr/bin/env python
"""
segregate_gtex.py
Abhi Nellore / December 16, 2015

Segregates junctions from any output file of combine_gtex.py by site from
SraRunInfo.csv. Requires that combine_gtex.py has already been executed.

The junction input file is read from stdin and takes the following form:
1. chromosome
2. start position if +; end position if -
3. end position if +; start position if -
4. strand (+ or -)
5. start motif (e.g., GT)
6. end motif (e.g., AG)
7. comma-separated list of sample indexes
8. comma-separated list of coverages

The input is SORTED by fields 1 and 2. It is obtained from the output of
combine_gtex.py by running the following commands.
gzip -cd first_pass_gtex_junctions.tsv.gz
| awk '$4 == "-" {holder=$3; $3=$2; $2=$holder; printf $1;
    for(i=2;i<=NF;i++){printf "\t" $i}} $4 == "+" {print $0}'
| sort -k1,1 -k2,2n 

Output files are written to some output directory --out in the same format,
but preceded by two fields:
1. Ensembl gene IDs for start/end positions of junction. Comma-separated; NA
if not available.
2. [Overlap group index].[junction index], where the overlap group index
corresponds to a unique group of introns with the same start position, and the
junction index is the index of a junction within a group
"""
import gzip
from collections import defaultdict
import os
import re
import sys

def split_iterator(stream):
    for line in stream:
        yield line.strip().split('\t')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output-dir', required=True,
        help='directory in which to write output files')
    parser.add_argument('--gtf', required=True,
        help='path to gtf.gz file; we used '
             'ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/'
             'Homo_sapiens.GRCh38.83.gtf.gz')
    parser.add_argument('--map', required=True,
        help='path to file mapping sample indexes to SRRs; this is '
             'samples.tsv output by combine_gtex.py')
    args = parser.parse_args()
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    # Map Ensembl chrs to iGenome chrs
    ensembl_to_igenome = {}
    with open(os.path.join(
                    containing_dir, 'DER_analysis/coverageMatrix/genomicState/'
                                    'hg38.ucsc.sizes.ensembl.gencode')
            ) as chrmap_stream:
        chrmap_stream.readline()
        for line in chrmap_stream:
            tokens = line.strip().split('\t')
            ensembl_to_igenome[tokens[2]] = tokens[0]
    bounds = defaultdict(list)
    bounds_to_gene = {}
    with gzip.open(args.gtf) as gtf_stream:
        for line in gtf_stream:
            if line[0] == '#': continue
            tokens = line.strip().split('\t')
            if tokens[0] not in ensembl_to_igenome: continue
            if tokens[2] != 'gene': continue
            gene = tokens[8].split('"')[1]
            start = int(tokens[3])
            end = int(tokens[4]) + 1 # so it's not inclusive
            assert gene not in bounds
            bound = (ensembl_to_igenome[tokens[0]], start, end)
            bounds[gene] = bound
            bounds_to_gene[bound] = gene
    sorted_bounds = sorted(bounds.items(), key=lambda x: x[1])
    # Ensure gene boundaries do not overlap
    for i in xrange(1, len(sorted_bounds)):
        if sorted_bounds[i][1][0] == sorted_bounds[i-1][1][0]:
            assert sorted_bounds[i][1][0] > sorted_bounds[i-1][1][1]
    chrom_bounds = defaultdict(list)
    for gene in bounds:
        chrom_bounds[bounds[gene][0]].extend([bounds[gene][1],
                                                bounds[gene][2]])
    for chrom in chrom_bounds:
        chrom_bounds[chrom].sort()
    # Map sample index to tissue
    index_to_run = {}
    run_to_index = {}
    with open(args.map) as map_stream:
        for line in map_stream:
            tokens = line.strip().split('\t')
            index_to_run[tokens[0]] = tokens[1]
            run_to_index[tokens[1]] = tokens[0]
    index_to_site = {}
    with open(os.path.join(containing_dir, 'SraRunInfo.csv')) as run_stream:
        run_stream.readline()
        for line in run_stream:
            line = line.strip()
            if not line: continue
            tokens = line.split(',')
            try:
                index_to_site[run_to_index[tokens[0]]] = tokens[42]
            except KeyError:
                pass
    output_handles = {}
    overlap_group = 0
    import itertools
    from bisect import bisect_left
    for key, group \
        in itertools.groupby(split_iterator(sys.stdin),
                                    lambda x: (x[0], x[1])):
        total_coverages = defaultdict(int)
        write_data = []
        for junction_index, tokens in enumerate(group):
            samples_and_coverages = zip(tokens[-2].split(','),
                                        tokens[-1].split(','))
            for sample, coverage in samples_and_coverages:
                total_coverages[sample] += int(coverage)
            junction = tokens[:6]
            start = int(junction[1])
            end = int(junction[2])
            index = bisect_left(chrom_bounds[junction[0]], start)
            try:
                start_gene = bounds_to_gene[(junction[0],
                    chrom_bounds[junction[0]][index-1],
                    chrom_bounds[junction[0]][index])]
            except (IndexError, KeyError):
                start_gene = 'NA'
            index = bisect_left(chrom_bounds[junction[0]], end)
            try:
                end_gene = bounds_to_gene[(junction[0],
                    chrom_bounds[junction[0]][index-1],
                    chrom_bounds[junction[0]][index])]
            except (IndexError, KeyError):
                end_gene = 'NA'
            sites = defaultdict(list)
            for sample, coverage in samples_and_coverages:
                sites[index_to_site[sample]].append((sample, coverage))
            write_data.append([junction, sites, start_gene, end_gene])
        for junction, sites, start_gene, end_gene in write_data:
            for site in sites:
                line_to_write = '\t'.join(
                                    [start_gene + ',' + end_gene,
                                        str(overlap_group)]
                                        + junction + [','.join(
                                                [el[0] for el
                                                    in sites[site]]
                                            ), ','.join(
                                                [el[1] for el
                                                    in sites[site]]
                                            ), ','.join(
                                                [str(
                                                    float(el[1])
                                                    / total_coverages[el[0]]
                                                    ) for el in sites[site]]
                                                )]
                                )
                try:
                    print >>output_handles[site], line_to_write
                except KeyError:
                    output_handles[site] = gzip.open(
                            os.path.join(args.output_dir, 
                                re.sub('[^a-zA-Z\d:]+', '_',
                                        site.lower().strip()
                                        ).strip('_')
                        ) + '_gtex_junctions.tsv.gz', 'w'
                    )
        overlap_group += 1