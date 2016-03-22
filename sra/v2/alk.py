#!/usr/bin/env python
"""
alk.py

Ranks projects by the range in the absolute difference between the
sample with the highest number of splice sites in ALK and the sample with the
lowest number of splice sites in ALK. Excludes all partially aligned samples
from consideration.
"""
import sys
import os
import gzip
from collections import defaultdict

if __name__ == '__main__':
    containing_dir = os.path.dirname(__file__)
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument(
            '--alk-junctions', type=str, required=False,
            default=os.path.join(containing_dir, 'alk_junctions.tsv.gz'),
            help=('junctions in alk gene; this should have been computed in '
                  'alk.sh')
        )
    parser.add_argument('--idmap', type=str, required=False,
            default=os.path.join(containing_dir,
                                    'intropolis.idmap.v2.hg38.tsv'),
            help='map from ids to SRA accession numbers'
        )
    parser.add_argument('--incomplete', type=str, required=False,
            default=os.path.join(containing_dir,
                                    'incomplete.tsv'),
            help=('output of incomplete.py; used to exclude samples from '
                  'consideration')
        )
    args = parser.parse_args()

sample_to_5p_splice_sites = defaultdict(set)
sample_to_3p_splice_sites = defaultdict(set)
project_to_sample = defaultdict(list)
index_to_sample = {}
incomplete_samples = set()
with open(args.incomplete) as incomplete_stream:
    incomplete_samples = set([line.strip().split('\t')[3]
                                for line in incomplete_stream.readlines()[1:]])
with open(args.idmap) as idmap_stream:
    for line in idmap_stream:
        tokens = line.strip().split('\t')
        if tokens[4] not in incomplete_samples:
            project_to_sample[tokens[1]].append(tokens[4])
            index_to_sample[tokens[0]] = tokens[4]
with gzip.open(args.alk_junctions) as alk_stream:
    for line in alk_stream:
        tokens = line.strip().split('\t')
        chrom, start, end, strand = tokens[0], tokens[1], tokens[2], tokens[3]
        for sample_index in tokens[-2].split(','):
            try:
                if strand == '-':
                    sample_to_5p_splice_sites[
                                        index_to_sample[sample_index]
                                    ].add((chrom, end))
                    sample_to_3p_splice_sites[
                                        index_to_sample[sample_index]
                                    ].add((chrom, start))
                elif strand == '+':
                    sample_to_5p_splice_sites[
                                        index_to_sample[sample_index]
                                    ].add((chrom, start))
                    sample_to_3p_splice_sites[
                                        index_to_sample[sample_index]
                                    ].add((chrom, end))
            except KeyError:
                continue
for project in project_to_sample:
    splice_site_counts = [(sample, len(sample_to_5p_splice_sites[sample])
                                    + len(sample_to_3p_splice_sites[sample]))
                                for sample in project_to_sample[project]]
    splice_site_counts.sort(key=lambda x: x[1], reverse=True)
    try:
        rank = splice_site_counts[0][1] - splice_site_counts[-1][1] 
    except IndexError:
        pass
    print '\t'.join([str(rank), project]
                    + [','.join(
                            (sample_and_count[0], str(sample_and_count[1]))
                        ) for sample_and_count in splice_site_counts])
