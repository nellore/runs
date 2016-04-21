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

def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return float(sum(sorted(lst)[quotient - 1:quotient + 1]) / 2)

def kmedian(to_cluster):
    to_cluster_size = len(to_cluster)
    if to_cluster_size == 1:
        return 0
    elif to_cluster_size == 2:
        return to_cluster[0][1] - to_cluster[1][1]
    else:
        to_cluster.sort(key=lambda x: x[1], reverse=True)
        counts = [el[1] for el in to_cluster]
        costs = [sum(
                [abs(el - counts[:i][i / 2])
                    for el in counts[:i]]
            ) + sum(
                [abs(el - counts[i:][(to_cluster_size - i) / 2])
                    for el in counts[i:]]
            ) for i in xrange(1, to_cluster_size)]
        partition_index = costs.index(min(costs)) + 1
        left_median = median(counts[:partition_index])
        right_median = median(counts[partition_index:])
        return left_median - right_median

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
totalexp = defaultdict(int)
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
        coverages = [int(token) for token in tokens[-1].split(',')]
        for i, sample_index in enumerate(tokens[-2].split(',')):
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
                totalexp[index_to_sample[sample_index]] += coverages[i]
            except KeyError:
                continue
for project in project_to_sample:
    splice_site_counts = [(sample, len(sample_to_5p_splice_sites[sample])
                                    + len(sample_to_3p_splice_sites[sample]))
                                for sample in project_to_sample[project]]
    splice_site_counts.sort(key=lambda x: x[1], reverse=True)
    rank = kmedian(splice_site_counts)
    if rank is not None:
        print '\t'.join([str(rank), project]
                        + [','.join(
                                (sample_and_count[0], str(sample_and_count[1]),
                                 str(totalexp[sample_and_count[0]]))
                            ) for sample_and_count in splice_site_counts])
