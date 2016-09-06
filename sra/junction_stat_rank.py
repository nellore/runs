#!/usr/bin/env python
"""
junction_stat_rank.py

Computes junction stats by sample.
Reads intropolis.v1.hg19.tsv.gz from stdin; writes tab-separated output:
1. sample id
2. project accession number
3. sample accession number
4. experiment accession number
5. run accession number
6. total number of junctions
7. number of unique junctions contributed by the sample
8. average number of samples expressing junction

We ran gzip -cd intropolis.v1.hg19.tsv.gz
        | pypy unique_and_total_junctions.py | sort -k8,8n
        >junction_stats.tsv

junction_stats.tsv in this directory is the output.

We also executed

awk '$6 >= 100000 {print $0 "\t" $7/$6}' junction_stats.tsv 
    | sort -k9,9gr >junction_stats_with_unique_to_total_ratio.tsv

to add an extra column giving the ratio of unique junctions contributed to
total junctions contributed in those samples with over 100,000 junctions.
"""
import argparse
import os
import sys
from collections import defaultdict

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
            '--idmap', type=str, required=False,
            default=os.path.join(os.path.dirname(__file__),
                                 'intropolis.idmap.v1.hg19.tsv'),
            help='path to intropolis v1 ID map'
        )
    args = parser.parse_args()
    idmap = {}
    with open(args.idmap) as idmap_stream:
        for line in idmap_stream:
            current_id = line.partition('\t')[0]
            idmap[current_id] = line.strip()
    (total_junctions, unique_junctions,
        total_samples) = (defaultdict(int), defaultdict(int), defaultdict(int))
    for line in sys.stdin:
        tokens = line.strip().split('\t')
        samples = tokens[-2].split(',')
        num_samples = len(total_samples)
        for sample in samples:
            total_junctions[sample] += 1
            total_samples[sample] += num_samples
        if len(samples) == 1:
            unique_junctions[sample] += 1
    for sample in total_samples:
        print '\t'.join([idmap[sample]] + map(str, [total_junctions[sample],
                                              unique_junctions[sample],
                                              float(total_samples[sample])
                                                / total_junctions[sample]]))
