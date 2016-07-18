#!/usr/bin/env python
"""
mismatch_means.py

Writes commands for finding means of mismatch bigWig tracks by weighting mean
mismatch bigWigs for different batches. Excludes samples for which _no_ reads
were downloaded, but includes samples for which at least half of reads were
downloaded. Relies on  manifest files, incomplete.tsv (which was generated by
incomplete.sh), coverage bigWigs from Rail-RNA run on 50K samples, and
wiggletools.
"""
import argparse
import glob

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
            '--incomplete', type=str, required=True,
            help=('path to incomplete.tsv, list of incompletely downloaded '
                  'samples')
        )
    parser.add_argument(
            '--manifest', type=str, required=True,
            help=('path to directory containing manifest files from Rail-RNA '
                  'run on 50K samples')
        )
    parser.add_argument(
            '--batch', type=str, required=True,
            help=('path to directory containing batch_* subdirectories from '
                  'Rail-RNA run on 50K SRA samples')
        )
    parser.add_argument(
            '--output', type=str, required=True,
            help='output directory'
        )
    parser.add_argument(
            '--wiggletools', type=str, required=True,
            help='path to wiggletools executable; use v1.1'
        )
    args = parser.parse_args()
    with open(args.incomplete) as incomplete_stream:
        incomplete_srrs = set(
                [tokens[3] for tokens in 
                    [line.strip().split('\t') for line in incomplete_stream]
                    if float(tokens[6]) < 0.5]
            )
    manifest_files = [os.path.join(
            args.manifest, 'sra_batch_{}.manifest'.format(i)
        ) for i in xrange(100)]
    manifest_sizes = []
    for manifest_file in manifest_files:
        with open(manifest_file) as manifest_stream:
            srrs = [tokens[0].partition(':')[2] for tokens in
                    [line.strip().split('\t') for line in manifest_stream]]
            manifest_sizes.append(
                    len(filter(lambda x: x not in incomplete_srrs, srrs))
                )
    total_samples = sum(manifest_sizes)
    weights = ['%.12f' % (float(manifest_size) / total_samples)
                for manifest_size in manifest_sizes]
    for base in 'ACGTN':
        'sum', 'scale {weight} {batch_dir}/batch_{batch_number}/coverage_bigwigs/'.format(weight)
