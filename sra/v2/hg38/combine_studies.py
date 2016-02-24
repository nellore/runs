#!/usr/bin/env python
"""
combine_studies.py

Combines GTEx and SRA junctions.
"""
import gzip
import os
from collections import defaultdict

def stream_to_list(stream):
    """ Converts each line in a stream to a list and yields it.

        stream: file handle

        Yield value: list of tab-separated values
    """
    for line in stream:
        yield line.strip().split('\t')

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--sra-junctions', type=str, nargs='+',
            required=True,
            help='space-separated gzipped SRA junction files IN BATCH ORDER'
        )
    parser.add_argument('--gtex-junctions', type=str, required=True,
            help='bid price for each m3.xlarge instance; this instance '
                 'type is used for preprocessing data'
        )
    parser.add_argument('--manifest-path', type=str, required=False,
            default=os.path.dirname(os.path.abspath(__file__)),
            help='path to manifest files for all of SRA job flow'
        )
    args = parser.parse_args()
    import tempfile
    temp_junction_path = tempfile.mkdtemp()
    import atexit
    import shutil
    atexit.register(shutil.rmtree, temp_junction_path)

    all_junctions = os.path.join(temp_junction_path, 'all_junctions.tsv.gz')
    batch_count = len(args.sra_junctions)
    index_to_index = [defaultdict(int) for _ in xrange(batch_count)]
    index_counter = 0
    for i in xrange(batch_count):
        old_index_counter = 0
        with open(os.path.join(
                        args.manifest_path,
                        'sra_batch_{}.manifest'.format(i)
                    )) as manifest_stream:
            for line in enumerate(manifest_stream):
                line = line.strip()
                if not line or line[0] == '#':
                    continue
                index_to_index[i][str(old_index_counter)] = str(index_counter)
                index_counter += 1
                old_index_counter += 1
    with gzip.open(all_junctions, 'w') as temp_stream:
        for i, sra_file in enumerate(sra_junctions):
            with gzip.open(sra_file) as sra_stream:
                for line in sra_file:
                    tokens = line.strip().split('\t')
                    tokens[-2] = [index_to_index[i][j]
                                    for j in tokens[-2].split(',')]
                    if int(tokens[2]) > int(tokens[1]):
                        print >>temp_stream, '\t'.join(
                                [tokens[0][:-1], tokens[1], tokens[2],
                                    tokens[0][-1], tokens[-2], tokens[-1], 's']
                            )
        with gzip.open(args.gtex_junctions) as gtex_stream:
            for line in gtex_stream:
                tokens = line.strip().split('\t')
                print >>temp_stream, '\t'.join(tokens[:4] + tokens[6:] + ['g'])
    import subprocess
    sorted_file = os.path.join(temp_junction_path,
                                'all_junctions.sorted.tsv.gz')
    subprocess.check_call('sort -k1,1 -k2,2n -k3,3n -T {} | gzip >{}'.format(
                                temp_junction_path,
                                sorted_file
                            ), shell=True, bufsize=-1)
    import itertools
    with gzip.open(sorted_file) as sorted_stream:
        for key, group in itertools.groupby(stream_to_list(sorted_stream),
                                            lambda x: x[:4]):
            gtex_samples = []
            gtex_coverages = []
            sra_samples = []
            sra_coverages = []
            for tokens in group:
                if tokens[-1] == 'g':
                    gtex_samples.append(tokens[-3])
                    gtex_coverages.append(tokens[-2])
                else:
                    sra_samples.append(tokens[-3])
                    sra_coverages.append(tokens[-2])
            print '\t'.join(key + [','.join(gtex_samples),
                                   ','.join(gtex_coverages),
                                   ','.join(sra_samples),
                                   ','.join(sra_coverages)])
