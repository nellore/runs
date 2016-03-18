#!/usr/bin/env python
"""
incomplete.py

Finds SRA samples that were incompletely downloaded and/or aligned by comparing
counts.tsv.gz files across batches with read counts in SraRunInfo.csv.

Tab-separated output fields:
1. run accession number
2. number of reads from SraRunInfo (i.e. mates for paired end samples)
3. number of reads Rail attempted to map
4. proportion of reads Rail attempted to map

Only runs for which Rail downloaded and aligned fewer than 100% of reads are
included.

We ran pypy incomplete.py | sort -k4,4er >incomplete.tsv
"""
import os
import csv
import gzip

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sra-dir', required=True,
        help='path to SRA output files; this is where the batch_* '
             'subdirectories are')
    args = parser.parse_args()
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    srr_to_read_count = {}
    with open(os.path.join(containing_dir, 'hg38',
                            'SraRunInfo.csv')) as run_stream:
        run_stream.readline()
        run_reader = csv.reader(run_stream, delimiter=',', quotechar='"')
        for tokens in run_reader:
            if not tokens or (len(tokens) == 1 and tokens[0] == ''): continue
            if tokens[15] == 'PAIRED':
                srr_to_read_count[tokens[0]] = int(tokens[3]) * 2
            elif tokens[15] == 'SINGLE':
                srr_to_read_count[tokens[0]] = int(tokens[3])
    for i in xrange(100):
        with gzip.open(
                os.path.join(args.sra_dir,
                                'batch_{}'.format(i),
                                'cross_sample_results'
                                'counts.tsv.gz')
            ) as count_stream:
            count_stream.readline()
            for line in count_stream:
                srr = line.partition('\t')[0]
                read_count = int(
                        line.rpartition('\t')[2][-1].partition(',')[0]
                    )
                if read_count == srr_to_read_count[srr]:
                    # We got this one right
                    continue
                else:
                    print '\t'.join([srr, srr_to_read_count[srr],
                                        read_count,
                                        float(read_count)
                                            / srr_to_read_count[srr]])
