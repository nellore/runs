#!/usr/bin/env python
"""
incomplete.py

Finds GTEx samples that were incompletely downloaded and/or aligned by comparing
counts.tsv.gz files across batches with read counts in SraRunInfo.csv.

Tab-separated output fields:
1. project accession number
2. sample accession number
3. experiment accession number
4. run accession number
5. number of reads from SraRunInfo (i.e. mates for paired end samples)
6. number of reads Rail attempted to map
7. proportion of reads Rail attempted to map
8. 1 if LIKELY paired-end; else 0
9. 1 if SRA misreported paired-end status

Only runs for which Rail downloaded and aligned fewer than 100% of reads are
included.

We ran 

    pypy incomplete.py --sra-dir /path/to/gtex/dir |
        (read -r; printf "%s\n" "$REPLY"; sort -k7,7nr)
        >incomplete.tsv

Above, --gtex-dir is the path to the directory with the batch_* subdirs.
"""
import os
import csv
import gzip

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtex-dir', required=True,
        help='path to GTEx output files; this is where the batch_* '
             'subdirectories are')
    args = parser.parse_args()
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    srr_to_read_count = {}
    srr_to_line = {}
    srr_to_paired_status = {}
    with open(os.path.join(containing_dir, 'SraRunInfo.csv')) as run_stream:
        run_stream.readline()
        run_reader = csv.reader(run_stream, delimiter=',', quotechar='"')
        for tokens in run_reader:
            if not tokens or (len(tokens) == 1 and tokens[0] == ''): continue
            srr_to_line[tokens[0]] = '\t'.join(
                    [tokens[20], tokens[24], tokens[10], tokens[0]]
                )
            if tokens[15] == 'PAIRED':
                srr_to_read_count[tokens[0]] = int(tokens[3]) * 2
                srr_to_paired_status[tokens[0]] = '1'
            elif tokens[15] == 'SINGLE':
                srr_to_read_count[tokens[0]] = int(tokens[3])
                srr_to_paired_status[tokens[0]] = '0'
            else:
                raise RuntimeError(
                        'Fail: {} is neither SINGLE nor PAIRED'.format(
                                                                    tokens[15]
                                                                )
                    )
    print '\t'.join(['project', 'sample', 'experiment', 'run',
                        'read count as reported by SRA', 'reads aligned',
                        'proportion of reads reported by SRA aligned',
                        'paired-end', 'SRA misreported paired-end'])
    for i in xrange(100):
        with gzip.open(
                os.path.join(args.sra_dir,
                                'batch_{}'.format(i),
                                'cross_sample_results',
                                'counts.tsv.gz')
            ) as count_stream:
            count_stream.readline()
            for line in count_stream:
                srr = line.partition('\t')[0].partition('_')[0]
                read_count = int(
                        line.rpartition('\t')[2].partition(',')[0]
                    )
                mislabeled = False
                try:
                    ratio = float(read_count) / srr_to_read_count[srr]
                    if ratio > 1:
                        '''Mislabeled sample; whp recorded as SINGLE
                        when PAIRED'''
                        srr_to_read_count[srr] *= 2
                        ratio = float(read_count) / srr_to_read_count[srr]
                        srr_to_paired_status[srr] = '1'
                        mislabeled = True
                    elif ratio == 0.5:
                        '''Mislabeled sample; whp recorded as PAIRED
                        when SINGLE'''
                        srr_to_read_count[srr] /= 2
                        ratio = float(read_count) / srr_to_read_count[srr]
                        srr_to_paired_status[srr] = '0'
                        mislabeled = True
                except ZeroDivisionError:
                    ratio = 'NA'
                if read_count == srr_to_read_count[srr]:
                    # We got this one right
                    continue
                else:
                    print '\t'.join(map(str, [srr_to_line[srr],
                                        srr_to_read_count[srr],
                                        read_count, ratio,
                                        srr_to_paired_status[srr],
                                        '1' if mislabeled
                                        else '0']))
