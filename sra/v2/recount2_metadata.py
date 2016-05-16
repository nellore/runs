#!/usr/bin/env python
"""
recount2_metadata.py

Dumps SRA metadata for Recount 2. Based on incomplete.py. Uses metadata from
SHARQ, developed by Carl Kingsford's group. Requires that AUC.sh has been
executed.

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
10. mapped read count
11. AUC
12. SHARQ tissue if available else NA
13: SHARQ cell type if available else NA
14. submission date (from biosample DB)
15. publication date (from biosample DB)
16. update date (from biosample DB)

We ran 

    pypy recount2_metadata.py --sra-dir /path/to/sra/dir |
        (read -r; printf "%s\n" "$REPLY"; sort -k7,7nr)
        >recount2_metadata.tsv

Above, --sra-dir is the path to the directory with the batch_* subdirs.
"""
import os
import csv
import gzip
from collections import defaultdict

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sra-dir', required=True,
        help='path to SRA output files; this is where the batch_* '
             'subdirectories are')
    args = parser.parse_args()
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    # Get mappings from srs to srr
    srs_to_srr = defaultdict(list)
    with open(os.path.join(containing_dir,
                            'intropolis.idmap.v2.hg38.tsv')) as idmap_stream:
        for line in idmap_stream:
            tokens = line.strip().split('\t')
            srs_to_srr[tokens[2]].append(tokens[4])
    srr_to_submission_date = defaultdict(lambda: 'NA')
    srr_to_publication_date = defaultdict(lambda: 'NA')
    srr_to_update_date = defaultdict(lambda: 'NA')
    # Get dates from biosample
    with open(os.path.join(containing_dir,
                            'hg38', 'biosample_tags.tsv')) as biosample_stream:
        biosample_stream.readline()
        for line in biosample_stream:
            tokens = line.strip().split('\t')
            for srr in srs_to_srr[tokens[9]]:
                srr_to_submission_date[srr] = tokens[10]
                srr_to_publication_date[srr] = tokens[11]
                srr_to_update_date[srr] = tokens[12]
    srr_to_tissue, srr_to_cell_type =  (defaultdict(lambda: 'NA'),
                                            defaultdict(lambda: 'NA'))
    with open(os.path.join(containing_dir,
                            'sra-all-fields-2015-9-13.txt')) as sharq_stream:
        sharq_stream.readline()
        sharq_reader = csv.reader(sharq_stream, delimiter=',', quotechar='"')
        for tokens in sharq_reader:
            srr_to_tissue[tokens[0]] = tokens[5]
            srr_to_cell_type[tokens[0]] = tokens[6]
    srr_to_auc = defaultdict(lambda: 'NA')
    with open(os.path.join(containing_dir, 'auc.tsv')) as auc_stream:
        for line in auc_stream:
            tokens = line.strip().split('\t')
            srr_to_auc[tokens[0].split('.')[0]] = str(
                            int(float(tokens[1].strip()))
                        )
    srr_to_read_count = {}
    srr_to_line = {}
    srr_to_paired_status = {}
    srr_to_GSM = defaultdict(lambda: 'NA')
    with open(os.path.join(containing_dir, 'hg38',
                            'SraRunInfo.csv')) as run_stream:
        run_stream.readline()
        run_reader = csv.reader(run_stream, delimiter=',', quotechar='"')
        for tokens in run_reader:
            if not tokens or (len(tokens) == 1 and tokens[0] == ''): continue
            srr_to_line[tokens[0]] = '\t'.join(
                    [tokens[20], tokens[24], tokens[10], tokens[0]]
                )
            if tokens[29]:
                srr_to_line[tokens[0]] = tokens[29]
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
                        'paired-end', 'SRA misreported paired-end',
                        'mapped read count', 'AUC',
                        'SHARQ tissue', 'SHARQ cell type',
                        'Biosample submission date',
                        'Biosample publication date',
                        'Biosample update date', 'GSM'])
    srr_to_mapped_reads = defaultdict(int)
    for i in xrange(100):
        with gzip.open(
                os.path.join(args.sra_dir,
                                'batch_{}'.format(i),
                                'cross_sample_results',
                                'counts.tsv.gz')
            ) as count_stream:
            count_stream.readline()
            for line in count_stream:
                srr = line.partition('\t')[0]
                tokens = line.strip().split('\t')
                read_count = int(tokens[-1].partition(',')[0])
                srr_to_mapped_reads[srr] = int(tokens[-2].partition(',')[0])
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
                print '\t'.join(map(str, [srr_to_line[srr],
                                    srr_to_read_count[srr],
                                    read_count, ratio,
                                    srr_to_paired_status[srr],
                                    '1' if mislabeled
                                    else '0', srr_to_mapped_reads[srr],
                                    srr_to_auc[srr],
                                    srr_to_tissue[srr],
                                    srr_to_cell_type[srr],
                                    srr_to_submission_date[srr],
                                    srr_to_publication_date[srr],
                                    srr_to_update_date[srr]]))
