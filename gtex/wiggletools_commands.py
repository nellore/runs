#!/usr/bin/env python
"""
wiggletools_commands.py
Abhi Nellore / December 15, 2015

Outputs wiggletools commands for computing mean bigwigs by tissue.
"""
import gzip
from collections import defaultdict
import os

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtex-dir', required=True,
        help='path to GTEx output files; this is where the batch_* '
             'subdirectories are')
    parser.add_argument('--wiggletools', required=True,
        help='path to wiggletools')
    parser.add_argument('--out', required=True,
        help='path to output directory')
    args = parser.parse_args()
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    # Create original sample index to new sample index map
    sample_name_to_bw = {}
    big_name_to_sample_name = {}
    batch_numbers = []
    for manifest in glob.glob(os.path.join(containing_dir, '*.manifest')):
        batch_number = int(manifest.partition('.')[0].rpartition('_')[2])
        batch_numbers.append(batch_number)
        with open(manifest) as manifest_stream:
            for j, line in enumerate(manifest_stream):
                line = line.strip()
                if not line: continue
                sample_name = line.partition('\t')[0].partition(':')[2]
                big_name_to_sample_name[line.rpartition('\t')[2]] = sample_name 
                sample_name_to_bw[sample_name] = os.path.join(
                        args.gtex_dir,
                        'batch_{}'.format(batch_number),
                        'coverage_bigwigs',
                        line.rpartition('\t')[2],
                        '.bw'
                    )
    sample_name_to_mapped_reads = {}
    for batch_number in batch_numbers:
        with gzip.open(os.path.join(
                            args.gtex_dir,
                            'batch_{}'.format(batch_number),
                            'cross_sample_results',
                            'counts.tsv.gz'
                        )) as count_stream:
            count_stream.readline()
            for line in count_stream:
                tokens = line.strip().split('\t')
                sample_name_to_mapped_reads[
                        big_name_to_sample_name[tokens[0]]
                    ] = int(tokens[-2].split(',')[0])
    sample_name_to_tissue = {}
    with open(os.path.join(containing_dir, 'SraRunInfo.csv')) as sra_stream:
        sra_stream.readline()
        for line in sra_stream:
            line = line.strip()
            if not line: continue
            if '_rep1' in line or '_rep2' in line: continue
            tokens = line.split('\t')
            sample_name_to_tissue[tokens[0]] = tokens[41]
    # Handle exceptions: some samples in SraRunInfo don't have tissues labeled
    sample_name_to_tissue['SRR1325138'] = 'Skin'
    sample_name_to_tissue['SRR1325690'] = 'Stomach'
    sample_name_to_tissue['SRR1397115'] = 'Esophagus'
    sample_name_to_tissue['SRR1405266'] = 'Esophagus'
    sample_name_to_tissue['SRR1467633'] = 'Skin'
    tissue_to_sample_names = defaultdict(list)
    for sample_name in sample_name_to_mapped_reads:
        tissue_to_sample_names[
                sample_name_to_tissue[sample_name]
            ].append(sample_name)
    os.makedirs(args.out)
    for tissue in tissue_to_sample_names:
        print ' '.join([args.wiggletools, 'sum'] + [
                    'scale {} {}'.format(
                            sample_name_to_mapped_reads[sample_name],
                            sample_name_to_bw[sample_name]
                        ) for sample_name in sample_name_to_mapped_reads
                ]) + '>{}'.format(os.path.join(args.out, tissue + '.mean.wig'))