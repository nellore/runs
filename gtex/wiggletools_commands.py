#!/usr/bin/env python
"""
wiggletools_commands.py
Abhi Nellore / December 15, 2015

Writes wiggletools commands for computing mean bigwigs by tissue. Each set
of commands is numbered. They should be executed in order; some commands
in successive files depend on commands from previous files.
"""
import gzip
from collections import defaultdict
import os
import glob
import sys

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtex-dir', required=True,
        help='path to GTEx output files; this is where the batch_* '
             'subdirectories are')
    parser.add_argument('--auc', required=True,
        help='path to file with bigwig AUCs for normalization; '
             'this script normalizes to library size 40 million 100-bp reads')
    parser.add_argument('--wiggletools', required=True,
        help='path to wiggletools')
    parser.add_argument('--max-bw', required=False,
        default=500,
        help='max number of bigwig files to process at a time')
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
                        line.rpartition('\t')[2] + 
                        '.bw'
                    )
    sample_name_to_auc = {}
    with open(args.auc) as auc_stream:
        for line in auc_stream:
            tokens = line.strip().split('\t')
            sample_name_to_auc[tokens[0].strip()] = float(tokens[1])
    sample_name_to_tissue = {}
    with open(os.path.join(containing_dir, 'SraRunInfo.csv')) as sra_stream:
        sra_stream.readline()
        for line in sra_stream:
            line = line.strip()
            if not line: continue
            if '_rep1' in line or '_rep2' in line: continue
            tokens = line.split(',')
            sample_name_to_tissue[tokens[0]] = tokens[41]
    # Handle exceptions: some samples in SraRunInfo don't have tissues labeled
    sample_name_to_tissue['SRR1325138'] = 'Skin'
    sample_name_to_tissue['SRR1325690'] = 'Stomach'
    sample_name_to_tissue['SRR1397115'] = 'Esophagus'
    sample_name_to_tissue['SRR1405266'] = 'Esophagus'
    sample_name_to_tissue['SRR1467633'] = 'Skin'
    tissue_to_sample_names = defaultdict(list)
    for sample_name in sample_name_to_auc:
        tissue_to_sample_names[
                sample_name_to_tissue[sample_name]
            ].append(sample_name)
    try:
        os.makedirs(args.out)
    except OSError as e:
        if 'File exists' not in e:
            raise
    file_handles = [
            open(os.path.join(args.out, 'wiggletools_commands_0'), 'w')
        ]
    for tissue in tissue_to_sample_names:
        file_index = 0
        divided_sample_names = [
                tissue_to_sample_names[tissue][i:i+args.max_bw]
                for i in xrange(0, len(tissue_to_sample_names[tissue]),
                                args.max_bw)
            ]
        # Remove a lonely sample
        if (len(divided_sample_names) >= 2
            and len(divided_sample_names[-1]) == 1):
            divided_sample_names[-2].append(divided_sample_names[-1][0])
            divided_sample_names = divided_sample_names[:-1]
        for i, sample_group in enumerate(divided_sample_names):
            command_to_print = ' '.join([args.wiggletools, 'sum'] + [
                        'scale {} {}'.format(
                                float(40000000) * 100
                                    / sample_name_to_auc[sample_name],
                                sample_name_to_bw[sample_name]
                            ) for sample_name in tissue_to_sample_names[tissue]
                    ]) + ' >{}'.format(
                                    os.path.join(args.out,
                                        tissue.replace(' ', '_') + '.sum.wig_'
                                        + str(i))
                                )
            try:
                print >>file_handles[i], command_to_print
            except IndexError:
                file_handles.append(
                    open(os.path.join(args.out,
                                        'wiggletools_commands_' + str(i)), 'w')
                )
                print >>file_handles[i], command_to_print
    next_index = len(file_handles)
    file_handles.append(
                open(os.path.join(args.out,
                                    'wiggletools_commands_'
                                    + str(next_index)), 'w')
            )
    import glob
    for tissue in tissue_to_sample_names:
        print >>file_handles[-1], ' '.join([args.wiggletools, 'scale',
                                str(1./len(tissue_to_sample_names[tissue])),
                                ('sum ' + os.path.join(args.out,
                                 tissue.replace(' ', '_') + '.sum.wig_*'))
                                if len(glob.glob(os.path.join(args.out,
                                 tissue.replace(' ', '_') + '.sum.wig_*')))
                                >= 2 else os.path.join(args.out,
                                 tissue.replace(' ', '_') + '.sum.wig_*'),
                                '>' + os.path.join(args.out,
                                    tissue.replace(' ', '_') + '.mean.wig')])
    file_handles.append(
                open(os.path.join(args.out,
                                    'wiggletools_commands_'
                                    + str(next_index+1)), 'w')
            )
    sample_count = len(sample_name_to_auc)
    print >>file_handles[-1], ' '.join([args.wiggletools, 'sum']
                                    + ['scale {} {}'.format(
                                    float(
                                        len(tissue_to_sample_names[tissue])
                                ) / sample_count,
                                    os.path.join(args.out,
                                    tissue.replace(' ', '_') + '.mean.wig')
                                ) for tissue in tissue_to_sample_names]) + (
                                        ' >{}'.format(
                                                os.path.join(args.out,
                                                'mean.wig')
                                            )
                                    )
    for file_handle in file_handles:
        file_handle.close()

