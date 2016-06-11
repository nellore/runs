#!/usr/bin/env python
"""
junctions_by_project.py

Takes as input junction files from GTEx (first_pass_gtex_junctions.tsv.gz
output by gtex/combine_gtex.py) and public SRA (intropolis.v2.hg38.tsv.gz
output by sra/v2/hg38/combine_sra.py) in the following
tab-separated format:

1. chromosome
2. start position (1-based, inclusive)
3. end position (1-based, inclusive)
4. strand (+ or -)
5. start motif (e.g., GT)
6. end motif (e.g., AG)
7. comma-separated list of sample indexes
8. comma-separated list of coverages

as well as samples.tsv output by gtex/combine_gtex.py (mapping GTEx samples to
sample IDs in first_pass_gtex_junctions.tsv.gz)
and intropolis.idmap.v2.hg38.tsv output by sra/v2/hg38/combine_sra.py

and for each SRA project accession writes:
A) A BED file assigning IDs to junctions. Columns are:
  1. chromosome
  2. start position (1-based, inclusive)
  3. end position (1-based, inclusive)
  4. ID (integer) -- this is the same across all projects
  5. the number 1000 (this is the score, and here just says use a dark color
                      in the UCSC genome browser)
  6. strand (+ or -)
B) A TSV file mapping junction IDs to sample IDs and their corresponding
    coverages. Columns are:
  1. junction ID
  2. comma-separated list of IDs of samples in which junction is found
  3. comma-separated list of corresponding numbers of reads found to map across
      junction
"""
import os
import gzip
import tempfile
import subprocess
import itertools
from collections import defaultdict
import shutil
import atexit

_gtex_project_id = 'SRP012682'

def stream_to_tokens(input_stream):
    """ Converts newline-separated TSV input stream into token generator

        input_stream: open file handle

        Yield value: tab-separated tokens from line
    """
    for line in input_stream:
        yield line.strip().split('\t')

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--gtex-junctions', type=str, required=True,
            help=('path to first-pass GTEx junctions file output by '
                  'combine_gtex.py')
        )
    parser.add_argument('--sra-junctions', type=str, required=True,
            help=('path to first-pass SRA junctions file output by '
                  'combine_sra.py')
        )
    parser.add_argument('--gtex-ids', type=str, required=True,
            help='path to GTEx IDs file output by combine_gtex.py'
        )
    parser.add_argument('--sra-ids', type=str, required=True,
            help='path to SRA IDs file output by combine_sra.py'
        )
    parser.add_argument('--output-dir', type=str, required=True,
            help=('directory in which to dump all output files; created if '
                  'it doesn\'t already exist')
        )
    parser.add_argument('--temp-dir', type=str, required=False,
            default=None,
            help=('where to store temporary files; None if Python should '
                  'decide')
        )
    parser.add_argument('--projects-per-batch', type=int, required=False,
            default=500,
            help='1/2 maximum number of open file handles'
        )
    args = parser.parse_args()

    # Create output dir if it doesn't exist
    try:
        os.makedirs(args.output_dir)
    except OSError as e:
        if os.path.isdir(args.output_dir): pass

    temp_dir = tempfile.mkdtemp(dir=args.temp_dir)
    atexit.register(shutil.rmtree, temp_dir)
    merged_junctions = os.path.join(temp_dir, 'merged_junctions.tsv')
    merged_sorted_junctions = os.path.join(temp_dir,
                                            'merged_sorted_junctions.tsv')
    subprocess.check_output(
        'gzip -cd {} | awk \'{{print $0 "\\tg"}}\' >>{}'.format(
                                                    args.gtex_junctions,
                                                    merged_junctions
                                                ),
                                executable='/bin/bash',
                                shell=True)
    subprocess.check_output(
        'gzip -cd {} | awk \'{{print $0 "\\ts"}}\' >>{}'.format(
                                                    args.sra_junctions,
                                                    merged_junctions
                                                ),
                                executable='/bin/bash',
                                shell=True)
    subprocess.check_output('sort {}-k1,1 -k2,2n -k3,3n {} >{}'.format(
                                                    '-T %s ' % args.temp_dir
                                                    if args.temp_dir
                                                    is not None else '', 
                                                    merged_junctions,
                                                    merged_sorted_junctions
                                                ),
                                executable='/bin/bash',
                                shell=True
        )

    '''Map SRA sample numbers to projects; note GTEx project number is always
    SRP012682, so we don't need to do the same thing for GTEx'''
    with open(args.sra_ids) as id_stream:
        id_to_srp = { tokens[0] : tokens[1] for tokens in
                        (line.strip().split('\t') for line in id_stream) }

    srps = id_to_srp.values()
    srp_batches = [srps[i:i+args.projects_per_batch]
                    for i in xrange(0, len(srps), args.projects_per_batch)]
    srp_batches[-1].append(_gtex_project_id)
    score = '1000'

    # Number by which to increment every GTEx sample id
    translation = max([int(sample_id) for sample_id in id_to_srp.keys()]) + 1

    # Print all sample ids
    with open(
                os.path.join(args.output_dir, 'sample_ids.tsv'), 'w'
            ) as sample_id_stream:
        with open(args.sra_ids) as sra_id_stream:
            for line in sra_id_stream:
                tokens = line.strip().split('\t')
                print >>sample_id_stream, '\t'.join(
                        [tokens[0], tokens[1], tokens[-1]]
                    )
        with open(args.gtex_ids) as gtex_id_stream:
            for line in gtex_id_stream:
                tokens = line.strip().split('\t')
                print >>sample_id_stream, '\t'.join(
                        [str(int(tokens[0]) + translation), _gtex_project_id,
                            tokens[1]]
                    )
    for srp_batch in srp_batches:
        # For storing file handles
        project_bed_handles = {}
        project_coverage_handles = {}
        junction_id = 0
        try:
            srp_batch_set = set(srp_batch)
            with open(merged_sorted_junctions) as merged_stream:
                for key, group in itertools.groupby(
                                    stream_to_tokens(merged_stream),
                                    lambda x: x[:6]
                                ):
                    junction_id_string = str(junction_id)
                    junction = '\t'.join(
                            [key[0], key[1], key[2],
                                junction_id_string, score, key[3]]
                        )
                    for tokens in group:
                        if tokens[-1] == 'g':
                            if _gtex_project_id not in srp_batch_set: continue
                            try:
                                print >>project_bed_handles[
                                                _gtex_project_id
                                            ], junction
                            except KeyError:
                                project_bed_handles[
                                        _gtex_project_id
                                    ] = open(
                                        os.path.join(
                                                args.output_dir,
                                                _gtex_project_id
                                                + '.junction_id.bed'
                                            ), 'w'
                                    )
                                print >>project_bed_handles[
                                                    _gtex_project_id
                                                ], junction
                            current_samples = [
                                        str(int(token) + translation)
                                        for token in tokens[-3].split(',')
                                    ]
                            try:
                                print >>project_coverage_handles[
                                        _gtex_project_id
                                    ], '\t'.join(
                                            [junction_id_string,
                                                ','.join(current_samples),
                                                tokens[-2]]
                                        )
                            except KeyError:
                                project_coverage_handles[
                                            _gtex_project_id
                                        ] = open(os.path.join(
                                                    args.output_dir,
                                                    _gtex_project_id
                                                + '.junction_coverage.tsv'
                                            ), 'w'
                                        )
                                print >>project_coverage_handles[
                                        _gtex_project_id
                                    ], '\t'.join(
                                            [junction_id_string,
                                                ','.join(current_samples),
                                                tokens[-2]]
                                        )
                        else:
                            # SRA
                            project_to_samples = defaultdict(list)
                            project_to_coverages = defaultdict(list)
                            for sample, coverage in zip(tokens[-3].split(','),
                                                        tokens[-2].split(',')):
                                project_to_samples[id_to_srp[sample]].append(
                                        sample
                                    )
                                project_to_coverages[id_to_srp[sample]].append(
                                        coverage
                                    )
                            for project in project_to_samples:
                                if project not in srp_batch_set: continue
                                try:
                                    print >>project_bed_handles[project], (
                                            junction
                                        )
                                except KeyError:
                                    project_bed_handles[project] = open(
                                            os.path.join(
                                                args.output_dir,
                                                project + '.junction_id.bed'
                                            ), 'w'
                                        )
                                    print >>project_bed_handles[project], (
                                            junction
                                        )
                                try:
                                    print >>project_coverage_handles[
                                            project
                                        ], '\t'.join(
                                             [junction_id_string,
                                                ','.join(
                                                    project_to_samples[project]
                                                ),
                                                ','.join(
                                                project_to_coverages[project]
                                            )]
                                        )
                                except KeyError:
                                    project_coverage_handles[
                                                project
                                            ] = open(
                                                os.path.join(
                                                    args.output_dir,
                                                    project
                                                + '.junction_coverage.tsv'
                                            ), 'w'
                                        )
                                    print >>project_coverage_handles[
                                            project
                                        ], '\t'.join(
                                             [junction_id_string,
                                                ','.join(
                                                    project_to_samples[project]
                                                ),
                                                ','.join(
                                                project_to_coverages[project]
                                            )]
                                        )
                    junction_id += 1
        finally:
            for project in project_bed_handles:
                project_bed_handles[project].close()
            for project in project_coverage_handles:
                project_coverage_handles[project].close()
