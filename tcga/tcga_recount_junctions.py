#!/usr/bin/env python
"""
tcga_recount_junctions.py

Adds TCGA junctions to recount; requires specifying where an old version of 
recount is so old junction IDs can be read.
"""
import atexit
import tempfile
import os
import subprocess
import glob
import gzip
import itertools

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--recount-dir', type=str, required=True,
            help=('path to the last version of recount from which junctions '
                  'should be ripped')
        )
    parser.add_argument('--tcga-junctions', type=str, required=True,
            help='path to first-pass TCGA junctions file output by '
                 'combine_tcga.py; should be first_pass_tcga_junctions.tsv.gz'
        )
    parser.add_argument('--tcga-ids', type=str, required=True,
            help='path to TCGA sample IDs file output by '
                 'combine_tcga.py; should be in same directory as '
                 'first_pass_tcga_junctions.tsv.gz')
    parser.add_argument('--output-dir', type=str, required=True,
            help=('directory in which to dump all output files; created if '
                  'it doesn\'t already exist')
        )
    parser.add_argument('--temp-dir', type=str, required=False,
            default=None,
            help=('where to store temporary files; None if Python should '
                  'decide')
        )
    parser.add_argument('--junction-limit', type=int, required=False,
            default=3000000,
            help=('number of junctions to read from BED files before dumping '
                  'to file')
        )
    args = parser.parse_args()

    # Create output dir if it doesn't exist
    try:
        os.makedirs(args.output_dir)
    except OSError as e:
        if os.path.isdir(args.output_dir): pass

    temp_dir = tempfile.mkdtemp(dir=args.temp_dir)
    atexit.register(shutil.rmtree, temp_dir)

    jx_dump_count = 0
    junctions = {}
    max_id = -1
    for project_path in glob.glob(os.path.join(args.recount_dir, '*')):
        project = project_path.rpartition('/')[2]
        with gzip.open(
                os.path.join(project_path,
                                project
                                    + '.junction_id_with_transcripts.bed.gz')
                        ) as junction_stream:
            for line in junction_stream:
                chrom, start, end, name = line.strip().split('\t')[:4]
                name = int(name.partition('|')[0])
                max_id = max(max_id, name)
                junctions[(chrom, start, end)] = name
        if len(junctions) > args.junction_limit:
            # Dump junctions and start over; this limits size of dictionary
            with open(os.path.join(
                            temp_dir, str(jx_dump_count) + '.junc'
                        ), 'w') as junction_stream:
                for junction in junctions:
                    print >>junction_stream, '\t'.join(
                                junction + (str(junctions[junction]),)
                            )
            junctions = {}
            jx_dump_count += 1
    if junctions:
        # Dump junctions and start over; this limits size of dictionary
        with open(os.path.join(
                        temp_dir, str(jx_dump_count) + '.junc'
                    ), 'w') as junction_stream:
            for junction in junctions:
                print >>junction_stream, '\t'.join(
                            junction + (junctions[junction],)
                        )

    # Merge junction files
    os.chdir(temp_dir)
    subprocess.check_call(
            ('set -exo pipefail; '
             'cat *.junc | sort -k1,1 -k2,2n -k3,3n -u | gzip >{before_tcga}; '
             'rm -f *.junc; '
             '(gzip -cd {before_tcga}; gzip -cd {tcga_junctions}'
             ' | awk \'{{$2 -= 1; $3 -= 1; printf $1; '
                       'for (i=2;i<=NF;i++) {{printf "\t" $i}} printf "\n"}}\''
             ') | sort -k1,1 -k2,2n -k3,3n | gzip >sorted_junctions.tsv.gz'
             ).format(
                    before_tcga=os.path.join(
                        args.output_dir, 'recount_junctions_before_tcga.tsv.gz'
                    ),
                    tcga_junctions=args.tcga_junctions
                ),
            executable='/bin/bash',
            shell=True
        )

    # Get sample ID offset and write final sample ID file
    with open(
                os.path.join(args.output_dir, 'sample_ids.tsv')
            ) as output_stream:
        with open(
                os.path.join(args.recount_dir, 'sample_ids.tsv')
            ) as id_stream:
            for line in id_stream:
                print >>output_stream, line,
                current_id = line.partition('\t')[0]
        id_offset = int(current_id) + 1
        with open(args.tcga_ids) as id_stream:
            for line in id_stream:
                current_id, gdc_uuid = line.strip().split('\t')
                print >>output_stream, '\t'.join(
                        str(int(current_id) + id_offset), 'TCGA', gdc_uuid
                    )

    # Write TCGA junction files; pick up junction IDs where they were left off
    current_id = max_id
    with gzip.open(
                'sorted_junctions.tsv.gz'
            ) as junction_stream, open(
                os.path.join(
                        args.output_dir,
                        'TCGA.junction_id_with_transcripts.bed'
                    ), 'w'
            ) as bed_stream, open(
                os.path.join(
                        args.output_dir,
                        'TCGA.junction_coverage.tsv'
                    ), 'w'
            ) as coverage_stream:
        for key, group in itertools.groupby(
                junction_stream, key=lambda x: x.split('\t')[:3]
            )
            group = [line.strip().split('\t') for line in list(group)]
            if len(group) == 2:
                # There's already a junction ID
                group.sort(key=len)
                junction_id = group[0][-1]
                tokens = group[1]
            else:
                assert len(group) == 1
                current_id += 1
                junction_id = current_id
            print >>bed_stream, '\t'.join(tokens[:3] + [
                    str(junction_id), '1000', tokens[4]
                ])
            print >>coverage_stream, '\t'.join(
                    [str(junction_id)] + [
                        ','.join(str(int(sample_id) + id_offset)
                                    for sample_id in tokens[-2].split(',')),
                        tokens[-1]
                    ]
                )
