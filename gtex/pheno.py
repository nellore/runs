#!/usr/bin/env bash
"""
pheno.py

Constructs GTEx table of phenotypes and sample stats (read length, etc.)
Requires SRA Tools. We used v2.5.4-1. The script also has dependencies
in the runs/gtex directory, including the manifest files, so don't move it.

Dumps pheno table to stdout.
"""
import os
from csv import reader

def blank_to_NA(it):
    """ Converts blank fields to NAs in some list it

        it: list

        Return value: list with blank fields represented as NAs
    """
    return [(el if el else 'NA') for el in it]

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--pheno-dir', type=str, required=True,
            help='path to directory containing GTEx phenotype data from '
               'dbGaP; files required are '
            'phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt, '
            'phs000424.v6.pht002741.v6.p1.GTEx_Sample.MULTI.txt, and '
            'phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt'
        )
    parser.add_argument('--fastq-dump', type=str, required=True,
            help='path to fastq-dump executable'
        )
    parser.add_argument('--secure-working-dir', type=str, required=True,
            help='path to secure working directory especially for working '
                 'with GTEx data; this is set up by vdb-config using a dbGaP '
                 'key (NGC) file as described at '
                 'https://www.youtube.com/watch?v=FjYO6Ys5cpc')
    parser.add_argument('--bigwig-out-dir', type=str, required=True,
            help='path to output directory from download.sh that contains '
                 'various batch_k files for k a batch number; this script '
                 'verifies that all sample bigwigs are present and includes '
                 'paths to sample bigwigs in the table it dumps')
    args = parser.parse_args()

    current_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(current_dir, 'SraRunInfo.csv')) as run_stream:
        run_reader = reader(run_stream)
        run_labels = run_stream.next()
        sample_to_run = {}
        run_to_sample = {}
        run_to_sra_info = {}
        for tokens in run_reader:
            if not tokens: break
            sample_to_run[tokens[25]] = tokens[0]
            run_to_sample[tokens[0]] = tokens[25]
            run_to_sra_info[tokens[0]] = tokens[1:]
    run_to_batch = {}
    run_to_bw_file = {}
    for batch_number, manifest in [
            (filename.split('_')[-1][:-9], filename)
            for filename in os.listdir(current_dir)
            if filename[-9:] == '.manifest'
        ]:
        with open(manifest) as manifest_stream:
            for line in manifest_stream:
                tokens = line.strip().split('\t')
                run = tokens[0][6:]
                bw_file = os.path.join(args.bigwig_out_dir,
                                        'batch_%s' % batch_number,
                                        'coverage_bigwigs',
                                        tokens[2] + '.bw')
                if not os.path.exists(bw_file):
                    raise RuntimeError(
                            'BigWig file {} was not found.'.format(bw_file)
                        )
                run_to_batch[run] = batch_number
                run_to_bw_file[run] = bw_file
    with open(os.path.join(args.pheno_dir,
            'phs000424.v6.pht002741.v6.p1.GTEx_Sample.MULTI.txt'
        )) as id_stream:
        for i in xrange(10): id_stream.readline()
        id_labels = id_stream.readline().strip().split('\t')
        id_labels.remove('BioSample Accession')
        sample_to_id = {}
        sample_to_SUBJID = {}
        sample_to_SAMPID = {}
        for line in id_stream:
            tokens = line.strip().split('\t')
            sample_to_id[tokens[2]] = tokens[:2] + tokens[3:]
            sample_to_SUBJID[tokens[2]] = tokens[3]
            sample_to_SAMPID[tokens[2]] = tokens[4]
    with open(os.path.join(args.pheno_dir,
            'phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt'
        )) as sample_stream:
        for i in xrange(10): sample_stream.readline()
        SAMPID_to_sample_pheno = {}
        sample_pheno_labels = sample_stream.readline().strip().split('\t')[2:]
        label_count = len(sample_pheno_labels)
        for line in sample_stream:
            if not line.strip(): continue
            tokens = line.strip().split('\t')
            SAMPID_to_sample_pheno[tokens[1]] = tokens[2:] + [''] * (
                                            label_count - len(tokens) + 2
                                        )
    with open(os.path.join(args.pheno_dir,
            'phs000424.v6.pht002742.v6.p1.c1.GTEx_Subject_Phenotypes.GRU.txt'
        )) as subject_stream:
        for i in xrange(10): subject_stream.readline()
        SUBJID_to_subject_pheno = {}
        subject_pheno_labels = subject_stream.readline().strip().split(
                                                                    '\t'
                                                                )[2:]
        label_count = len(subject_pheno_labels)
        for line in subject_stream:
            if not line.strip(): continue
            tokens = line.strip().split('\t')
            SUBJID_to_subject_pheno[tokens[1]] = tokens[2:] + [''] * (
                                            label_count - len(tokens) + 2
                                        )
    # Grab read lengths
    run_to_mate_length = {}
    '''
    for i, run in enumerate(run_to_batch):
        fastq_dump_command = (
                'set -exo pipefail; cd {secure_working_dir}; '
                '{fastq_dump} --split-spot -I --stdout -X 1 {run}'
            ).format(secure_working_dir=args.secure_working_dir,
                        fastq_dump=args.fastq_dump,
                        run=sra_run)
        single_read = subprocess.check_output(fastq_dump_command,
                                                shell=True,
                                                executable='/bin/bash').split(
                                                                        '\n'
                                                                    )
        run_to_mate_length[run] = len(single_read[2])
        print >>sys.stderr, \
            'Finished grabbing read lengths for {} samples.'.format(i + 1)'''
    for run in run_to_batch:
        run_to_mate_length[run] = 'filler'
    print '\t'.join(['Run', 'MateLength', 'RailRnaBatchNumber', 'BigWigPath']
                        + id_labels + sample_pheno_labels
                        + subject_pheno_labels)
    for run in run_to_batch:
        sample = run_to_sample[run]
        id_data = sample_to_id[sample]
        samp_data = SAMPID_to_sample_pheno[sample_to_SAMPID[sample]]
        subj_data = SUBJID_to_subject_pheno[sample_to_SUBJID[sample]]
        assert len(id_data) == len(id_labels)
        assert len(samp_data) == len(sample_pheno_labels)
        assert len(subj_data) == len(subject_pheno_labels)
        print '\t'.join([run, run_to_mate_length[run], run_to_batch[run],
                            run_to_bw_file[run]]
                    + blank_to_NA(sample_to_id[sample])
                    + blank_to_NA(
                            SAMPID_to_sample_pheno[sample_to_SAMPID[sample]]
                        )
                    + blank_to_NA(
                            SUBJID_to_subject_pheno[sample_to_SUBJID[sample]]
                        )
                )