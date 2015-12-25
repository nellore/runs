#!/usr/bin/env python
"""
create_runs.py

Creates scripts for running on human Illumina RNA-seq from SRA in batches on
EMR. Use Rail-RNA v0.1.7a to reproduce results.
Usage: cat human_illumina_sra.txt | python create_batches.py [options]
"""
import random
import sys
import os

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--s3-bucket', type=str, required=True,
            help=('path to S3 bucket in which preprocessed data and junction '
                  'data will be dumped')
        )
    parser.add_argument('--region', type=str, required=True,
            help='AWS region in which to run job flows'
        )
    parser.add_argument('--c3-2xlarge-bid-price', type=float, required=False,
            default=0.13,
            help='bid price for each c3.2xlarge instance; this instance '
                 'type is used for preprocessing data'
        )
    parser.add_argument('--c3-8xlarge-bid-price', type=float, required=False,
            default=0.33,
            help='bid price for each c3.8xlarge instance; this instance '
                 'type is used for aligning data'
        )
    args = parser.parse_args()
    lines = sys.stdin.read().strip().split('\n')
    random.seed(lines[0])
    random.shuffle(lines)
    manifests = [lines[i:i+500] for i in xrange(0, len(lines), 500)]
    # Merge last two batches
    manifests[-2] += manifests[-1]
    manifests.pop()
    for i, manifest in enumerate(manifests):
        manifest_filename = 'sra_batch_{}_sample_size_{}.txt'.format(
                                                                i, len(manifest)
                                                            )
        prep_output_dir = os.path.join(
                            args.s3_bucket,
                            'sra_batch_{}_sample_size_{}_prep'.format(
                                                            i, len(manifest)
                                                        )
                        )
        align_output_dir = os.path.join(
                            args.s3_bucket,
                            'sra_batch_{}_sample_size_{}_itn'.format(
                                                            i, len(manifest)
                                                        )
                        )
        with open('sra_batch_{}_sample_size_{}_prep.sh'.format(
                                                            i, len(manifest)
                                                        ), 'w') as prep_stream:
            prep_stream.write(
    """#!/usr/bin/env bash
    rail-rna prep elastic -m {manifest} -o {output_dir} -c 20 --region {region} --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge --no-consistent-view --master-instance-bid-price {bid_price} --core-instance-bid-price {bid_price} --do-not-check-manifest --skip-bad-records --max-task-attempts 8
    """.format(manifest=manifest_filename, output_dir=prep_output_dir, region=args.region, 
                bid_price=args.c3_2xlarge_bid_price)
    )
        with open('sra_batch_{}_sample_size_{}_itn.sh'.format(
                                                            i, len(manifest)
                                                        ), 'w') as align_stream:
            align_stream.write(
    """#!/usr/bin/env bash
    rail-rna align elastic -m {manifest} -i {input_dir} -o {output_dir} -a hg19 --region {region} -c 60 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --core-instance-bid-price {bid_price} --master-instance-bid-price {bid_price} --no-consistent-view --deliverables itn --max-task-attempts 6
    """.format(manifest=(manifest_filename if i != 9 else 
                            'sra_batch_9_sample_size_500_old.txt'),
                input_dir=prep_output_dir,
                output_dir=align_output_dir, region=args.region, 
                bid_price=args.c3_8xlarge_bid_price)
    )
