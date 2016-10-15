#!/usr/bin/env python
"""
gen.py

Uses tcga_file_list.tsv to creates dummy manifest files and scripts for running
Rail-RNA on TCGA RNA-seq data hosted on S3 courtesy of Seven Bridges Genomics.
tcga_file_list.tsv was obtained by running tcga_file_list.py in this
directory; see README.md for details. Preprocess scripts use
true_manifest.py, which reads a dummy manifest file with CGC storage paths and
uses the CGC API to replace them with S3 URLs.

We ran

python gen.py --s3-bucket s3://sb-rail-rna-mapreduce --region us-east-1
    --c3-2xlarge-bid-price 0.50 --c3-8xlarge-bid-price 1.90
    --prep-stack-names dbgap-us-east-1a dbgap-us-east-1b dbgap-us-east-1d
    dbgap-us-east-1e --align-stack-names dbgap-us-east-1a dbgap-us-east-1b
    dbgap-us-east-1d dbgap-us-east-1e
    --cgc-auth-token /path/to/cgc_authorization_token.txt

and used Rail-RNA v0.2.4a. Note that /path/to/cgc_authorization_token.txt is
the path to a text file with a single line of text: the CGC authorization token
which we obtained at https://cgc.sbgenomics.com/account/#developer. After
the manifests, preprocess, and align scripts were generated, we modified some
so we could launch into VPC public subnets in different availability zones.
We created a total of five VPCs using the instructions at
http://docs.rail.bio/dbgap/ . To simulate our setup, first use the
CloudFormation template.

https://github.com/nellore/rail/blob/v0.2.1/src/cloudformation/dbgap.template

and name the stack "dbgap". Then create another four VPCs in different
availability zones in us-east-1 using

https://github.com/nellore/rail/blob/v0.2.1/src/cloudformation/
    dbgap_minus_cloudtrail.template

, and name them "dbgap-us-east-1a", "dbgap-us-east-1b", "dbgap-us-east-1c", and
"dbgap-us-east-1d". Each subnet accommodates up to ~65k IPs and is associated
with a given availability zone. We executed each prep_tcga_batch_<index>.sh
script and waited for the job flow to finish before executing the corresponding
align_tcga_batch_<index>.sh script. Sometimes, we changed the availability
zone to which we submitted a job flow to minimize cost. (At a given time, the 
market price in one availability zone may be lower than in another.)
"""
import random
import sys
import os
from itertools import cycle
import re

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--s3-bucket', type=str, required=True,
            help=('path to S3 bucket in which preprocessed data and junction '
                  'data will be dumped; should be a secure bucket created '
                  'by following the instructions at '
                  'http://docs.rail.bio/dbgap/; ours was s3://rail-dbgap')
        )
    parser.add_argument('--region', type=str, required=True,
            help='AWS region in which to run job flows; we used us-east-1'
        )
    parser.add_argument('--c3-2xlarge-bid-price', type=float, required=False,
            default=0.20,
            help='bid price for each c3.xlarge instance; this instance '
                 'type is used for preprocessing data'
        )
    parser.add_argument('--c3-8xlarge-bid-price', type=float, required=False,
            default=1.20,
            help='bid price for each c3.8xlarge instance; this instance '
                 'type is used for aligning data'
        )
    parser.add_argument('--prep-stack-names', type=str, required=False,
            default='dbgap', nargs='+',
            help='stack name(s) for prep job flow; cycle through them'
        )
    parser.add_argument('--align-stack-names', type=str, required=False,
            default='dbgap', nargs='+',
            help='stack name(s) for align job flow; cycle through them'
        )
    parser.add_argument('--seed', type=int, required=False,
            default=78943,
            help=('seed for random number generator; random.shuffle is used '
                  'to shuffle the TCGA samples before dividing them up into '
                  '--batch-count batches')
        )
    parser.add_argument('--tcga-file-list', type=str, required=False,
            default=os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    'tcga_file_list.tsv'
                ),
            help='path to tcga_file_list.tsv generated by tcga_file_list.py'
        )
    parser.add_argument('--cgc-auth-token', type=str, required=True,
            help='path to text file containing CGC authorization token; '
        )
    parser.add_argument('--batch-count', type=int, required=False,
            default=30,
            help='number of batches to create; batches are designed to be '
                 'of approximately equal size'
        )
    args = parser.parse_args()
    manifest_lines = []
    with open(args.tcga_file_list) as tcga_file_list_stream:
        tcga_file_list_stream.readline() # header line
        for line in tcga_file_list_stream:
            tokens = line.strip().split('\t')
            manifest_lines.append('\t'.join(
                    [tokens[1], '0', tokens[0]] # use GDC UUID as sample name
                ))
    random.seed(args.seed)
    random.shuffle(manifest_lines)
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # Write all manifest files
    manifest_files = [[] for i in xrange(args.batch_count)]
    for i, manifest_index in enumerate(cycle(range(args.batch_count))):
        try:
            manifest_files[manifest_index].append(manifest_lines[i])
        except IndexError:
            # No more manifest lines
            break
    for i, manifest_file in enumerate(manifest_files):
        with open('tcga_batch_{}.manifest'.format(i), 'w') as manifest_stream:
            for line in manifest_file:
                print >>manifest_stream, line
    # Write all prep and align scripts
    prep_stack_name_cycle = cycle(args.prep_stack_names)
    align_stack_name_cycle = cycle(args.align_stack_names)
    for i in xrange(args.batch_count):
        with open('prep_tcga_batch_{}.sh'.format(i), 'w') as prep_stream:
            print >>prep_stream, (
"""#!/usr/bin/env bash
DIR="$( cd "$( dirname "${{BASH_SOURCE[0]}}" )" && pwd )"
# based on http://stackoverflow.com/questions/4632028/how-to-create-a-temporary-directory
WORKDIR=$(mktemp -d)

# deletes the temp directory
function cleanup {{
  rm -rf $WORKDIR
  echo "Deleted temp working directory $WORKDIR"
}}

# register the cleanup function to be called on the EXIT signal
trap cleanup EXIT
cat {manifest_file} | python $DIR/true_manifest.py \
--cgc-auth-token {cgc_auth_token} >$WORKDIR/{manifest_file}
rail-rna prep elastic -m $WORKDIR/{manifest_file} --profile dbgap \
--secure-stack-name {stack_name} \
--core-instance-type c3.2xlarge --master-instance-type c3.2xlarge \
-o {s3_bucket}/tcga_prep_batch_{batch_number} \
-c 48 --core-instance-bid-price {core_price} \
--master-instance-bid-price {core_price} -f \
--max-task-attempts 10 --skip-bad-records --do-not-check-manifest \
--name TCGA_prep_batch_{batch_number}_job_flow
"""
                ).format(manifest_file='tcga_batch_{}.manifest'.format(i),
                            s3_bucket=args.s3_bucket,
                            batch_number=i,
                            core_price=args.c3_2xlarge_bid_price,
                            stack_name=next(prep_stack_name_cycle),
                            cgc_auth_token=args.cgc_auth_token)
        with open('align_tcga_batch_{}.sh'.format(i), 'w') as align_stream:
            print >>align_stream, '#!/usr/bin/env bash'
            print >>align_stream, (
                    'DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"'
                )
            print >>align_stream, (
                    'rail-rna align elastic -m $DIR/{manifest_file} '
                    '--profile dbgap --secure-stack-name {stack_name} '
                    '--core-instance-type c3.8xlarge '
                    '--master-instance-type c3.8xlarge '
                    '--junction-criteria .01,-1 '
                    '-c 60 --core-instance-bid-price {core_price} '
                    '--master-instance-bid-price {core_price} '
                    '-i {s3_bucket}/tcga_prep_batch_{batch_number} '
                    '-o {s3_bucket}/tcga_align_batch_{batch_number} '
                    '-a hg38 -f -d jx,tsv,bed,bw,idx '
                    '--max-task-attempts 6 '
                    '--name TCGA_align_batch_{batch_number}_job_flow'
                ).format(manifest_file='tcga_batch_{}.manifest'.format(i),
                            s3_bucket=args.s3_bucket,
                            batch_number=i,
                            core_price=args.c3_8xlarge_bid_price,
                            stack_name=next(align_stack_name_cycle))
