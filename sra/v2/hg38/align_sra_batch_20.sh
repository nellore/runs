#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna align elastic -m $DIR/sra_batch_20.manifest --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge -c 80 --core-instance-bid-price 1.7 --master-instance-bid-price 1.7 -i s3://rail-sra-hg38/sra_prep_batch_20 -o s3://rail-sra-hg38/sra_align_batch_20 -a hg38 -f -d jx,tsv,bed,bw,idx --max-task-attempts 6 --region us-east-1
