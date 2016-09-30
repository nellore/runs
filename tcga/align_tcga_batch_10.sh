#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna align elastic -m $DIR/tcga_batch_10.manifest --profile dbgap --secure-stack-name dbgap --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge -c 80 --core-instance-bid-price 1.2 --master-instance-bid-price 1.2 -i s3://dummy-bucket/tcga_prep_batch_10 -o s3://dummy-bucket/tcga_align_batch_10 -a hg38 -f -d jx,tsv,bed,bw,idx --max-task-attempts 6
