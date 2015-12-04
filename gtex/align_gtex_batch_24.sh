#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna align elastic -m $DIR/gtex_batch_24.manifest --profile dbgap --secure-stack-name dbgap-2 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge -c 80 --core-instance-bid-price 1.2 --master-instance-bid-price 1.2 -i s3://dbgap-stack-361204003210/gtex_prep_batch_24 -o s3://dbgap-stack-361204003210/gtex_align_batch_24 -a hg38 -f -d jx,tsv,bed,bw,idx --max-task-attempts 6
