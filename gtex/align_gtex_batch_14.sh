#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna align elastic -m $DIR/gtex_batch_14.manifest --profile dbgap --secure-stack-name dbgap -o s3://dbgap-stack-361204003210 --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge -c 80 --core-instance-bid-price 1.2 --master-instance-bid-price 1.2 -i s3://dbgap-stack-361204003210/gtex_prep_batch_14 -o s3://dbgap-stack-361204003210/gtex_align_batch_14 -a hg38
