#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna align elastic -m $DIR/tcga_batch_18.manifest --profile dbgap --secure-stack-name dbgap-us-east-1d --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --junction-criteria .01,-1 -c 80 --core-instance-bid-price 1.9 --master-instance-bid-price 1.9 -i s3://sb-rail-rna-mapreduce/tcga_prep_batch_18 -o s3://sb-rail-rna-mapreduce/tcga_align_batch_18 -a hg38 -f -d jx,tsv,bed,bw,idx --max-task-attempts 6 --name "TCGA align batch 18 job flow"
