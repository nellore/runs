#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna align elastic -m $DIR/sra_batch_89.manifest --core-instance-type c3.8xlarge --master-instance-type c3.8xlarge --junction-criteria .01,-1 -c 80 -e --core-instance-bid-price 1.7 --master-instance-bid-price 1.7 -i s3://rail-sra-hg38/sra_prep_batch_89 -o s3://rail-sra-hg38/sra_results_batch_89 -a hg38 -f -d jx,tsv,bed,bw,idx --max-task-attempts 6 --region us-east-1 --ec2-key-name useast1
