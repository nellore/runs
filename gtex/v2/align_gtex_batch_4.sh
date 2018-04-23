#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna align elastic -m $DIR/gtex_batch_4.manifest --profile cwilks_dbgap --secure-stack-name dbgap-east-1c --name gtex_rail_align_4_dbgap-east-1c --core-instance-type c4.8xlarge --master-instance-type c4.8xlarge -c 160 --core-instance-bid-price 1.0 --master-instance-bid-price 1.0 -i s3://langmeadlab-private-gtex/gtex_prep_batch_4 -o s3://langmeadlab-private-gtex/gtex_align_batch_4 -a hg38 -f -d jx,tsv,bed,bw,idx --max-task-attempts 6 --use-ebs --ebs-volume-type gp2 --ebs-volumes-per-instance 6 --ebs-gb 250 --service-role EMR_dbgap --instance-profile EMR_EC2_dbgap --emr-debug
