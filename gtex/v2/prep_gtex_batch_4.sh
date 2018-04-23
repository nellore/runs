#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna prep elastic -m $DIR/gtex_batch_4.manifest --profile cwilks_dbgap --secure-stack-name dbgap --name gtex_rail_prep_4_dbgap --dbgap-key /home/cwilks/.dbgap/gtex_prj_8716.ngc --core-instance-type m3.xlarge --master-instance-type m3.xlarge -o s3://langmeadlab-private-gtex/gtex_prep_batch_4 -c 20 --core-instance-bid-price 1.0 --master-instance-bid-price 1.0 -f --max-task-attempts 6 --service-role EMR_dbgap --instance-profile EMR_EC2_dbgap --emr-debug
