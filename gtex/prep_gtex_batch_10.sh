#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna prep elastic -m $DIR/gtex_batch_10.manifest --profile dbgap --secure-stack-name dbgap --dbgap-key /Users/eterna/gtex/prj_8716.ngc --core-instance-type m3.xlarge --master-instance-type m3.xlarge -o s3://dbgap-stack-361204003210/gtex_prep_batch_10 -c 20 --core-instance-bid-price 0.25 --master-instance-bid-price 0.25 -f --max-task-attempts 6
