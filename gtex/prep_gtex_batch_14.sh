#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
rail-rna prep elastic -m $DIR/gtex_batch_14.manifest --profile dbgap --secure-stack-name dbgap --dbgap-key /Users/eterna/gtex/prj_8716.ngc --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge -o s3://dbgap-stack-361204003210/gtex_prep_batch_14 -c 63 --core-instance-bid-price 0.25 --master-instance-bid-price 0.25
