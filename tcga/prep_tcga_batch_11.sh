#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# based on http://stackoverflow.com/questions/4632028/how-to-create-a-temporary-directory
WORKDIR=$(mktemp -d)

# deletes the temp directory
function cleanup {
  rm -rf $WORKDIR
  echo "Deleted temp working directory $WORKDIR"
}

# register the cleanup function to be called on the EXIT signal
trap cleanup EXIT
cat tcga_batch_11.manifest | python $DIR/true_manifest.py --cgc-auth-token /Users/eterna/cgcauth.txt >$WORKDIR/tcga_batch_11.manifest
rail-rna prep elastic -m $WORKDIR/tcga_batch_11.manifest --profile dbgap --secure-stack-name dbgap-us-east-1e --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge -o s3://sb-rail-rna-mapreduce/tcga_prep_batch_11 -c 63 --core-instance-bid-price 0.9 --master-instance-bid-price 0.9 -f --max-task-attempts 6 --skip-bad-records
