#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# based on http://stackoverflow.com/questions/4632028/how-to-create-a-temporary-directory
WORKDIR=`mktemp -d -p "$DIR"`

# deletes the temp directory
function cleanup {
  rm -rf "$WORKDIR"
  echo "Deleted temp working directory $WORKDIR"
}

# register the cleanup function to be called on the EXIT signal
trap cleanup EXIT
cat tcga_batch_17.manifest | python $DIR/true_manifest.py >$WORKDIR/tcga_batch_17.manifest
rail-rna prep elastic -m $WORKDIR/tcga_batch_17.manifest --profile dbgap --secure-stack-name dbgap --core-instance-type c3.2xlarge --master-instance-type c3.2xlarge -o s3://dummy-bucket/gtex_prep_batch_17 -c 63 --core-instance-bid-price 0.2 --master-instance-bid-price 0.2 -f --max-task-attempts 6
