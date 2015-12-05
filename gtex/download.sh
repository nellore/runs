#!/usr/bin/env bash
# Downloads results across all GTEx alignment batches and places them
# in subdirectories batch_k of $2 for k a batch number
# Requires AWS CLI
# $1: S3 bucket containing results; this was dbgap-stack-361204003210 for us
# $2: destination directory
# $3: AWS CLI profile to use; this was "dbgap" for us
S3BUCKET=$1
OUTPUT=$2
PROFILE=$3
mkdir -p $OUTPUT
cd $OUTPUT
for i in $(aws s3 ls $S3BUCKET --profile $PROFILE | tr -s '[:blank:]' '\t' | cut -f3 | grep align | rev | awk 'substr($1,2,1) ~ /^[0-9]+$/ {print substr($1,2)}' | rev)
do
	BATCH=$(echo $i | awk -F '_' '{print $3 "_" $4}')
	mkdir -p $BATCH
	cd $BATCH
	aws s3 cp $S3BUCKET/$i ./ --profile $PROFILE --recursive
	cd ..
done