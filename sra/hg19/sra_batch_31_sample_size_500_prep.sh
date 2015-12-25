#!/usr/bin/env bash
rail-rna prep elastic -m sra_batch_31_sample_size_500.txt -o s3://rail-eu-west-1/sra_batch_31_sample_size_500_prep -c 20 --region eu-west-1 --master-instance-type c3.2xlarge --core-instance-type c3.2xlarge --no-consistent-view --master-instance-bid-price 0.13 --core-instance-bid-price 0.13 --ec2-key-name raileuw1 --do-not-check-manifest --skip-bad-records --max-task-attempts 8
