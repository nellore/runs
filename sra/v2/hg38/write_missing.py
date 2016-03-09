#!/usr/bin/env python
"""
write_missing.py

A bug in some job flows (which used s3cmd 1.5.2) prevented some outputs from
being uploaded to S3. To assess which files were missing, while logs were still
stored on S3, we ran

for er in $(for i in {0..99}; do a=$(echo "s3://rail-sra-hg38/sra_results_batch_${i}.logs/$(aws s3 ls s3://rail-sra-hg38/sra_results_batch_${i}.logs/ | awk '{print $2}')task-attempts/"); b=$(aws s3 ls ${a} | grep "_0019" | awk '{print $2}'); c=$(aws s3 ls ${a} | grep "_0022" | awk '{print $2}'); echo $a$b; echo $a$c; done); do aws s3 cp --recursive ${er} ./; for k in $(find . -name stderr.gz 2>/dev/null); do echo $k; gzip -cd $k | grep "Skipping that file."; done; rm -rf *; done | gzip >~/res.txt.gz

in an empty directory, and then renamed the result missing_for_grep.txt.gz
which is in this directory.

This script writes which files are missing from which batches. To use it, run

gzip -cd missing_for_grep.txt.gz | pypy write_missing.py >missing.tsv

missing.tsv does not include missing files for samples that were simply
excluded from analysis. See NOTES for those.

Tab-separated output:
1. batch number
2. name of missing file
"""
import sys
import re
print 'batch\tfile'
for line in sys.stdin:
    try:
        batch_number = re.search(r'batch_(\d+)', line).groups()[0]
    except AttributeError:
        pass
    if 'ERROR' in line:
        print '\t'.join([batch_number, line.split('/')[-1].partition(':')[0]])
