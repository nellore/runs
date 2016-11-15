#!/usr/bin/env bash
curl --request POST --header "Content-Type: application/json" --data @tcga_query.json 'https://gdc-api.nci.nih.gov/v0/legacy/files' | gzip -9 >tcga_metadata.json.gz
python cgc_tcga_metadata.py Sample | gzip -9 >sample.tsv.gz
python cgc_tcga_metadata.py Case | gzip -9 >case.tsv.gz
