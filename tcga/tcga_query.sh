#!/usr/bin/env bash
curl --request POST --header "Content-Type: application/json" --data @tcga_query.json 'https://gdc-api.nci.nih.gov/v0/legacy/files' | gzip -9 >tcga_metadata.json.gz
python cgc_tcga_metadata.py Sample | gzip -9 >sample.tsv.gz
python cgc_tcga_metadata.py Case | gzip -9 >case.tsv.gz
python cgc_tcga_metadata.py Slide | gzip -9 >slide.tsv.gz
python cgc_tcga_metadata.py Analyte | gzip -9 >analyte.tsv.gz
python cgc_tcga_metadata.py Portion | gzip -9 >portion.tsv.gz
python cgc_tcga_metadata.py Aliquot | gzip -9 >aliquot.tsv.gz
python cgc_tcga_metadata.py DrugTherapy | gzip -9 >drug_therapy.tsv.gz
python cgc_tcga_metadata.py RadiationTherapy | gzip -9 >radiation_therapy.tsv.gz
python cgc_tcga_metadata.py FollowUp | gzip -9 >follow_up.tsv.gz
python cgc_tcga_metadata.py File | gzip -9 >file.tsv.gz
