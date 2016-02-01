#!/usr/bin/env python
"""
transform.py

Integrates metadata accumulated across files.
"""

srs_to_biosample = {}
with open('parsed_biosample_data.tsv') as parsed_stream:
    biosample_header = '\t'.join([
			'biosample_' + field.strip('"').replace(' ', '_') for field in parsed_stream.readline().strip().split('\t')
		])
    for line in parsed_stream:
        srs_to_biosample[line.strip().split('\t')[0]] = line.strip()

srr_to_srs = {}
with open('index_to_SRA_accession.tsv') as index_stream:
    for line in index_stream:
        line = line.strip().split('\t')
        srr_to_srs[line[-1]] = line[2]

srr_to_gsm = {}
with open('SRRToGSM.tsv') as index_stream:
    for line in index_stream:
        line = line.strip().split('\t')
        srr_to_gsm[line[0]] = line[1]

srr_to_geo_raw = {}
with open('geo_pheno.csv') as geo_stream_c, \
	open('geo_pheno.txt') as geo_stream_t: # From Jaffe
	geo_raw_header = '\t'.join([
			'geo_raw_' + field.strip('"').replace(' ', '_') for field in geo_stream_c.readline().strip().split(',')[1:]
		])
	geo_stream_t.readline()
	geo_raw_len = geo_raw_header.count('\t') + 1
	for line in geo_stream_c:
		line2 = geo_stream_t.readline()
		real_line = [el if el == line2[i] else '\x1f' for i, el in enumerate(line)]
		tokens = line.strip().split('\x1f')
		srr = tokens[0].strip('"')
		srr_to_geo_raw[srr] = '\t'.join([token.strip('"').replace('\t', ' ') for token in tokens[1:]])

gsm_to_geo_processed = {} 
with open('importantGEOmetadata.tsv') as geo_stream: # From Jaffe's soft files, processed with extract_geo.py
	geo_proc_header = '\t'.join([
			'geo_processed_' + field.strip('"').replace(' ', '_') for field in geo_stream.readline().strip().split('\t')
		])
	geo_proc_len = geo_proc_header.count('\t') + 1
	for line in geo_stream:
		tokens = line.strip().split('\t')
		tokens[0] = tokens[0][:-5]
		gsm_to_geo_processed[tokens[0]] = '\t'.join([token.strip('"') for token in tokens])

srr_to_sra = {}
with open('all_illumina_sra_for_human.txt') as sra_stream:
	sra_header = '\t'.join(['sra_' + token.replace(' ', '_') for token in sra_stream.readline().strip().split('\t')[1:]])
	for line in sra_stream:
		tokens = line.strip().split('\t')
		srr_to_sra[tokens[0]] = '\t'.join(tokens[1:])

print '\t'.join(['run_accession', biosample_header, geo_proc_header, geo_raw_header, sra_header])
for srr in srr_to_srs:
    print '\t'.join([srr, srs_to_biosample[srr_to_srs[srr]],
    	gsm_to_geo_processed[srr_to_gsm[srr]] if srr in srr_to_gsm else '\t'.join(['NA']*geo_proc_len),
    	srr_to_geo_raw[srr] if srr in srr_to_geo_raw else '\t'.join(['NA']*geo_raw_len),
    	srr_to_sra[srr]])
