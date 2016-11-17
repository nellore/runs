"""
merge_tables.py

Merges all CGC metadata tables output by cgc_tcga_metadata.py.
"""
import gzip
from collections import defaultdict

def notavailable(token):
    if token == 'Not available':
        return 'NA'
    return token

if __name__ == '__main__':
    (file_to_analyte, file_to_portion, file_to_case,
        file_to_aliquot, file_to_sample) = (defaultdict(lambda: 'NA'),
        defaultdict(lambda: 'NA'), defaultdict(lambda: 'NA'),
        defaultdict(lambda: 'NA'), defaultdict(lambda: 'NA'))
    with gzip.open('file.tsv.gz') as file_stream:
        file_header = ['filename'] + ['FILE_' + token for token in
                        file_stream.readline().strip().split('\t')
                        if token != 'gdcfile_uuid']
        file_to_line = defaultdict(lambda: ['NA']*(len(file_header) - 1))
        for line in file_stream:
            tokens = line.strip().split('\t')
            gdc_uuid = tokens[14]
            file_to_line[gdc_uuid] = tokens[:14] + tokens[15:]
            file_to_analyte[gdc_uuid] = tokens[-1]
            file_to_portion[gdc_uuid] = tokens[-3]
            file_to_case[gdc_uuid] = tokens[-6]
            file_to_aliquot[gdc_uuid] = tokens[-10]
            file_to_sample[gdc_uuid] = tokens[4]
    with gzip.open('sample.tsv.gz') as sample_stream:
        sample_header = ['SAMPLE_' + token for token in 
                            sample_stream.readline().strip().split('\t')[1:]]
        sample_to_line = defaultdict(lambda: ['NA']*len(sample_header))
        for line in sample_stream:
            tokens = line.strip().split('\t')
            sample_to_line[tokens[0]] = tokens[2:]
    with gzip.open('case.tsv.gz') as case_stream:
        case_header = ['CASE_' + token for token in
                            case_stream.readline().strip().split('\t')]
        case_to_line = defaultdict(lambda: ['NA']*len(case_header))
        for line in case_stream:
            tokens = line.strip().split('\t')
            case_to_line[tokens[0]] = tokens[1:]
    portion_to_slide = defaultdict(lambda: 'NA')
    with gzip.open('slide.tsv.gz') as slide_stream:
        slide_header = ['SLIDE_' + token for token in
                            slide_stream.readline().strip().split('\t')
                            if token != 'portion']
        slide_to_line = defaultdict(lambda: ['NA']*len(slide_header))
        for line in slide_stream:
            tokens = line.strip().split('\t')
            portion_to_slide[tokens[6]] = tokens[0]
            slide_to_line[tokens[0]] = tokens[1:6] + tokens[7:]
    with gzip.open('portion.tsv.gz') as portion_stream:
        portion_header = ['PORTION_' + token for token in
                            portion_stream.readline().strip().split('\t')]
        portion_to_line = defaultdict(lambda: ['NA']*len(portion_header))
        for line in portion_stream:
            tokens = line.strip().split('\t')
            portion_to_line[tokens[0]] = tokens[1:]
    case_to_drug_therapy = defaultdict(lambda: 'NA')
    with gzip.open('drug_therapy.tsv.gz') as drug_therapy_stream:
        drug_therapy_header = [
                    'DRUG_THERAPY_' + token for token in
                    drug_therapy_stream.readline().strip().split('\t')[1:]
                ]
        drug_therapy_to_line = defaultdict(lambda:
                                            ['NA']*len(drug_therapy_header))
        for line in drug_therapy_stream:
            tokens = line.strip().split('\t')
            case_to_drug_therapy[tokens[1]] = tokens[0]
            drug_therapy_to_line[tokens[0]] = tokens[2:]
    case_to_radiation_therapy = defaultdict(lambda: 'NA')
    with gzip.open('radiation_therapy.tsv.gz') as radiation_therapy_stream:
        radiation_therapy_header = [
                    'RADIATION_THERAPY_' + token for token in
                    radiation_therapy_stream.readline().strip().split('\t')[1:]
                ]
        radiation_therapy_to_line = defaultdict(lambda:
                ['NA']*len(radiation_therapy_header)
            )
        for line in radiation_therapy_stream:
            tokens = line.strip().split('\t')
            case_to_radiation_therapy[tokens[1]] = tokens[0]
            radiation_therapy_to_line[tokens[0]] = tokens[2:]
    case_to_follow_up = defaultdict(lambda: 'NA')
    with gzip.open('follow_up.tsv.gz') as follow_up_stream:
        follow_up_header = [
                    'FOLLOW_UP_' + token for token in
                    follow_up_stream.readline().strip().split('\t')[1:]
                ]
        follow_up_to_line = defaultdict(lambda: ['NA']*len(follow_up_header))
        for line in follow_up_stream:
            tokens = line.strip().split('\t')
            case_to_follow_up[tokens[1]] = tokens[0]
            follow_up_to_line[tokens[0]] = tokens[2:]
    print '\t'.join(['gdc_uuid'] + file_header + sample_header + case_header
                       + portion_header + slide_header + drug_therapy_header
                       + radiation_therapy_header + follow_up_header)
    for gdc_uuid in file_to_line:
        print '\t'.join(
                [notavailable(el) for el in [gdc_uuid] + file_to_line[gdc_uuid]
                + sample_to_line[file_to_sample[gdc_uuid]]
                + case_to_line[file_to_case[gdc_uuid]]
                + portion_to_line[file_to_portion[gdc_uuid]]
                + slide_to_line[portion_to_slide[file_to_portion[gdc_uuid]]]
                                + drug_therapy_to_line[
                        case_to_drug_therapy[file_to_case[gdc_uuid]]
                    ]
                + radiation_therapy_to_line[
                        case_to_radiation_therapy[file_to_case[gdc_uuid]]
                    ]
                + follow_up_to_line[
                        case_to_follow_up[file_to_case[gdc_uuid]]
                    ]]
            )
