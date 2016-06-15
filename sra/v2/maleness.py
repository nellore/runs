#!/usr/bin/env python
"""
maleness.py

See arguments for requirements. Output:
1. project accession number
2. sample accession number
3. experiment accession number
4. run accession number
5. number of _annotated_ chrY junctions expressed
"""
from collections import defaultdict
import gzip

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--intropolis', required=True,
        help='path to intropolis.v2.hg38.tsv.gz, which contains all '
             'junctions from the second run of Rail')
    parser.add_argument('--ids', required=True,
        help='path to intropolis.idmap.v2.hg38.tsv, which maps sample ids '
             'from intropolis to SRA accession numbers')
    parser.add_argument('--annotated', required=True,
        help='path to annotated junctions; this is annotated_junctions.tsv.gz '
             'and covers all annotations defined in annotation_definition.md')
    args = parser.parse_args()

    with gzip.open(args.annotated) as annotated_stream:
        annotated_junctions = [
                line.strip().split('\t') for line in annotated_stream
            ]

    annotated_chrY_junctions = set(
            [tuple(junction[1:]) for junction in annotated_junctions
                if junction[0] == 'chrY']
        )

    male_samples = defaultdict(int)
    with gzip.open(args.intropolis) as intropolis_stream:
        for line in intropolis_stream:
            (chrom, pos, end_pos, strand,
                _, _, samples, _) = line.strip().split('\t')
            if chrom != 'chrY': continue
            if (pos, end_pos, strand) in annotated_chrY_junctions:
                for sample in samples.split(','):
                    male_samples[sample] += 1

    id_to_accession = {}
    with open(args.ids) as id_stream:
        for line in id_stream:
            tokens = line.strip().split('\t')
            id_to_accession[tokens[0]] = '\t'.join(tokens[1:])

    for sample, count in sorted(male_samples.items(),
                                key=lambda x: x[1],
                                reverse=True):
        print '\t'.join([id_to_accession[sample], str(count)])
