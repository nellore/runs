"""
add_knowngene.py

Appends transcript IDs to name field of a junction ID BED file output by
junctions_by_project.py: one is a comma-separated list of IDs of transcripts in
which the junction is found, and the other is a comma-separated list of
corresponding gene IDs.

Requires that GTF has been downloaded via the UCSC Table Browser at 
j ; we downloaded this GTF
on June 13, 2016 as knownGene_hg38.gtf.gz, which is in the same dir as this
script. ALSO CORRECTS JUNCTION COORDINATES FROM junctions_by_project.py so
they're zero-based!
"""
from collections import defaultdict
import gzip

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtf', required=True,
        help='path to UCSC knownGene GTF.gz file for hg38')
    parser.add_argument('--beds', required=True, nargs='+',
        help='paths to BED files to process')
    args = parser.parse_args()
    # Get junctions and associated gene/transcript IDs
    exons = defaultdict(set)
    with gzip.open(args.gtf) as gtf_stream:
        for line in gtf_stream:
            if line[0] == '#': continue
            tokens = line.strip().split('\t')
            if tokens[2].lower() != 'exon': continue
            '''key: transcript_id
               value: (rname, exon start (1-based), exon end (1-based))

            transcript_id in token 12 is decked with " on the left and "; on
            the right; kill them in key below.
            '''
            attribute = tokens[-1].split(';')
            id_index = [i for i, name in enumerate(attribute)
                        if 'transcript_id' in name]
            assert len(id_index) == 1, ('More than one transcript ID specified; '
                                        'offending line: %s') % line 
            id_index = id_index[0]
            attribute[id_index] = attribute[id_index].strip()
            quote_index = attribute[id_index].index('"')
            exons[attribute[id_index][quote_index+1:-1]].add(
                    (tokens[0], int(tokens[3]), int(tokens[4]), tokens[6])
                )

    junctions = defaultdict(list)
    donors = defaultdict(list)
    acceptors = defaultdict(list)
    for transcript_id in exons:
        exons_from_transcript = sorted(list(exons[transcript_id]))
        strand = exons_from_transcript[0][3]
        # Recall that GTF is end-inclusive, and so is STAR's junctions.txt
        for i in xrange(1, len(exons_from_transcript)):
            if exons_from_transcript[i][0] == exons_from_transcript[i-1][0]:
                # Kill any introns 3 bases or smaller
                if (exons_from_transcript[i][1]
                    - exons_from_transcript[i-1][2]) < 5:
                    continue
                junction = (exons_from_transcript[i][0],
                            exons_from_transcript[i-1][2],
                            exons_from_transcript[i][1] - 2)
                donor = (exons_from_transcript[i][0],
                         exons_from_transcript[i-1][2])
                acceptor = (exons_from_transcript[i][0],
                            exons_from_transcript[i][1] - 2)
                if strand == '-':
                    donor, acceptor = acceptor, donor
                donors[donor].append(transcript_id)
                acceptors[acceptor].append(transcript_id)
                junctions[junction].append(transcript_id)

    # Postprocess
    junctions_write = defaultdict(lambda: 'NA')
    donors_write = defaultdict(lambda: 'NA')
    acceptors_write = defaultdict(lambda: 'NA')
    for junction in junctions:
        junctions_write[junction] = ';'.join(junctions[junction])
    for donor in donors:
        donors_write[donor] = ';'.join(donors[donor])
    for acceptor in acceptors:
        acceptors_write[acceptor] = ';'.join(acceptors[acceptor])
    # Modify BED files
    for bed in args.beds:
        assert bed.endswith('.bed.gz'), (
                '{} does not have extension ".bed.gz".'.format(bed)
            )
        modified_bed_filename = bed[:-7] + '_with_transcripts.bed.gz'
        with gzip.open(modified_bed_filename, 'w') as output_stream, gzip.open(
                bed
            ) as input_stream:
            for line in input_stream:
                tokens = line.strip().split('\t')
                junction = (tokens[0], int(tokens[1]) - 1, int(tokens[2]) - 1)
                strand = tokens[5]
                donor = junction[:2]
                acceptor = (junction[0], junction[2])
                if strand == '-':
                    donor, acceptor = acceptor, donor
                print >>output_stream, '\t'.join(
                            list(map(str, junction))
                            + [''.join([tokens[3],
                                        '|D:',
                                        donors_write[donor],
                                        '|A:',
                                        acceptors_write[acceptor],
                                        '|J:',
                                        junctions_write[junction]])]
                            + tokens[4:]
                        )
