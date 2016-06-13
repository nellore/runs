"""
add_knowngene.py

Appends transcript IDs to name field of a junction ID BED file output by
junctions_by_project.py: one is a comma-separated list of IDs of transcripts in
which the junction is found, and the other is a comma-separated list of
corresponding gene IDs.

Requires that knownGene.gtf has been downloaded via the UCSC Table Browser at 
http://genome.ucsc.edu/cgi-bin/hgTables?command=start .

Reads junction ID BEDs from 
"""
from collections import defaultdict
import sys

if __name__ == '__main__':
	import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gtf', required=True,
        help='path to UCSC knownGene GTF file for hg38')
    parser.add_argument('--beds', required=True, nargs='+',
        help='paths to BED files to process')
    args = parser.parse_args()
    # Get junctions and associated gene/transcript IDs
    with open(args.gtf) as gtf_stream:
    	exons = defaultdict(set)
    for line in sys.stdin:
        if line[0] == '#': continue
        tokens = line.strip().split('\t')
        if tokens[2].lower() != 'exon': continue
        '''key: transcript_id
           value: (rname, exon start (1-based), exon end (1-based))

        transcript_id in token 12 is decked with " on the left and "; on the
        right; kill them in key below.
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
                (tokens[0], int(tokens[3]), int(tokens[4]))
            )

    junctions = defaultdict(list)
    for transcript_id in exons:
        exons_from_transcript = sorted(list(exons[transcript_id]))
        # Recall that GTF is end-inclusive, and so is STAR's junctions.txt
        for i in xrange(1, len(exons_from_transcript)):
            if exons_from_transcript[i][0] == exons_from_transcript[i-1][0]:
                # Kill any introns 4 bases or smaller
                if (exons_from_transcript[i][1]
                    - exons_from_transcript[i-1][2]) < 5:
                    continue
                junction = (exons_from_transcript[i][0],
                            exons_from_transcript[i-1][2] + 1,
                            exons_from_transcript[i][1] - 1)
                junctions[junction].append(transcript_id)

    # Modify BED files
    for bed in args.beds:
        assert bed.endswith('.bed'), (
                '{} does not have extension ".bed".'.format(bed)
            )
        modified_bed_filename = bed[:-4] + '_with_transcript.bed'
        with open(modified_bed_filename, 'w') as output_stream, open(
                bed
            ) as input_stream:
            for line in input_stream:
                tokens = line.strip().split('\t')
                junction = (tokens[0], int(tokens[1]), int(tokens[2]))
                if junction in junctions:
                    print >>output_stream, '\t'.join(
                            tokens[:3]
                            + [tokens[4] + ':' + ';'.join(
                                    [str(el) for el in junctions[junction]]
                                )] + tokens[5:]
                        )
                else:
                    print >>output_stream, '\t'.join(
                            tokens[:3] + [tokens[4] + ':NA'] + tokens[5:]
                        )
