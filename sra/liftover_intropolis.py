#!/usr/bin/env python
"""
liftover_intropolis.py

Lifts intropolis over to hg38 and writes output including both hg19 coordinates
and hg38 coordinates.

Requires 

http://hgdownload.cse.ucsc.edu/goldenpath/mm10/
    liftOver/hg19ToHg38.over.chain.gz

liftOver executable available from
    https://genome-store.ucsc.edu/products/

intropolis.v1.hg19.tsv.gz .

Writes to stdout. We ran

pypy liftover_intropolis.py
    --liftover /path/to/liftOver
    --chain /path/to/hg19ToHg38.over.chain
    --intropolis /path/to/intropolis.v1.hg19.tsv.gz
    | gzip >intropolis.v1.hg19_with_hg38_liftover.tsv.gz

Tab-separated output fields
1. hg19 chrom
2. hg19 start (1-based, inclusive)
3. hg19 end (1-based, inclusive)
4. hg19 strand
5. left motif (e.g., GT)
6. right motif (e.g., AG)
7. comma-separated list of indexes of samples from
    intropolis in which junction was found
8. comma-separated list of numbers of reads in corresponding samples from 
    field 7 overlapping junction
9. hg38 chrom or NA if liftover unsuccessful
10. hg38 start or NA
11. hg38 end or NA
12. hg38 strand or NA
"""
import tempfile
import gzip
import shutil
import atexit
import subprocess

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser.add_argument('--liftover', type=str, required=True,
            help=('path to liftOver executable available from '
                  'https://genome-store.ucsc.edu/products/')
        )
    parser.add_argument('--chain', type=str, required=True,
            help=('path to unzipped liftover chain; this should be '
                  'hg19ToHg38.over.chain')
        )
    parser.add_argument('--intropolis', type=str, required=True,
            help='path to intropolis.v1.hg19.tsv.gz'
        )
    parser.add_argument('--temp-dir', type=str, required=False,
            default=None,
            help='where to store temporary files; defaults to TMPDIR'
        )
    args = parser.parse_args()
    temp_dir = tempfile.mkdtemp(dir=args.temp_dir)
    to_liftover = os.path.join(args.temp_dir, 'to_liftover.bed')
    temp_hg38 = os.path.join(temp_dir, 'hg38.bed')
    temp_hg19 = os.path.join(temp_dir, 'hg19.bed')
    with open(temp_hg19, 'w') as hg19_stream, gzip.open(
            args.intropolis
        ) as input_stream:
        for i, line in enumerate(input_stream):
            tokens = line.strip().split('\t')
            chrom, strand, start, end = (
                    tokens[0][:-1], tokens[0][-1], str(int(tokens[1]) - 1),
                    tokens[2]
                ) # zero-based, half-open coordinates for BED
            # Tack original junction onto junction name
            junction_name = ';'.join([str(i), chrom, start, end, strand])
            print >>hg19_stream, '{}\t{}\t{}\tinfo_{}\t1\t{}'.format(
                    chrom, start, end, junction_name, strand
                )
    # Convert junctions from hg19 to hg38
    liftover_process = subprocess.call(' '.join([
                                            args.liftover,
                                            temp_hg19,
                                            args.chain,
                                            temp_hg38,
                                            os.path.join(
                                                    temp_dir,
                                                    'unmapped.bed'
                                                )
                                        ]),
                                        shell=True,
                                        executable='/bin/bash'
                                    )
    to_sort = os.path.join(temp_dir, 'intropolis_and_liftover.tsv.gz')
    with gzip.open(to_sort, 'w') as both_stream:
        with open(temp_hg38) as hg38_stream:
            for line in hg38_stream:
                chrom, start, end, name, score, strand = line.strip().split(
                                                                        '\t'
                                                                    )[:6]
                (_, hg19_chrom, hg19_start,
                        hg19_end, hg19_strand) = name.split(';')
                hg19_start = int(hg19_start), int(start)
                print >>both_stream, '\t'.join(
                                [hg19_chrom, str(hg19_start + 1), hg19_end,
                                    hg19_strand, chrom, str(start + 1),
                                    end, strand, 'FAKE']
                            )
        with gzip.open(args.intropolis) as intropolis_stream:
            for line in intropolis_stream:
                print >>both_stream, line,
    sorted_together = os.path.join(temp_dir, 'sorted_together.tsv.gz')
    subprocess.check_call(
            'gzip -cd {} | sort -k1,1 -k2,2n -k3,3n | gzip >{}'.format(
                    to_sort, sorted_together
                ), shell=True, bufsize=-1
        )
    with gzip.open(sorted_together) as sorted_stream:
        junction_1_tokens = sorted_stream.readline().strip().split('\t')
        junction_2_tokens = sorted_stream.readline().strip().split('\t')
        while True:
            if junction_1_tokens[:4] == junction_2_tokens[:4]:
                # Liftover available
                if len(junction_1_tokens) > len(junction_2_tokens):
                    hg38_tokens = junction_1_tokens
                    hg19_tokens = junction_2_tokens
                else:
                    hg38_tokens = junction_2_tokens
                    hg19_tokens = junction_1_tokens
                print '\t'.join(hg19_tokens + hg38_tokens[:4])
                junction_1_tokens = sorted_stream.readline().strip().split(
                                                                        '\t'
                                                                    )
                junction_2_tokens = sorted_stream.readline().strip().split(
                                                                        '\t'
                                                                    )
            else:
                '''Liftover not available for junction 1, but have to check
                junction 2 against next junction.'''
                print '\t'.join(junction_1_tokens + ['NA'] * 4)
                junction_1_tokens = junction_2_tokens
                junction_2_tokens = sorted_stream.readline().strip()
                if not junction_2_tokens:
                    # End of file; print junction 1 tokens and sign out
                    print '\t'.join(junction_1_tokens + ['NA'] * 4)
                    break
                junction_2_tokens = junction_2_tokens.split('\t')
