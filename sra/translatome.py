#!/usr/bin/env python
"""
translatome.py

Performs cross-species liftover of junctions from translatome study SRP031883
obtained by Rail-RNA in translatome.sh from mm10 to hg19
(stored in mm10_translatome_junctions.tsv.gz); then identifies
which of these junctions are in intropolis.v1.hg19.tsv.gz but aren't annotated
(i.e., in annotated_junctions.tsv.gz).

Requires

http://hgdownload.cse.ucsc.edu/goldenpath/mm10/
    liftOver/mm10ToHg19.over.chain.gz

liftOver executable available from
    https://genome-store.ucsc.edu/products/,

intropolis.v1.hg19.tsv.gz

mm10_translatome_junctions.tsv.gz (in this directory)

and annotated_junctions.tsv.gz (in this directory).

Stats are written to stderr; we store them in
    translatome_stats.txt. We store mm10 regions that do not
    map to hg19 in unmapped_mm10.bed .

From the runs/sra directory, we ran

pypy translatome.py
    --chain /path/to/hg38ToHg19.over.chain
    --liftover /path/to/liftOver
    --unmapped unmapped_hg38.bed 2>translatome_stats.txt
    | sort -k1,1 -k2,2n -k3,3n | gzip
    >translatome_mm10_to_hg19_junctions.tsv.gz

Tab-separated output is unannotated lifted-over junctions in intropolis:
1. lifted-over hg19 chrom
2. lifted-over hg19 start position (1-based, inclusive)
3. lifted-over hg19 end position (1-based, inclusive)
4. lifted-over strand
5. lifted-over 5' motif (e.g., GT)
6. lifted-over 3' motif (e.g., AG)
7. comma-separated list of indexes of samples from
    intropolis in which junction was found
8. comma-separated list of numbers of reads in corresponding samples from 
    field 7 overlapping junction
9. mm10 chrom
10. mm10 start position
11. mm10 end position
12. mm10 strand
13. comma-separated list of indexes of translatome samples in which junction
    was found; see translatome.manifest for these samples in order, indexed
    from 0 thru 13
14. comma-separated list of numbers of reads in corresponding samples from
    field 9 overlapping junction
"""
import gzip
import shutil
import atexit
import subprocess

if __name__ == '__main__':
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--liftover', type=str, required=True,
            help=('path to liftOver executable available from '
                  'https://genome-store.ucsc.edu/products/')
        )
    parser.add_argument('--chain', type=str, required=True,
            help=('path to unzipped liftover chain; this should be '
                  'hg38ToHg19.over.chain')
        )
    parser.add_argument('--intropolis', type=str, required=True,
            help=('path to intropolis.v1.hg19.tsv.gz')
        )
    parser.add_argument('--unmapped', type=str, required=True,
            help='BED in which unmapped junctions should be stored'
        )
    args = parser.parse_args()
    temp_dir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, temp_dir)
    current_dir = os.path.abspath(os.path.dirname(__file__))
    # Read annotated junctions
    annotated_junctions = set()
    with gzip.open(
                os.path.join(current_dir, 'annotated_junctions.tsv.gz')
            ) as annotation_stream:
        for line in annotation_stream:
            tokens = line.strip().split('\t')
            annotated_junctions.add(
                    (tokens[0], str(int(tokens[1]) - 1), tokens[2], tokens[3])
                ) # zero-based, half-open
    # Convert translatome junctions from mm10 to hg19
    temp_mm10 = os.path.join(temp_dir, 'mm10.bed')
    temp_hg19 = os.path.join(temp_dir, 'hg19.bed')
    with open(temp_mm10, 'w') as mm10_stream, gzip.open(
            os.path.join(current_dir, 'mm10_translatome_junctions.tsv.gz')
        ) as input_stream:
        for i, line in enumerate(input_stream):
            tokens = line.strip().split('\t')
            chrom, strand, start, end = (
                    tokens[0][:-1], tokens[0][-1], str(int(tokens[1]) - 1),
                    tokens[2]
                ) # zero-based, half-open coordinates
            # Tack original junction onto junction name
            junction_name = ';'.join([str(i), chrom, start, end, strand,
                                        tokens[3], tokens[4]])
            print >>mm10_stream, '{}\t{}\t{}\tinfo_{}\t1\t{}'.format(
                    chrom, start, end, junction_name, strand
                )
    liftover_process = subprocess.call(' '.join([
                                            args.liftover,
                                            temp_mm10,
                                            args.chain,
                                            temp_hg19,
                                            args.unmapped
                                        ]),
                                        shell=True,
                                        executable='/bin/bash'
                                    )
    to_sort = os.path.join(temp_dir, 'intropolis_and_translatome.tsv.gz')
    with open(temp_hg19) as hg19_stream:
        with gzip.open(to_sort, 'w') as both_stream:
            for line in hg19_stream:
                chrom, start, end, name, score, strand = line.strip().split(
                                                                        '\t'
                                                                    )[:6]
                if (chrom, start, end, strand) not in annotated_junctions:
                    # Only care about unannotated variation
                    (_, mm10_chrom, mm10_start, mm10_end, mm10_strand,
                        mm10_samples, mm10_coverages) = name.split(';')
                    start, mm10_start = int(start), int(mm10_start)
                    if end - start >= 4:
                        # Only write lifted-over introns >= 4 bases long
                        print >>both_stream, '\t'.join(
                                [chrom, str(start + 1),
                                    end, strand, mm10_chrom,
                                    str(mm10_start + 1), mm10_end,
                                    mm10_strand, mm10_samples, mm10_coverages]
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
        last_junction, last_tokens = None, None
        for line in sorted_stream:
            tokens = line.strip().split('\t')
            junction = tokens[:4]
            if junction == last_junction:
                if len(last_tokens) > len(tokens):
                    mm10_tokens = last_tokens
                    hg19_tokens = tokens
                else:
                    mm10_tokens = tokens
                    hg19_tokens = last_tokens
                print '\t'.join()