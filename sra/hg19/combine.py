#!/usr/bin/env python
"""
combine.py
Abhi Nellore / June 28, 2015

Reformats sorted and merged "collected intron" file from SRA batch runs
into a file with one line per intron. Batch index + original sample indexes
are remapped as follows.

new sample index = 500 * batch index + old sample index.

The original introns are in a format where the start position is 1-based and
inclusive and the end position is 1-based and exclusive. The format of each
intron is changed so the end position is 1-based and inclusive. An extra
field is added to indicate whether the intron is GT-AG, GC-AG, or AT-AC, so
a Bowtie index basename must be added at the command line.

Also writes a file mapping new indexes to SRA accession numbers.
"""
import sys
import itertools
import glob

# This code is taken from bowtie_index.py in Rail-RNA
import os
import struct
import mmap
from operator import itemgetter
from collections import defaultdict
from bisect import bisect_right

class BowtieIndexReference(object):
    """
    Given prefix of a Bowtie index, parses the reference names, parses the
    extents of the unambiguous stretches, and memory-maps the file containing
    the unambiguous-stretch sequences.  get_stretch member function can
    retrieve stretches of characters from the reference, even if the stretch
    contains ambiguous characters.
    """

    def __init__(self, idx_prefix):

        # Open file handles
        if os.path.exists(idx_prefix + '.3.ebwt'):
            # Small index (32-bit offsets)
            fh1 = open(idx_prefix + '.1.ebwt', 'rb')  # for ref names
            fh3 = open(idx_prefix + '.3.ebwt', 'rb')  # for stretch extents
            fh4 = open(idx_prefix + '.4.ebwt', 'rb')  # for unambiguous sequence
            sz, struct_unsigned = 4, struct.Struct('I')
        else:
            raise RuntimeError('No Bowtie index files with prefix "%s"' % idx_prefix)

        #
        # Parse .1.bt2 file
        #
        one = struct.unpack('<i', fh1.read(4))[0]
        assert one == 1

        ln = struct_unsigned.unpack(fh1.read(sz))[0]
        line_rate = struct.unpack('<i', fh1.read(4))[0]
        lines_per_side = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]
        ftab_chars = struct.unpack('<i', fh1.read(4))[0]
        _ = struct.unpack('<i', fh1.read(4))[0]

        nref = struct_unsigned.unpack(fh1.read(sz))[0]
        # get ref lengths
        reference_length_list = []
        for i in xrange(nref):
            reference_length_list.append(struct.unpack('<i', fh1.read(sz))[0])

        nfrag = struct_unsigned.unpack(fh1.read(sz))[0]
        # skip rstarts
        fh1.seek(nfrag * sz * 3, 1)

        # skip ebwt
        bwt_sz = ln // 4 + 1
        line_sz = 1 << line_rate
        side_sz = line_sz * lines_per_side
        side_bwt_sz = side_sz - 8
        num_side_pairs = (bwt_sz + (2*side_bwt_sz) - 1) // (2*side_bwt_sz)
        ebwt_tot_len = num_side_pairs * 2 * side_sz
        fh1.seek(ebwt_tot_len, 1)

        # skip zOff
        fh1.seek(sz, 1)

        # skip fchr
        fh1.seek(5 * sz, 1)

        # skip ftab
        ftab_len = (1 << (ftab_chars * 2)) + 1
        fh1.seek(ftab_len * sz, 1)

        # skip eftab
        eftab_len = ftab_chars * 2
        fh1.seek(eftab_len * sz, 1)

        refnames = []
        while True:
            refname = fh1.readline()
            if len(refname) == 0 or ord(refname[0]) == 0:
                break
            refnames.append(refname.split()[0])
        assert len(refnames) == nref

        #
        # Parse .3.bt2 file
        #
        one = struct.unpack('<i', fh3.read(4))[0]
        assert one == 1

        nrecs = struct_unsigned.unpack(fh3.read(sz))[0]

        running_unambig, running_length = 0, 0
        self.recs = defaultdict(list)
        self.offset_in_ref = defaultdict(list)
        self.unambig_preceding = defaultdict(list)
        length = {}

        ref_id, ref_namenrecs_added = 0, None
        for i in xrange(nrecs):
            off = struct_unsigned.unpack(fh3.read(sz))[0]
            ln = struct_unsigned.unpack(fh3.read(sz))[0]
            first_of_chromosome = ord(fh3.read(1)) != 0
            if first_of_chromosome:
                if i > 0:
                    length[ref_name] = running_length
                ref_name = refnames[ref_id]
                ref_id += 1
                running_length = 0
            assert ref_name is not None
            self.recs[ref_name].append((off, ln, first_of_chromosome))
            self.offset_in_ref[ref_name].append(running_length)
            self.unambig_preceding[ref_name].append(running_unambig)
            running_length += (off + ln)
            running_unambig += ln

        length[ref_name] = running_length
        assert nrecs == sum(map(len, self.recs.itervalues()))

        #
        # Memory-map the .4.bt2 file
        #
        ln_bytes = (running_unambig + 3) // 4
        self.fh4mm = mmap.mmap(fh4.fileno(), ln_bytes, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

        # These are per-reference
        self.length = length
        self.refnames = refnames

        # To facilitate sorting reference names in order of descending length
        sorted_rnames = sorted(self.length.items(),
                               key=lambda x: itemgetter(1)(x), reverse=True)
        self.rname_to_string = {}
        self.string_to_rname = {}
        for i, (rname, _) in enumerate(sorted_rnames):
            rname_string = ('%012d' % i)
            self.rname_to_string[rname] = rname_string
            self.string_to_rname[rname_string] = rname
        # Handle unmapped reads
        unmapped_string = ('%012d' % len(sorted_rnames))
        self.rname_to_string['*'] = unmapped_string
        self.string_to_rname[unmapped_string] = '*'

        # For compatibility
        self.rname_lengths = self.length

    def get_stretch(self, ref_id, ref_off, count):
        """
        Return a stretch of characters from the reference, retrieved
        from the Bowtie index.

        @param ref_id: name of ref seq, up to & excluding whitespace
        @param ref_off: offset into reference, 0-based
        @param count: # of characters
        @return: string extracted from reference
        """
        assert ref_id in self.recs
        # Account for negative reference offsets by padding with Ns
        N_count = min(abs(min(ref_off, 0)), count)
        stretch = ['N'] * N_count
        count -= N_count
        if not count: return ''.join(stretch)
        ref_off = max(ref_off, 0)
        starting_rec = bisect_right(self.offset_in_ref[ref_id], ref_off) - 1
        assert starting_rec >= 0
        off = self.offset_in_ref[ref_id][starting_rec]
        buf_off = self.unambig_preceding[ref_id][starting_rec]
        # Naive to scan these records linearly; obvious speedup is binary search
        for rec in self.recs[ref_id][starting_rec:]:
            off += rec[0]
            while ref_off < off and count > 0:
                stretch.append('N')
                count -= 1
                ref_off += 1
            if count == 0:
                break
            if ref_off < off + rec[1]:
                # stretch extends through part of the unambiguous stretch
                buf_off += (ref_off - off)
            else:
                buf_off += rec[1]
            off += rec[1]
            while ref_off < off and count > 0:
                buf_elt = buf_off >> 2
                shift_amt = (buf_off & 3) << 1
                stretch.append(
                    'ACGT'[(ord(self.fh4mm[buf_elt]) >> shift_amt) & 3]
                )
                buf_off += 1
                count -= 1
                ref_off += 1
            if count == 0:
                break
        # If the requested stretch went past the last unambiguous
        # character in the chromosome, pad with Ns
        while count > 0:
            count -= 1
            stretch.append('N')
        return ''.join(stretch)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bowtie-idx', type=str, required=True,
        help='Path to Bowtie index basename')
    parser.add_argument('--fix-batch-9', action='store_const',
        const=True, default=False,
        help='Uses old manifest file for sample indexes from '
             'batch 9. Addresses how during our batch runs on '
             'all of SRA, we used a manifest file with a '
             'sample eliminated to preprocess the data but '
             'a different manifest file with that sample '
             'included to search for introns. See NOTES '
             'for more information.')
    args = parser.parse_args()
    # Write index-to-accession file
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    reversed_complements = {
            ('CT', 'AC') : ('GT', 'AG'),
            ('CT', 'GC') : ('GC', 'AG'),
            ('GT', 'AT') : ('AT', 'AC')
        }
    filenames = [glob.glob(
                            os.path.join(
                                containing_dir,
                                'sra_batch_%d_sample_size*.txt' % i
                            )
                        )[0] for i in xrange(43)]
    if args.fix_batch_9:
        # Use old manifest file with 500 samples
        filenames[9] = os.path.join(containing_dir,
                            'sra_batch_9_sample_size_500_old.txt')
    else:
        # Use new manifest file with 499 samples
        filenames[9] = os.path.join(containing_dir,
                            'sra_batch_9_sample_size_500.txt')
    with open('index_to_SRA_accession.tsv', 'w') as output_stream:
        for i, filename in enumerate(filenames):
            with open(filename) as input_stream:
                for j, line in enumerate(input_stream):
                    tokens = line.strip().split('\t')
                    if args.fix_batch_9 \
                        and 'SRP000941_SRS306616_SRX190128_SRR651690-1-1' \
                        in tokens:
                        # ignore sample because it wasn't found
                        continue
                    print >>output_stream, (str(i * 500 + j) + '\t'
                            + '\t'.join(tokens[-1][:-4].split('_')))

    reference_index = BowtieIndexReference(args.bowtie_idx)
    for intron, lines in itertools.groupby(
                            sys.stdin, key=lambda x: x.split('\t')[1:4]
                        ):
        chrom = intron[0][:-1]
        start = int(intron[1])
        end = int(intron[2]) - 1
        start_motif = reference_index.get_stretch(chrom, start - 1, 2)
        end_motif = reference_index.get_stretch(chrom, end - 2, 2)
        if intron[0][-1] == '-':
            start_motif, end_motif = reversed_complements[
                                            (start_motif, end_motif)
                                        ]
        pairs = []
        for line in lines:
            tokens = line.strip().split('\t')
            batch_index = int(tokens[0])
            coverages = tokens[-1].split(',')
            sample_indexes = [500 * batch_index + int(original_index)
                                for original_index in tokens[-2].split(',')]
            pairs.extend(zip(sample_indexes, coverages))
        pairs.sort(key=lambda x: x[0])
        print '\t'.join([chrom, str(start), str(end), intron[0][-1],
                            start_motif, end_motif,
                            ','.join([str(pair[0]) for pair in pairs]),
                            ','.join([str(pair[1]) for pair in pairs])])