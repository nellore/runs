#!/usr/bin/env python
"""
combine_sra.py
Abhi Nellore / March 10, 2015

Reformats sorted and merged "collected junctions" files from SRA batch runs
into files with one line per junction. Three files are written:

intropolis.v2.hg38.tsv.gz / intropolis.2pass.v2.hg38.tsv.gz
Coverages reported are after first- or second-pass alignment
1. chromosome
2. start position
3. end position
4. strand (+ or -)
5. start motif (e.g., GT)
6. end motif (e.g., AG)
7. comma-separated list of sample indexes
8. comma-separated list of coverages

intropolis.allpasses.v2.hg38.tsv.gz
Puts coverages after first- and second-pass alignment in same file
1. chromosome
2. start position
3. end position
4. strand (+ or -)
5. start motif (e.g., GT)
6. end motif (e.g., AG)
7. comma-separated list of sample indexes
8. comma-separated list of coverages after first-pass alignment (some can be
                                                                    zero)
9. comma-separated list of coverages after second-pass alignment (some can be
                                                                    zero)
"""
import sys
import itertools
import glob
import gzip
import tempfile
import atexit
import subprocess
import shutil

# This code is taken from bowtie_index.py in Rail-RNA
import os
import struct
import mmap
from operator import itemgetter
from collections import defaultdict
from bisect import bisect_right
import csv

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
        help='path to Bowtie index basename')
    parser.add_argument('--sra-dir', required=True,
        help='path to SRA output files; this is where the batch_* '
             'subdirectories are')
    parser.add_argument('--output-dir', required=True,
        help='directory in which to write output files')
    parser.add_argument('--temp-dir', required=False,
        default=None,
        help='where temporary files should be stored')
    parser.add_argument('--sort', required=False,
        default='sort',
        help='path to sort executable')
    args = parser.parse_args()
    containing_dir = os.path.dirname(os.path.realpath(__file__))
    # Create original sample index to new sample index map
    original_index_to_final_index = {}
    final_index_to_sample_name = {}
    sample_name_to_final_index = {}
    batch_numbers = []
    i = 0
    for manifest in glob.glob(os.path.join(containing_dir, '*.manifest')):
        batch_number = int(manifest.partition('.')[0].rpartition('_')[2])
        batch_numbers.append(batch_number)
        with open(manifest) as manifest_stream:
            j = 0
            for line in manifest_stream:
                line = line.strip()
                if line[0] == '#' or not line: continue
                sample_name = line.partition('\t')[0].partition(':')[2]
                original_index_to_final_index[(batch_number, j)] = i
                final_index_to_sample_name[i] = sample_name
                sample_name_to_final_index[sample_name] = i
                i += 1
                j += 1
    sample_name_to_line = {}
    with open(os.path.join(containing_dir, 'SraRunInfo.csv')) as run_stream:
        run_stream.readline()
        run_reader = csv.reader(run_stream, delimiter=',', quotechar='"')
        for tokens in run_reader:
            if not tokens or (len(tokens) == 1 and tokens[0] == ''): continue
            sample_name_to_line[tokens[0]] = '\t'.join(
                    [tokens[20], tokens[24], tokens[10], tokens[0]]
                )
    with open(
            os.path.join(args.output_dir, 'intropolis.idmap.v2.hg38.tsv'), 'w'
        ) as sample_stream:
        for i in sorted(final_index_to_sample_name.keys()):
            print >>sample_stream, '{}\t{}'.format(
                    i, sample_name_to_line[final_index_to_sample_name[i]]
                )
    # Junction are already sorted, so we can advance counters and print
    first_pass_handles = [
        gzip.open(
            os.path.join(args.sra_dir, 'batch_%d' % batch_number,
                         'cross_sample_results', 'first_pass_junctions.tsv.gz')
        ) for batch_number in batch_numbers]
    second_pass_handles = [
        gzip.open(
            os.path.join(args.sra_dir, 'batch_%d' % batch_number,
                         'cross_sample_results', 'junctions.tsv.gz')
        ) for batch_number in batch_numbers]
    column_to_final_index = {}
    for i, batch_number in enumerate(batch_numbers):
        for j, name in enumerate(
                        second_pass_handles[i].readline().strip().split('\t')
                    ):
            column_to_final_index[
                    (batch_number, j)
                ] = sample_name_to_final_index[name.partition('_')[0]]
    if args.temp_dir is not None:
        temp_dir = tempfile.mkdtemp(dir=args.temp_dir)
    else:
        temp_dir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, temp_dir)
    terminate = False
    temp_file = os.path.join(temp_dir, 'temp.tsv')
    import shutil
    with open(temp_file, 'w') as temp_stream:
        while not terminate:
            terminate = True
            for i, batch_number in enumerate(batch_numbers):
                tokens = first_pass_handles[i].readline().strip().split('\t')
                if tokens[0] != '':
                    terminate = False
                    sample_indexes = [original_index_to_final_index[
                                            (batch_number, int(original_index))
                                        ] for original_index
                                        in tokens[3].split(',')]
                    junction = (tokens[0][:-1], int(tokens[1]),
                                    int(tokens[2]), tokens[0][-1])
                    print >>temp_stream, '{}\t{}\t{}\t{}\t{}\t{}\t0'.format(
                            *(junction + (','.join([str(sample_index)
                                                    for sample_index
                                                    in sample_indexes]),
                                          tokens[4]))
                        )
                tokens = second_pass_handles[i].readline().strip().split('\t')
                if tokens[0] != '':
                    terminate = False
                    junction = tokens[0].split(';')
                    junction = (junction[0], int(junction[2]),
                                     int(junction[3]), junction[1])
                    sample_indexes = [column_to_final_index[
                                                    (batch_number, column)
                                                ] for column, coverage
                                        in enumerate(tokens[1:])
                                        if coverage != '0']
                    coverages = [coverage for coverage in tokens[1:]
                                        if coverage != '0']
                    print >>temp_stream, \
                        '{}\t{}\t{}\t{}\t{}\t{}\t1'.format(
                            *(junction + (','.join([str(sample_index)
                                                    for sample_index
                                                    in sample_indexes]),
                                          ','.join(coverages)))
                        )
    sort_process = subprocess.check_call(args.sort + ' -k1,1 -k2,2n -k3,3n '
                                            + (('-T ' + temp_dir + ' ')
                                                if args.temp_dir else '')
                                            + temp_file + ' >'
                                            + temp_file + '.sorted', 
                                            shell=True,
                                            executable='/bin/bash')
    import itertools
    reference_index = BowtieIndexReference(args.bowtie_idx)
    reversed_complements = {
            ('CT', 'AC') : ('GT', 'AG'),
            ('CT', 'GC') : ('GC', 'AG'),
            ('GT', 'AT') : ('AT', 'AC')
        }
    try:
        os.makedirs(args.output_dir)
    except OSError as e:
        if 'File exists' not in e: 
            raise
    with open(temp_file + '.sorted') as temp_stream, gzip.open(
            os.path.join(args.output_dir, 'intropolis.v2.hg38.tsv.gz'), 'w'
        ) as first_pass_stream, gzip.open(
            os.path.join(args.output_dir,
                            'intropolis.2pass.v2.hg38.tsv.gz'), 'w'
        ) as second_pass_stream, gzip.open(
            os.path.join(args.output_dir,
                            'intropolis.allpasses.v2.hg38.tsv.gz'), 'w'
        ) as consolidated_stream:
        for key, group in itertools.groupby(temp_stream,
                                            key=lambda x: x.split('\t')[:4]):
            first_pass_coverages, second_pass_coverages = [], []
            first_pass_dict, second_pass_dict = defaultdict(int), \
                defaultdict(int)
            for line in group:
                tokens = line.strip().split('\t')
                sample_indexes = [int(el) for el in tokens[-3].split(',')]
                coverages = [int(el) for el in tokens[-2].split(',')]
                together = zip(sample_indexes, coverages)
                if tokens[-1] == '0':
                    first_pass_coverages.extend(together)
                    for sample_index, coverage in together:
                        first_pass_dict[sample_index] = coverage
                else:
                    assert tokens[-1] == '1'
                    second_pass_coverages.extend(together)
                    for sample_index, coverage in together:
                        second_pass_dict[sample_index] = coverage
            start_motif = reference_index.get_stretch(key[0], int(key[1]) - 1,
                                                        2)
            end_motif = reference_index.get_stretch(key[0], int(key[2]) - 2, 2)
            if key[3] == '-':
                start_motif, end_motif = reversed_complements[
                                               (start_motif, end_motif)
                                            ]
            first_pass_coverages.sort()
            second_pass_coverages.sort()
            all_sample_indexes = sorted(
                    list(set(first_pass_dict.keys() + second_pass_dict.keys()))
                )
            consolidated_first_pass_coverages = ','.join(
                    [str(first_pass_dict[sample_index]) for sample_index
                            in all_sample_indexes]
                )
            consolidated_second_pass_coverages = ','.join(
                    [str(second_pass_dict[sample_index]) for sample_index
                            in all_sample_indexes]
                )
            assert first_pass_coverages
            print >>first_pass_stream, '\t'.join(
                key + [start_motif, end_motif,
                ','.join([str(el[0]) for el in first_pass_coverages]),
                ','.join([str(el[1]) for el in first_pass_coverages])])
            if second_pass_coverages:
                print >>second_pass_stream, '\t'.join(
                    key + [start_motif, end_motif,
                    ','.join([str(el[0]) for el in second_pass_coverages]),
                    ','.join([str(el[1]) for el in second_pass_coverages])])
            print >>consolidated_stream, '\t'.join(
                    key + [start_motif, end_motif,
                    ','.join([str(el) for el in all_sample_indexes]),
                    consolidated_first_pass_coverages,
                    consolidated_second_pass_coverages]
                )
