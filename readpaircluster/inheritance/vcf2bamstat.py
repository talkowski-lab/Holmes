#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys
import pysam
from rpc import ReadPair
from array import array

def is_smaller_chrom(chrA, chrB, le=True):
    # Numeric comparison, if possible
    if chrA.isdigit() and chrB.isdigit():
        if le:
            return int(chrA) <= int(chrB)
        else:
            return int(chrA) < int(chrB)

    # String comparison for X/Y
    elif not chrA.isdigit() and not chrB.isdigit():
        if le:
            return chrA <= chrB
        else:
            return chrA < chrB

    # Numeric is always less than X/Y
    else:
        return chrA.isdigit()


class callRP(ReadPair):
    @classmethod
    def from_vcf_record(cls, record, sample, prog, read_len=151):
        read1 = pysam.AlignedRead()
        read2 = pysam.AlignedRead()

        svtype = record.INFO['SVTYPE']
        if svtype not in ['TRA', 'BND']:
            chrB = chrA
            posB = record.INFO['END']
        else:
            # Delly format
            if svtype == 'TRA':
                chrB = record.INFO['CHR2']
                posB = record.INFO['END']
            # Lumpy format
            elif svtype == 'BND':
                if record.ID.endswith('_2'):
                    return None
                chrB = record.ALT[0].chr
                posB = record.ALT[0].pos

        read1.reference_start = posA
        read2.reference_start = posB

        read1.query_name = '{0}__{1}__{2}'.format(sample, prog, ID)
        read2.query_name = '{0}__{1}__{2}'.format(sample, prog, ID)

        read1.flag = 1
        read2.flag = 1

        read1.mapping_quality = 37
        read2.mapping_quality = 37

        read1.query_sequence = ('N' * read_len).encode('utf-8')
        read1.query_qualities = array('B', [0 for x in range(read_len)])
        read2.query_sequence = ('N' * read_len).encode('utf-8')
        read2.query_qualities = array('B', [0 for x in range(read_len)])

        return ReadPair(read1, read2, chrA, chrB)
        


    @classmethod
    def from_vcf(cls, line, sample, prog, read_len=151):
        """
        Construct a ReadPair object from a VCF record

        Parameters
        ----------
        line : str
            A line from a VCF file

        Returns
        -------
        ReadPair
        """
        if line.startswith('#'):
            raise Exception('Header line not allowed')

        read1 = pysam.AlignedRead()
        read2 = pysam.AlignedRead()

        data = line.strip().split()
        chrA = data[0]
        posA = int(data[1])
        ID = data[2]
        ref = data[3]
        alt = data[4]
        qual = data[5]
        filt = data[6]
        INFO = {k: v for k, v in [f.split('=') for f in data[7].split(';') if '=' in f]}
        FORMAT = data[8]
        gt = data[9:]

        # chrA = chrom
        # posA = pos
        svtype = INFO['SVTYPE'] 
        if svtype not in ['TRA', 'BND']:
            chrB = chrA
            posB = INFO['END']
        else:
            # Delly format
            if svtype == 'TRA':
                chrB = INFO['CHR2']
                posB = INFO['END']
            # Lumpy format
            elif svtype == 'BND':
                if ID.endswith('_2'):
                    return None
                chrB, posB = alt.strip('[]N').split(':')
        posB = int(posB)

        if chrA == chrB:
            posA, posB = (min(posA, posB), max(posA, posB))
        elif not is_smaller_chrom(chrA, chrB):
            chrA, chrB = (chrB, chrA)
            posA, posB = (posB, posA)

        read1.reference_start = posA
        read2.reference_start = posB

        read1.query_name = '{0}__{1}__{2}'.format(sample, prog, ID)
        read2.query_name = '{0}__{1}__{2}'.format(sample, prog, ID)

        read1.flag = 1
        read2.flag = 1

        read1.mapping_quality = 37
        read2.mapping_quality = 37

        read1.query_sequence = ('N' * read_len).encode('utf-8')
        read1.query_qualities = array('B', [0 for x in range(read_len)])
        read2.query_sequence = ('N' * read_len).encode('utf-8')
        read2.query_qualities = array('B', [0 for x in range(read_len)])

        return ReadPair(read1, read2, chrA, chrB)


def VCFParser(fin, sample, prog):
    for line in fin:
        if line.startswith('#'):
            continue
        if line.startswith('GL') or line.startswith('MT'):
            continue
        yield callRP.from_vcf(line, sample, prog)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument('vcf', type=argparse.FileType('r'))
    parser.add_argument('sample')
    parser.add_argument('prog')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        nargs='?', default=sys.stdout)

    args = parser.parse_args()

    for rp in VCFParser(args.vcf, args.sample, args.prog):
        if rp is not None:
            args.fout.write(rp.as_bamstat() + '\n')


if __name__ == '__main__':
    main()
