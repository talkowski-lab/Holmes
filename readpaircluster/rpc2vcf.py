#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import sys
import argparse
import re
from collections import Counter
import numpy as np
import pysam
from rpc import ReadPair

META = """##fileformat=VCFv4.2
##source=RPC
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromosome in interchromosomal events">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">
##FORMAT=<ID=MAPQ,Number=1,Type=Float,Description="Average MapQ of reads in cluster">
##FORMAT=<ID=UNIQ,Number=1,Type=Float,Description="Cluster uniqueness = (# unique mapping positions) / (# total reads)">
##FORMAT=<ID=ALTREF,Number=1,Type=Float,Description="Alt/ref ratio = (2 * cluster size) / ((2 * cluster size) + (local coverage at A) + (local coverage at B))">
##FORMAT=<ID=GCOV,Number=1,Type=Float,Description="Global coverage ratio = (cluster size) / (median library coverage)">
##FORMAT=<ID=MAPQA,Number=1,Type=Float,Description="Average MapQ of reads at A">
##FORMAT=<ID=MAPQB,Number=1,Type=Float,Description="Average MapQ of reads at B">
##FORMAT=<ID=UNIQA,Number=1,Type=Float,Description="Uniqueness of mapping positions at A">
##FORMAT=<ID=UNIQB,Number=1,Type=Float,Description="Uniqueness of mapping positions at B">
##FORMAT=<ID=ALTREFA,Number=1,Type=Float,Description="Alt/ref ratio = (cluster size) / (cluster size + local coverage at A)">
##FORMAT=<ID=ALTREFB,Number=1,Type=Float,Description="Alt/ref ratio = (cluster size) / (cluster size + local coverage at B)">
"""

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n"


def RPCParser(filename=None, fsock=None, prefix=''):
    """
    Parse read pairs from RPC-formatted output.

    Parameters
    ----------
    filename : str
    fsock : handle

    Yields
    ------
    cname : str
        cluster identifier
    cluster : list of ReadPair
    """

    if fsock:
        rpcfile = fsock
    elif filename:
        rpcfile = open(filename)
    else:
        raise Exception('Must specify filename or fsock')

    blank_exp = re.compile(r'^\s*$')
    cluster = []
    cname = ''

    CHROMS = [str(x) for x in range(1, 23)] + 'X Y'.split()

    for line in rpcfile:
        if blank_exp.match(line):
            cname = '_'.join([prefix, cname]).strip('_')

            if cluster[0].chrA in CHROMS and cluster[0].chrB in CHROMS:
                yield RPCluster(cluster, cname)
            cluster = []
        else:
            cluster.append(ReadPair.from_rpc(line))
            cname = line.split()[0]

    cname = '_'.join([prefix, cname]).strip('_')
    if len(cluster) > 0:
        if cluster[0].chrA in CHROMS and cluster[0].chrB in CHROMS:
            yield RPCluster(cluster, cname)


class RPCluster(object):
    def __init__(self, pairs, name):

        self.pairs = pairs

        if len(set(rp.chrA for rp in self.pairs)) > 1:
            raise Exception('Mismatched chromosomes')
        if len(set(rp.chrB for rp in self.pairs)) > 1:
            raise Exception('Mismatched chromosomes')

        self.CHROM = self.pairs[0].chrA
        self.chrB = self.pairs[0].chrB

        self.strands = Counter([pair.strands for pair in self.pairs])
        self.POS, self.END = self._get_pos()
        self.ID = name
        self.REF = 'N'
        self.QUAL = '.'
        self.FILTER = '.'

        self.PE = len(self.pairs)

        def _uniq(starts):
            return len(set(starts)) / len(starts)

        self.uniqA = _uniq([pair.read1.reference_start for pair in self.pairs])
        self.uniqB = _uniq([pair.read2.reference_start for pair in self.pairs])
        self.UNIQ = _uniq([pair.read1.reference_start for pair in self.pairs] +
                          [pair.read2.reference_start for pair in self.pairs])

        self.mapqA = np.mean([pair.read1.mapq for pair in self.pairs])
        self.mapqB = np.mean([pair.read2.mapq for pair in self.pairs])
        self.MAPQ = np.mean([pair.read1.mapq for pair in self.pairs] +
                            [pair.read2.mapq for pair in self.pairs])

    @property
    def ALT(self):
        if self.SVTYPE == 'BND':
            if len(self.strands.keys()) == 1:
                strand = list(self.strands.keys())[0]
                if strand == '++':
                    alt = 'N]{chrB}:{end}]'
                elif strand == '+-':
                    alt = 'N[{chrB}:{end}['
                elif strand == '-+':
                    alt = ']{chrB}:{end}]N'
                elif strand == '--':
                    alt = '[{chrB}:{end}[N'
            else:
                alt = '<{chrB}:{end}>'
            return alt.format(chrB=self.chrB, end=self.END)
        else:
            return '<' + self.SVTYPE + '>'

    def _get_pos(self):
        posA = [rp.read1.reference_start for rp in self.pairs]
        plusA = [rp.read1.reference_start for rp in self.pairs
                 if rp.strands.startswith('+')]
        minusA = [rp.read1.reference_start for rp in self.pairs
                  if rp.strands.startswith('-')]

        posB = [rp.read2.reference_start for rp in self.pairs]
        plusB = [rp.read2.reference_start for rp in self.pairs
                 if rp.strands.startswith('+')]
        minusB = [rp.read2.reference_start for rp in self.pairs
                  if rp.strands.startswith('-')]

        svtype = self.SVTYPE

        if svtype == 'DEL':
            POS, END = max(posA), min(posB)

        elif svtype == 'DUP':
            POS, END = min(posA), max(posB)

        elif svtype == 'INV':
            # Average ++ and -- strand positions
            POS = np.mean([max(plusA), min(minusA)])
            END = np.mean([max(plusB), min(minusB)])

        elif svtype == 'BND':
            if len(self.strands.keys()) == 1:
                strand = list(self.strands.keys())[0]
                if strand == '++':
                    POS, END = max(posA), max(posB)
                elif strand == '+-':
                    POS, END = max(posA), min(posB)
                elif strand == '-+':
                    POS, END = min(posA), max(posB)
                elif strand == '--':
                    POS, END = min(posA), min(posB)

            else:
                As = []
                if plusA:
                    As.append(max(plusA))
                if minusA:
                    As.append(min(minusA))
                POS = np.mean(As)

                Bs = []
                if plusB:
                    Bs.append(max(plusB))
                if minusB:
                    Bs.append(min(minusB))
                END = np.mean(Bs)

        return int(POS), int(END)

    @property
    def SVTYPE(self):
        if self.chrB != self.CHROM:
            return 'BND'

        if len(self.strands.keys()) == 1:
            strand = list(self.strands.keys())[0]
            if strand == '++':
                return 'BND'
            elif strand == '--':
                return 'BND'
            elif strand == '+-':
                return 'DEL'
            elif strand == '-+':
                return 'DUP'
            else:
                raise Exception('Invalid strand')
        else:
            strands = list(self.strands.keys())
            if len(strands) > 2:
                return 'BND'
            elif '++' in strands and '--' in strands:
                return 'INV'
            else:
                return 'BND'

    @property
    def INFO(self):
        strands = sorted(self.strands.items(), key=lambda tup: tup[0])
        strands = ['{0}:{1}'.format(s, c) for s, c in strands]
        strands = ','.join(strands)

        if self.SVTYPE == 'BND':
            INFO = ('SVTYPE={svtype};SVLEN={svlen};CHR2={chr2};END={end};'
                    'STRANDS={strands};IMPRECISE')
            INFO = INFO.format(svtype=self.SVTYPE, svlen=(self.POS-self.END),
                               end=self.END, strands=strands, chr2=self.chrB)
        else:
            INFO = ('SVTYPE={svtype};SVLEN={svlen};END={end};'
                    'STRANDS={strands};IMPRECISE')
            INFO = INFO.format(svtype=self.SVTYPE, svlen=(self.POS-self.END),
                               end=self.END, strands=strands)
        return INFO

    def vcf(self, sample, mean_cov, cov_matrix):
        gcov_ratio = self.PE / mean_cov

        header = next(cov_matrix.header).decode('utf-8')
        cov_idx = header.split().index(sample)

        ref = 0
        localA = next(cov_matrix.fetch(self.CHROM, self.POS, self.POS+1))
        localA = int(localA.split()[cov_idx])
        ref = ref + localA
        localA = self.PE / float(self.PE + localA)

        localB = next(cov_matrix.fetch(self.chrB, self.END, self.END+1))
        localB = int(localB.split()[cov_idx])
        ref = ref + localB
        localB = self.PE / float(self.PE + localB)

        ALTREF = (2 * self.PE) / float(2 * self.PE + ref)

        record = '\t'.join('{CHROM} {POS} {ID} {REF} {ALT} {QUAL} {FILTER} '
                           '{INFO} {FORMAT} {CALL}'.split())

        FORMAT = ('GT:PE:MAPQ:UNIQ:ALTREF:GCOV:'
                  'MAPQA:MAPQB:UNIQA:UNIQB:ALTREFA:ALTREFB')

        CALL = ('{GT}:{PE}:{MAPQ:.01f}:{UNIQ:.03f}:{ALTREF:.03f}:{GCOV:.03f}:'
                '{MAPQA:.01f}:{MAPQB:.01f}:'
                '{UNIQA:.03f}:{UNIQB:.03f}:'
                '{ALTREFA:.03f}:{ALTREFB:.03f}')

        CALL = CALL.format(GT='./.',
                           PE=self.PE,
                           MAPQ=self.MAPQ,
                           UNIQ=self.UNIQ,
                           ALTREF=ALTREF,
                           GCOV=gcov_ratio,
                           MAPQA=self.mapqA,
                           MAPQB=self.mapqB,
                           UNIQA=self.uniqA,
                           UNIQB=self.uniqB,
                           ALTREFA=localA,
                           ALTREFB=localB)

        record = record.format(CHROM=self.CHROM,
                               POS=self.POS,
                               ID=self.ID,
                               REF=self.REF,
                               ALT=self.ALT,
                               QUAL=self.QUAL,
                               FILTER=self.FILTER,
                               INFO=self.INFO,
                               FORMAT=FORMAT,
                               CALL=CALL)
        return record


def main():
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument('-r', '--rpcfiles', nargs=4,
                        help='del dup inv tloc')
    parser.add_argument('vcf', type=argparse.FileType('w'),
                        nargs='?', default=sys.stdout)
    parser.add_argument('-s', '--sample')
    parser.add_argument('-m', '--cov_matrix', type=pysam.TabixFile)
    parser.add_argument('-c', '--gcov', type=int)
    args = parser.parse_args()

    args.vcf.write(META)
    args.vcf.write(HEADER.format(SAMPLE=args.sample))

    svtypes = 'del dup inv tloc'.split()

    for svtype, rpcfile in zip(svtypes, args.rpcfiles):
        prefix = '_'.join([args.sample, svtype])
        for cluster in RPCParser(filename=rpcfile, prefix=prefix):
            record = cluster.vcf(args.sample, args.gcov, args.cov_matrix)
            args.vcf.write(record + '\n')


if __name__ == '__main__':
    main()
