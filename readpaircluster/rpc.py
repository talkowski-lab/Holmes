#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""
Perform single linkage clustering on a set of aligned read pairs.
"""

import sys
import argparse
from io import open, IOBase
from collections import deque
from itertools import combinations
from array import array
from operator import itemgetter
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
import pysam
from Bio.SeqIO.QualityIO import _phred_to_sanger_quality_str as phred_ascii
phred_score = dict((char, score) for score, char in phred_ascii.items())


class ReadPair(object):
    """
    Simplified representation of a read pair.

    Parameters
    ----------
    read1 : pysam.AlignedSegment
        First read in pair.
    read2 : pysam.AlignedSegment
        Second read in pair.
    chrA : str
        Reference contig which read1 mapped to.
        Equivalent to pysam.AlignmentFile.getrname(read1.reference_id).
    chrB : str
        Reference contig which read2 mapped to.
        Equivalent to pysam.AlignmentFile.getrname(read2.reference_id).

    Attributes
    ----------
    read1 : pysam.AlignedSegment
        First read in pair.
    read2 : pysam.AlignedSegment
        Second read in pair.
    chrA : str
        Reference contig which read1 mapped to.
        Equivalent to pysam.AlignmentFile.getrname(read1.reference_id).
    chrB : str
        Reference contig which read2 mapped to.
        Equivalent to pysam.AlignmentFile.getrname(read2.reference_id).
    """
    def __init__(self, read1, read2, chrA, chrB):
        self.read1 = read1
        self.read2 = read2
        self.chrA = chrA
        self.chrB = chrB

    @classmethod
    def from_bamstat(cls, line):
        """
        Construct a ReadPair object from a bamstat-formatted string.

        Parameters
        ----------
        line : str
            A line from a file in bamstat format.

        Returns
        -------
        ReadPair
        """
        read1 = pysam.AlignedRead()
        read2 = pysam.AlignedRead()

        data = line.strip().split()
        read1.query_name = data[0].encode('utf-8')
        read2.query_name = data[0].encode('utf-8')
        chrA = data[1]
        read1.reference_start = int(data[2])
        read1.flag = int(data[3])
        chrB = data[4]
        read2.reference_start = int(data[5])
        read2.flag = int(data[6])
        read1.mapping_quality = int(data[8])
        read2.mapping_quality = int(data[10])
        read1.query_sequence = data[11].encode('utf-8')
        read1.query_qualities = array('B', [phred_score[q] for q in data[12]])
        read2.query_sequence = data[13].encode('utf-8')
        read2.query_qualities = array('B', [phred_score[q] for q in data[14]])

        return ReadPair(read1, read2, chrA, chrB)

    @classmethod
    def from_rpc(cls, line):
        """
        Construct a ReadPair object from a rpc-formatted string.

        Parameters
        ----------
        line : str
            A line from a file in bamstat format.

        Returns
        -------
        ReadPair
        """
        read1 = pysam.AlignedRead()
        read2 = pysam.AlignedRead()

        data = line.strip().split()
        read1.query_name = data[2].encode('utf-8')
        read2.query_name = data[2].encode('utf-8')
        chrA = data[3]
        read1.reference_start = int(data[5])
        read1.flag = int(data[11])
        chrB = data[4]
        read2.reference_start = int(data[6])
        read2.flag = int(data[12])
        read1.mapping_quality = int(data[7])
        read2.mapping_quality = int(data[8])
        read1.query_sequence = data[13].encode('utf-8')
        read1.query_qualities = array('B', [phred_score[q] for q in data[15]])
        read2.query_sequence = data[14].encode('utf-8')
        read2.query_qualities = array('B', [phred_score[q] for q in data[16]])

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

        fields = data[7].split(';')
        INFO = {k: v for k, v in [f.split('=') for f in fields if '=' in f]}

        FORMAT = data[8]
        gt = data[9:]

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

    def __eq__(self, other):
        return (self.chrA == other.chrA and
                self.chrB == other.chrB and
                self.read1.reference_start == other.read1.reference_start and
                self.read2.reference_start == other.read2.reference_start and
                self.read1.query_name == other.read1.query_name)

    def __hash__(self):
        return id(self)

    def is_clusterable_with(self, other, cutoff):
        """
        Test if first reads in two ReadPairs are within clustering distance.

        A batch of read pairs is read in from the input file while the first
        read in the next pair within clustering distance of the previous.

        Parameters
        ----------
        other : ReadPair
            ReadPair to test against.
        cutoff : int
            Maximum clustering distance.

        Returns
        -------
        bool
            True if the distance between the first reads in each pair is less
            than the cutoff.
        """
        self_pos1 = self.read1.reference_start
        other_pos1 = other.read1.reference_start

        return (self.chrA == other.chrA and
                self.chrB == other.chrB and
                abs(self_pos1 - other_pos1) < cutoff)

    def clusters_with(self, other, cutoff, min_mapq):
        """
        Test if two ReadPairs pass the requirements to cluster together.

        Clustering requirements:
            1) Both pairs must map to the same pair of reference contigs.
               i.e., self.read1_contig == other.read1_contig and
                     self.read2_contig == other.read2_contig
            2) The pairs of first reads and second reads must fall within the
               clustering cutoff distance.
               i.e., abs(self.read1_pos - other.read1_pos) < cutoff and
                     abs(self.read2_pos - other.read2_pos) < cutoff
            3) The two pairs cannot be duplicates. (are we sure?)
            3) All four reads must have at least the required minimum
               mapping quality.

        Parameters
        ----------
        other : ReadPair
            ReadPair to test against.
        cutoff : int
            Maximum distance at which two pairs will be clustered. Each pair
            of equivalent reads must be within this distance.
        min_mapq : int
            Minimum mapping quality for a pair to be clustered. Both reads in
            both pairs must pass this threshold for the pairs to be clustered
            together.

        Returns
        -------
        bool
            True if the two pairs pass the clustering requirements.
        """

        # TODO: remove chrB, read1 comparisons as unnecessary
        return (self != other and
                self.chrB == other.chrB and
                abs(self.read1.reference_start - other.read1.reference_start) < cutoff and
                abs(self.read2.reference_start - other.read2.reference_start) < cutoff)
                # (self.read1.reference_start != other.read1.reference_start or
                 # self.read2.reference_start != other.read2.reference_start) and
                # TODO: filter mapq out in get_candidates

    def is_in(self, tabixfile):
        """
        Test if either read in the pair overlaps a region in the tabixfile.

        Parameters
        ----------
        tabixfile : pysam.TabixFile
            List of regions to check for overlap.

        Returns
        -------
        bool
            True if either read in the pair overlaps a region in the tabixfile.
        """
        posA = self.read1.reference_start
        lenA = self.read1.query_length
        posB = self.read2.reference_start
        lenB = self.read2.query_length

        return ((self.chrA.encode('utf-8') in tabixfile.contigs and
                 any(tabixfile.fetch(str(self.chrA), posA, posA + lenA - 1))) or
                (self.chrB.encode('utf-8') in tabixfile.contigs and
                 any(tabixfile.fetch(str(self.chrB), posB, posB + lenB - 1))))

    def as_rpc(self):
        """Return a string representation in rpc format."""

        qualA = ''.join([phred_ascii[q] for q in self.read1.query_qualities])
        qualB = ''.join([phred_ascii[q] for q in self.read2.query_qualities])
        
        return '\t'.join([self.read1.query_name,
                          self.chrA,
                          self.chrB,
                          str(self.read1.reference_start),
                          str(self.read2.reference_start),
                          str(self.read1.mapping_quality),
                          str(self.read2.mapping_quality),
                          str(self.read1.query_length),
                          str(self.read2.query_length),
                          str(self.read1.flag),
                          str(self.read2.flag),
                          str(self.read1.query_sequence),
                          str(self.read2.query_sequence),
                          qualA,
                          qualB])

    def as_bamstat(self):
        """Return a string representation in rpc format."""

        qualA = ''.join([phred_ascii[q] for q in self.read1.query_qualities])
        qualB = ''.join([phred_ascii[q] for q in self.read2.query_qualities])

        return '\t'.join([self.read1.query_name,
                          self.chrA,
                          str(self.read1.reference_start),
                          str(self.read1.flag),
                          self.chrB,
                          str(self.read2.reference_start),
                          str(self.read2.flag),
                          str(self.read1.query_length),
                          str(self.read1.mapping_quality),
                          str(self.read2.query_length),
                          str(self.read2.mapping_quality),
                          str(self.read1.query_sequence),
                          qualA,
                          str(self.read2.query_sequence),
                          qualB])

    def __str__(self):
        return '\t'.join([self.read1.query_name,
                          self.chrA,
                          self.chrB,
                          str(self.read1.reference_start),
                          str(self.read2.reference_start)])


def ReadPairParser(handle):
    """
    Iterate over bamstat records and yield ReadPair objects.

    Parameters
    ----------
    handle : file
        File of read pairs in bamstat format.

    Yields
    ------
    ReadPair
    """
    if isinstance(handle, IOBase):
        for line in handle:
            yield ReadPair.from_bamstat(line)
    elif isinstance(handle, pysam.AlignmentFile):
        pass
        # for readA, readB in SamPairParser(handle):
            # chrA = handle.getrname(readA.tid)
            # chrB = handle.getrname(readB.tid)
            # yield ReadPair(readA, readB, chrA, chrB)
    else:
        # sys.stderr.write('Unsupported file format: %s\n' % filetype)
        sys.stderr.write('Unsupported file format: \n')
        sys.exit(1)


# Should mapq filtering happen here? May split clusters
def get_candidate_pairs(pairs, cutoff, excluded_regions=None,
                        excluded_file=None, min_mapq=-1):
    """
    Find batches of ReadPairs eligible for clustering.

    Requires input sorted by chromosome and position of first read in each
    pair. Pairs are collected while the first read in the next pair in the
    parser is within the maximum clustering distance of the first read of the
    previous pair.

    Parameters
    ----------
    pairs : ReadPairParser
        Iterator over sorted ReadPairs.
    cutoff : int
        Maximum clustering distance.
    excluded_regions : pysam.TabixFile or NoneType
        Regions to exclude from clustering. Any read pair that overlaps with
        a region is omitted.

    Yields
    ------
    deque of ReadPair
    """
    candidates = deque()
    prev = None

    excluded = 0
    failed_mapq = 0
    failed_eq = 0
    permitted = 0
    total = 0
    for pair in pairs:
        total += 1
        if prev is not None and pair == prev:
            failed_eq += 1
            continue
        if excluded_regions is not None and pair.is_in(excluded_regions):
            excluded += 1
            if excluded_file is not None:
                excluded_file.write(pair.as_bamstat() + '\n')
            continue

        if pair.read1.mapping_quality <= min_mapq or pair.read2.mapping_quality <= min_mapq:
            failed_mapq += 1
            continue

        if prev is None or prev.is_clusterable_with(pair, cutoff):
            permitted += 1
            candidates.append(pair)
        else:
            permitted += 1
            yield candidates
            candidates = deque([pair])

        prev = pair

    yield candidates

    def log_count(count, msg):
        if total != 0:
            pct = (float(count) / total) * 100
        else:
            pct = 0.0

        msg = msg + ' %d read pairs out of %d total (%.1f%%)\n'
        #  sys.stderr.write(msg % (count, total, pct))

    log_count(excluded, "Blacklisted")
    log_count(failed_mapq, "Mapq failed")
    log_count(permitted, "Clustered")

def cluster_candidates(candidates, cutoff, min_mapq, min_size):
    """
    Perform single linkage clustering on a batch of ReadPairs.

    Parameters
    ----------
    candidates : set of ReadPair
        ReadPairs to cluster.
    cutoff : int
        Maximum clustering distance. Both sets of corresponding reads in two
        pairs must fall within this distance in order for the pairs to be
        clustered together.
    min_mapq : int
        Minimum mapping quality for a read pair to be included in clustering.
        Both reads in pair must exceed this threshold.
    min_size : int
        Minimum cluster size that will be reported.

    Returns
    -------
    list of set of ReadPair
        List of clustered ReadPairs.
    """

    n = len(candidates)
    # Permit clusters of size 1
    G = sparse.eye(n, dtype=np.uint16, format='lil')

    for p1, p2 in combinations(range(n), 2):
        if candidates[p1].clusters_with(candidates[p2], cutoff, min_mapq):
            G[p1, p2] = 1

    n_comp, comp_list = csgraph.connected_components(G, connection='weak')
    cluster_names = np.arange(n_comp)
    cluster_sizes = np.array([np.sum(comp_list == i) for i in cluster_names])

    # Remove clusters with less than minimum size
    cluster_names = cluster_names[np.where(cluster_sizes >= min_size)]

    # TODO: check speed vs creating new sub-deques
    #  clusters = [deque()] * len(cluster_names)
    clusters = deque()
    for cname in cluster_names:
        cluster = itemgetter(*np.where(comp_list == cname)[0])(candidates)
        # Cluster size 1
        if not isinstance(cluster, tuple):
            clusters.append([cluster])
        else:
            clusters.append(list(cluster))

    return clusters


def rpc_format(cluster, cluster_id):
    """
    Converts set of read pairs to RPC entry.

    Arguments
    ---------
    cluster : set of ReadPair
        ReadPairs to format.
    cluster_id : str
        Cluster identifier. Usually the index of the cluster in the list of
        clusters reported.

    Returns
    -------
    str
        Newline delimited representation of the cluster. Each line is a tab
        delimited representation of a ReadPair, prefixed by `cluster_id` and
        the size of the cluster.
    """
    return '\n'.join([cluster_id + '\t' +
                      str(len(cluster)) + '\t' +
                      pair.as_rpc() for pair in cluster])


def rpc(pairs, cutoff, min_mapq, min_size, excluded_regions, excluded_file):
    """
    Clusters read pairs.

    Parameters
    ----------
    fin : file or pysam.AlignmentFile
        File of ReadPairs in bamstat format. (samfile support coming soon)
    cutoff : int
        Maximum clustering distance.
    min_mapq : int
        Minimum mapping quality for a read to be included.
    min_size : int
        Minimum size for a cluster to be reported.
    excluded_regions : pysam.TabixFile
        List of regions to exclude from clustering.

    Yields
    ------
    set of ReadPair
        A cluster of read pairs.
    """
    for candidates in get_candidate_pairs(pairs, cutoff, excluded_regions,
                                          excluded_file):
        clusters = cluster_candidates(candidates, cutoff, min_mapq, min_size)

        # Sort clusters internally by first read's position,
        # then sort clusters by first pair's first read's position
        clusters = [sorted(c, key=lambda p: (p.read1.reference_start, p.read1.query_name)) for c in clusters]
        clusters = sorted(clusters, key=lambda c: c[0].read1.reference_start)

        for cluster in clusters:
            yield cluster


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # TODO: allow bam input
    # group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('bam', type=pysam.Samfile,
                       # help='Discordant pairs reported by samblaster')
    # group.add_argument('fin', type=argparse.FileType('r'),
                       # default=sys.stdin,
                       # help='Bamstat formatted input file.')

    parser.add_argument('fin', type=open,
                        default=sys.stdin,
                        help='Bamstat formatted input file.')
    parser.add_argument('fout', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='RPC formatted output file.')

    parser.add_argument('-x', '--exclude', dest='excluded_regions',
                        metavar='BED.GZ',
                        type=pysam.TabixFile, default=None,
                        help='Tabix indexed bed file of regions to exclude '
                        'from clustering.')

    parser.add_argument('--save-exclude', dest='excluded_file',
                        metavar='TXT',
                        type=argparse.FileType('w'), default=None,
                        help='File to save excluded pairs in bamstat format.')

    # TODO: compute distance as median insert + 7 * MAD
    # parser.add_argument('insert_metrics')
    parser.add_argument('-d', '--dist', metavar='CUTOFF', type=int,
                        default=3000,
                        help='Max distance between reads for two pairs to be '
                        'clustered together.')

    parser.add_argument('-q', '--qual', metavar='MAPQ', type=int,
                        default=-1,
                        help='Minimum quality score for a read to be '
                        'included.')

    # TODO: compute 7% haploid by default; allow user specification
    # parser.add_argument('coverage_metrics')
    parser.add_argument('-s', '--size', metavar='PAIRS', type=int,
                        default=3,
                        help='Minimum number of read pairs for a cluster to '
                        'be reported.')
    
    # almost certainly not worth the trouble
    # parser.add_argument('-c', '--cluster', metavar='ALGORITHM',
                        # choices=['SLINK', 'DBSCAN'],
                        # default='SLINK',
                        # help='Clustering algorithm to use.')

    args = parser.parse_args()

    i = 1
    pairs = ReadPairParser(args.fin)
    for cluster in rpc(pairs, args.dist, args.qual, args.size,
                       args.excluded_regions, args.excluded_file):
        args.fout.write(rpc_format(cluster, str(i)) + '\n\n')
        i += 1


if __name__ == '__main__':
    main()
