#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""
Convert a coordinate sorted bam to readPairCluster's input format.

Designed for use with the discordant bam produced by samblaster. Assumes all
pairs present are discordant, and produces four output files based on pair
orientation (del: +/-, dup: -/+, inv: +/+ or -/-, tloc: mismatched
chromosomes). Assumes no multiply mapped reads, i.e., each pair ID appears only
twice.
"""

import argparse
import logging
from collections import deque, defaultdict
from itertools import combinations
import pysam


def get_svtype(read):
    """
    Classify type of discordant read based on orientation.

    Parameters
    ----------
    read : pysam.AlignedRead
        Read to classify

    Returns
    -------
    svtype : str
        One of 'del', 'dup', 'inv', or 'tloc'
    """

    if read.tid != read.rnext:
        return 'tloc'
    if read.is_reverse == read.mate_is_reverse:
        return 'inv'
    if read.pos < read.pnext:
        if read.is_reverse:
            return 'dup'
        else:
            return 'del'
    else:
        if read.is_reverse:
            return 'del'
        else:
            return 'dup'


def rpc_format(read1, read2, samfile):
    """
    Merge two Pysam reads into a line formatted for readPairCluster.

    Parameters
    ----------
    read1 : pysam.AlignedRead
        First read in pair

    read2 : pysam.AlignedRead
        Second read in pair

    samfile : pysam.Samfile
        Samfile that the reads came from. Necessary to convert each read's
        reference tid back to a chromosome name.

    Returns
    -------
    rpcformat : str
        A line of readPairCluster input. Does not include trailing newline.

    Raises
    ------
    Exception
        If the two reads came from different pairs.
    """

    if read1.qname != read2.qname:
        raise Exception('Trying to merge two reads from different pairs')

    if isinstance(read1.qual, bytes):
        qual1 = read1.qual.decode()
        qual2 = read2.qual.decode()
    elif isinstance(read1.qual, str):
        qual1 = read1.qual
        qual2 = read2.qual
    else:
        msg = 'Unrecognized pysam quality type %s' % type(read1.qual)
        raise Exception(msg)

    rpcformat = [read1.qname,
                 samfile.getrname(read1.tid),
                 read1.pos,
                 read1.flag,
                 samfile.getrname(read2.tid),
                 read2.pos,
                 read2.flag,
                 read1.rlen,
                 read1.mapq,
                 read2.rlen,
                 read2.mapq,
                 read1.seq,
                 qual1,
                 read2.seq,
                 qual2]

    return('\t'.join([str(x) for x in rpcformat]))

# Instead, just return pairs. Categorize on other end
# def SamPairParser(samfile):
    # """
    # Find mates of each read in coordinate sorted samfile and return pair,
    # maintaining coordinate sort order of first read in pair.

    # Parameters
    # ----------
    # samfile : pysam.Samfile
        # Samfile to parse

    # Yields
    # ------
    # (read1, read2) : (pysam.AlignedSegment, pysam.AlignedSegment)
        # Pair of reads (sorted by read1 coordinate position)
    # """

    # # Read IDs are stored in a deque to maintain sort order of input
    # read_ids = deque()

    # # Reads themselves are stored in a dictionary, keyed by their IDs
    # reads = defaultdict(list)

    # # At end of each chromosome, iterate through deque, yielding pairs
    # curr_tid = None

    # # Helper - yield all cached pairs
    # def yield_curr_reads():
        # unmatched_ids = deque()

        # for read_id in read_ids:
            # try:
                # read1, read2 = reads[read_id]
                # yield read1, read2

                # # Clean up
                # reads.pop(read_id)

            # # Unmatched read IDs at the end of a chromosome should be
            # # translocation pairs - save for possible mates on next chromosome
            # except ValueError:
                # if len(reads[read_id]) == 1:
                    # unmatched_ids.append(read_id)
                # else:
                    # logging.warning('Read ID %s corresponds to more than '
                                    # 'two sequence alignments. Omitting '
                                    # 'from output.' % read_id)
        # return unmatched_ids

    # for read in samfile:
        # # Initialize chromosome tracking
        # if curr_tid is None:
            # curr_tid = read.tid

        # # Write out all cached pairs after completing each chromosome
        # if read.tid != curr_tid:
            # read_ids = yield_curr_reads()
            # curr_tid = read.tid

        # # Exclude pairs whose reads map to the same position
        # if read.tid == read.rnext and read.pos == read.pnext:
            # continue

        # # If first read in pair, add its ID to the queue
        # if len(reads[read.qname]) == 0:
            # read_ids.append(read.qname)

        # # Save both reads in pair to dict
        # reads[read.qname].append(read)

    # yield_curr_reads()


# Instead, just return pairs. Categorize on other end
def parse_pairs(samfile):
    """
    Find mates of each read in samfile and segregate by svtype.

    Maintains coordinate sort order of first read in pair.

    Parameters
    ----------
    samfile : pysam.Samfile
        Samfile to convert

    Returns
    -------
    (deques, dicts) : (dict, dict)
        Tuple of read pair data for each chromosome.

        Each dictionary has the four svtypes as keys. deques[svtype] is a
        collections.deque containing read pair IDs in coordinate sorted order.
        dicts[svtype] is a dictionary linking these IDs to a pair of
        pysam.AlignedRead objects. To parse, popleft() IDs out of the deque and
        access read pairs in the dictionary via the ID.

        deques['tloc'] and dicts['tloc'] are both dictionaries with all
        chromosome reference id's as keys. deques['tloc'][tid] and
        dicts['tloc'][tid] have the same structure as deques[svtype] and
        dicts[svtype]
    """

    # Read IDs are stored in a deque to maintain sort order of input
    deques = {}
    deques['del'] = deque()
    deques['dup'] = deque()
    deques['inv'] = deque()

    # Reads themselves are stored in a dictionary, keyed by their IDs
    dicts = {}
    dicts['del'] = defaultdict(list)
    dicts['dup'] = defaultdict(list)
    dicts['inv'] = defaultdict(list)

    # Store separate deque and dict for each chromosome combination
    deques['tloc'] = {}
    dicts['tloc'] = {}
    for tid_pair in combinations(range(samfile.nreferences), 2):
        deques['tloc'][tid_pair] = deque()
        dicts['tloc'][tid_pair] = defaultdict(list)

    curr_tid = None

    for read in samfile:
        # Initialize chromosome tracking
        # if curr_tid is None:
            # curr_tid = read.tid

        # TODO: Write out all cached pairs after completing each chromosome
        # For tlocs: write out all combinations of (0..(curr - 1)) and curr
        # if read.tid != curr_tid:
            # yield deques, dicts
            # curr_tid = read.tid

        # Exclude pairs that map to the same position
        if read.tid == read.rnext and read.pos == read.pnext:
            continue

        svtype = get_svtype(read)

        if svtype != 'tloc':
            curr_deque = deques[svtype]
            curr_dict = dicts[svtype]
        else:
            tid_pair = tuple(sorted([read.tid, read.rnext]))
            curr_deque = deques[svtype][tid_pair]
            curr_dict = dicts[svtype][tid_pair]

        # If first read in pair, add its ID to the queue
        if len(curr_dict[read.qname]) == 0:
            curr_deque.append(read.qname)

        # Save both reads in pair to dict
        curr_dict[read.qname].append(read)

    return deques, dicts


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('samfile', type=pysam.Samfile)
    parser.add_argument('prefix')
    parser.add_argument('-n', '--name', default=None)
    args = parser.parse_args()

    logging.basicConfig(filename='disc_to_rpc.log', level=logging.INFO)

    if args.name is None:
        name = ''
    else:
        name = '.' + args.name

    fouts = {}
    fouts['del'] = open('%s.del.pairs%s.txt' % (args.prefix, name), 'w')
    fouts['dup'] = open('%s.dup.pairs%s.txt' % (args.prefix, name), 'w')
    fouts['inv'] = open('%s.inv.pairs%s.txt' % (args.prefix, name), 'w')
    fouts['tloc'] = open('%s.tloc.pairs%s.txt' % (args.prefix, name), 'w')

    deques, dicts = parse_pairs(args.samfile)
    for svtype in ['del', 'dup', 'inv']:
        while len(deques[svtype]) > 0:
            read_id = deques[svtype].popleft()
            try:
                read1, read2 = dicts[svtype][read_id]
            except ValueError:
                if len(dicts[svtype][read_id]) == 1:
                    logging.warning('Read ID %s does not have a mate. '
                                    'Omitting from output.' % read_id)
                else:
                    logging.warning('Read ID %s corresponds to more than two '
                                    'sequence alignments. Omitting from output.'
                                    % read_id)
                continue

            # Filter pairs with mapq=0 on both sides
            if read1.mapping_quality == 0 and read2.mapping_quality == 0:
                continue

            line = rpc_format(read1, read2, args.samfile)
            fouts[svtype].write(line + '\n')

    svtype = 'tloc'
    for tid_pair in combinations(range(args.samfile.nreferences), 2):
        while len(deques[svtype][tid_pair]) > 0:
            read_id = deques[svtype][tid_pair].popleft()

            try:
                read1, read2 = dicts[svtype][tid_pair][read_id]
            except ValueError:
                if len(dicts[svtype][tid_pair][read_id]) == 1:
                    logging.warning('Read ID %s does not have a mate. '
                                    'Omitting from output.' % read_id)
                else:
                    logging.warning('Read ID %s corresponds to more than two '
                                    'sequence alignments. Omitting from output.'
                                    % read_id)
                continue

            line = rpc_format(read1, read2, args.samfile)
            fouts[svtype].write(line + '\n')


if __name__ == '__main__':
    main()
