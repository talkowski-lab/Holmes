#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import heapq
from collections import deque
import pysam
from svcf.rpc import RPC
from svcf.svcall import SVCallCluster
from svcf.svfile import SVFile


class FileFormatError(Exception):
    """Improper file format"""


class SVCluster(RPC):
    def __init__(self, svfiles, dist=500, excluded_regions=None, svtype=None):
        nodes = heapq.merge(*svfiles)

        self.samples = sorted(set([svfile.sample for svfile in svfiles]))
        self.sources = sorted(set([svfile.source for svfile in svfiles]))

        # TODO: add size filter
        def _filter(nodes):
            for node in nodes:
                if svtype:
                    node.svtype = svtype
                if (not node.secondary) and (node.is_allowed_chrom):
                    yield node

        super().__init__(_filter(nodes), dist, 1, excluded_regions)

    def cluster(self):
        for cluster in super().cluster():
            cluster = [call for call in cluster if not call.secondary]
            if len(cluster) > 0:
                yield SVCallCluster(self.samples, self.sources, cluster)


def parse_filelist(filelist):
    """
    filelist : File
    """

    svfiles = deque()
    samples = set()
    sources = set()

    for line in filelist:
        data = line.strip().split()
        if len(data) != 3:
            msg = 'Expected three fields in filelist, found %d' % len(data)
            raise FileFormatError(msg)

        filename, sample, source = data

        svfiles.append(SVFile(filename, sample=sample, source=source))
        samples.add(sample)
        sources.add(source)

    return svfiles, sorted(samples), sorted(sources)


def main():
    parser = argparse.ArgumentParser(
        description='Intersect SV called by cluster-based algorithms.')
    parser.add_argument('filelist',
                        type=argparse.FileType('r'),
                        help='List of SV file, sample, source program. '
                        'Tab delimited.')
    parser.add_argument('--svcf',
                        type=argparse.FileType('w'), default=None,
                        help='Output SVCF.')
    parser.add_argument('--links',
                        type=argparse.FileType('w'), default=None,
                        help='Output links file. Connects merged IDs back to '
                        'variant IDs from original SV files.')
    parser.add_argument('--inherit',
                        type=argparse.FileType('w'), default=None,
                        help='Output inheritance file.')
    parser.add_argument('-t', '--svtype',
                        default=None,
                        help='Specify SV type of input files instead of '
                        'attempting to parse. Requires all files in list to '
                        'be of same svtype.')
    parser.add_argument('-d', '--dist',
                        type=int, default=500,
                        help='Maximum clustering distance. Suggested to use '
                        'max of median + 7*MAD over samples. [500]')
    parser.add_argument('-p', '--prefix',
                        default='MERGED',
                        help='Prefix for merged variant IDs. [MERGED]')
    parser.add_argument('-x', '--excluded', metavar='BED.GZ',
                        type=pysam.TabixFile, default=None,
                        help='Tabix indexed bed of blacklisted regions. Any '
                        'SV with a breakpoint falling inside one of these '
                        'regions is filtered from output.')
    parser.add_argument('--no-header', dest='header',
                        action='store_false', default=True,
                        help="Don't print file headers")
    args = parser.parse_args()

    svfiles, samples, sources = parse_filelist(args.filelist)

    inherit_header = '\t'.join('chrA posA posB chrB name svtype freq sample source status'.split())

    svcf_header = '\t'.join('chrA posA posB chrB name svtype freq samples FORMAT total'.split())
    svcf_header = svcf_header + '\t' + '\t'.join(samples)

    links_header = '\t'.join('variant svtype sample source source_ID source_svtype'.split())

    if args.header:
        args.inherit.write(inherit_header + '\n')
        args.svcf.write(svcf_header + '\n')
        args.links.write(links_header + '\n')

    svc = SVCluster(svfiles, dist=args.dist, excluded_regions=args.excluded)
    for i, cluster in enumerate(svc.cluster()):
        # TODO: auto add chromosome to name
        cluster.name = args.prefix + '_' + str(i + 1)

        if args.inherit:
            inherit_str = cluster.inherit
            if inherit_str:
                args.inherit.write(inherit_str + '\n')

        if args.svcf:
            args.svcf.write(cluster.svcf + '\n')
        if args.links:
            args.links.write(cluster.links + '\n')


if __name__ == '__main__':
    main()
