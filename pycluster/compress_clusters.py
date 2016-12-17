#!/usr/bin/env python


import argparse
import os
import glob
from pycluster import cluster as cl, parse_clusters


def compress_clusters(clusterfile, bedpe, prefix, samples=None, split=True):
    for cluster_list in parse_clusters(clusterfile, prefix):
        if split:
            for cluster in cluster_list:
                if samples is None:
                    bedpe.write("%s\n" % str(cluster))
                else:
                    bedpe.write('%s\n' % cluster.str_with_samples(samples))
        else:
            cluster = cl.merge(cluster_list)
            if samples is None:
                bedpe.write("%s\n" % str(cluster))
            else:
                bedpe.write('%s\n' % cluster.str_with_samples(samples))


def main():
    parser = argparse.ArgumentParser(
        description="Compresses each cluster in readpaircluster output to a "
        "single line in a bedpe file.")
    parser.add_argument('rpc', type=argparse.FileType('r'),
                        help="Input readpaircluster file")
    parser.add_argument('bedpe', type=argparse.FileType('w'),
                        help="Output bedpe file")
    parser.add_argument('prefix', help="Prefix for cluster IDs. Each cluster "
                        "will be assigned an ID consisting of this prefix "
                        "plus the rpc cluster ID number. "
                        "Generally SAMPLE_SVTYPE.")
    parser.add_argument('-s', '--samples', type=argparse.FileType('r'),
                        default=None)
    parser.add_argument('--no-split', help="Don't split clusters",
                        dest='split', action='store_false', default=True)
    args = parser.parse_args()

    if args.samples is not None:
        samples = args.samples.read().rstrip().split('\n')

        header = ('#chrA\tstartA\tendA\tchrB\tstartB\tendB\tname\tscore\t'
                  'strandA\tstrandB\tsize\tqualA\tqualB\tuniqA\tuniqB\t' +
                  'num_samples\tnum_probands\t' + '\t'.join(samples))
        args.bedpe.write(header + '\n')

    else:
        samples = None

    compress_clusters(args.rpc, args.bedpe, args.prefix, samples, args.split)


if __name__ == '__main__':
    main()
