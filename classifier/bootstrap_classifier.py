#!/usr/bin/env python


import argparse
import os
import sys
import re
import glob
import pandas as pd
from pycluster import scrape_bamstat


def bootstrap_classify(sample_list, fout, header=False, metrics=None,
                       cutoff_t=0.07, cutoff_min=4, mad=7):
    """
    Get cluster_distance, cluster_cutoff, and global_coverage for a sample
    """

    # Skip header
    if header:
        sample_list.next()

    # header_fmt = "%s\t" * 4 + "%s\n"
    # header = header_fmt % ('#sample', 'bamstat_directory', 'cluster_distance',
                           # 'cluster_cutoff', 'global_coverage')
    # fout.write(header)

    fout_fmt = "%s\t%s\t%d\t%d\t%d\n"

    if metrics is not None:
        metrics['cluster_dist'] = metrics.Median_Insert + mad * metrics.Insert_MAD

    for line in sample_list:
        sample, bstat_dir = line.rstrip().split()[0:2]

        # Calculate cluster_cutoff as 7% of global haploid coverage
        if metrics is not None:
            coverage = metrics['Hap_Phys_Cov'][sample]
            cluster_dist = metrics['cluster_dist'][sample]
        elif os.path.isfile('%s/stats.file' % bstat_dir):
            statfile = open("%s/stats.file" % bstat_dir)
            coverage, _, _, _, _ = scrape_bamstat(statfile)
            coverage = int(coverage)

            # TODO: use del_size returned by scrape_bamstat()
            dist_exp = re.compile(r'_d(\d+)_')
            del_fname = glob.iglob("%s/*deletion_*q-1*.txt" % bstat_dir).next()
            cluster_dist = int(dist_exp.search(del_fname).group(1))
        else:
            raise Exception('No stats file in bamstat directory and '
                            'no metrics table provided.')

        # Do we want to round or take ceiling?
        cutoff = max(int(round(cutoff_t * coverage)), cutoff_min)


        # Currently writing out to file
        # TODO: just call merge_clusters
        fout.write(fout_fmt % (sample, bstat_dir, cluster_dist, cutoff,
                               coverage))


def main():
    parser = argparse.ArgumentParser(
        description='Generate input file for run_classify.sh')
    parser.add_argument('sample_list', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin,
                        help="List of sample IDs and corresponding "
                        "bamstat directories. [stdin]")
    parser.add_argument('fout', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="Output file [stdout]")
    parser.add_argument('--header', action='store_true',
                        help="Sample list has header [false]")
    parser.add_argument('--metrics', '-m', type=argparse.FileType('r'),
                        default=None,
                        help='Table of sample metrics')
    parser.set_defaults(header=False)
    args = parser.parse_args()

    if args.metrics is not None:
        header = -1
        for line in args.metrics:
            if line.startswith('#'):
                header += 1
            else:
                break

        metrics = pd.read_table(args.metrics.name, header=header, index_col=0)
    else:
        metrics = None

    bootstrap_classify(args.sample_list, args.fout, args.header, metrics)


if __name__ == '__main__':
    main()
