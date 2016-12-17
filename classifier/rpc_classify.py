#!/usr/bin/env python


from __future__ import division
import argparse
import sys
import pysam
from pycluster import cluster as cl, parse_clusters
# import pdb

def rpc_classify(cfile, prefix, insert_handles, gcov_dict,
                 cluster_bedpe, event_bedpe):
    """
    Classify a set of rpcluster output.

    Args:
        cfile (file)          : rpcluster ouput file
        prefix (str)          : Prefix to include in all breakpoint IDs,
                                generally cohort name
        insert_handles (dict) : Insert coverage files
                                {sampleID (str) : icov_f (pysam.Tabixfile)}
        gcov_dict (dict)      : Global coverage
                                {sampleID (str) : gcov (int)}
        cluster_bedpe (file)  : Output file for raw breakpoint clusters
        event_bedpe (file)    : Output file for merged breakpoints (events)
    """

    # contigs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
               # 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               # 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
               # 'chrX', 'chrY']
    # contigs = contigs + [contig[3:] for contig in contigs]

    # pdb.set_trace()
    for cluster_list in parse_clusters(cfile, prefix):
        # We currently do not have insert coverage over unplaced/random
        # contigs. Skip any clusters that map to one.
        # cltest = cluster_list[0]
        # if cltest.chrA not in contigs or cltest.chrB not in contigs:
            # continue

        # mapq and uniq are calculated in cluster instantiation
        # Calculate local/global coverage, then write out
        max_gcov = 0
        max_localA = 0
        max_localB = 0

        for cluster in cluster_list:
            gcov = cluster.global_cov(gcov_dict)
            localA, localB = cluster.local_cov(insert_handles)
            cluster_bedpe.write("%s\t%.3f\t%.3f\t%.3f\t%d\t[%s]\n" %
                                (str(cluster), gcov, localA, localB,
                                 len(cluster.samples.keys()),
                                 ",".join(cluster.samples.keys())))

            if gcov > max_gcov:
                max_gcov = gcov
            if localA > max_localA:
                max_localA = localA
            if localB > max_localB:
                max_localB = localB

        event = cl.merge(cluster_list)
        call, branch = event.classify(max_gcov, max_localA, max_localB)
        event_bedpe.write("%s\t%.3f\t%.3f\t%.3f\t%d\t[%s]\t%s\t%d\n" %
                          (str(event), max_gcov, max_localA, max_localB,
                           len(event.samples.keys()),
                           ",".join(event.samples.keys()),
                           call, branch))


def main():
    parser = argparse.ArgumentParser(
        description="SV classifier. Computes mapq, uniqueness, global cov, "
        "and local cov for every cluster output by rpcluster. Decision tree "
        "then classifies based on these metrics. Outputs to bedpe format.")
    parser.add_argument('rpcluster', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin,
                        help="rpCluster output file [stdin]")

    parser.add_argument('event_bedpe', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="Event output file [stdout]")

    parser.add_argument('sample_list', type=argparse.FileType('r'),
                        help="List of sample IDs. Column format: "
                        "sample_ID bstat_dir cl_dist cutoff "
                        "global_cov insert_cov_file")

    parser.add_argument('--cluster-bedpe', type=argparse.FileType('w'),
                        help="Cluster output file")

    parser.add_argument('cluster_name',
                        help="Prefix for cluster IDs")

    args = parser.parse_args()

    # Skip header
    # TODO: add --no-header flag
    # args.sample_list.next()

    cluster_header = ("chrA\tstartA\tendA\tchrB\tstartB\tendB\t"
                      "name\tscore\tstrandA\tstrandB\t"
                      "size\tmapqA\tmapqB\tuniqA\tuniqB\t"
                      "global\tlocalA\tlocalB\tpoly_count\tsamples\n")

    event_header = ("%s\tcall\tbranch\n" % cluster_header.rstrip())

    args.cluster_bedpe.write(cluster_header)
    args.event_bedpe.write(event_header)

    # Create coverage dictionaries
    # gcov_dict = {sample_ID: int(sample_haploid_coverage)}
    # insert_handles = {sample_ID: pysam.Tabixfile(sample_insert_coverage)}
    gcov_dict = {}
    insert_handles = {}
    for line in args.sample_list:
        if line.startswith("#"):
            continue

        sample = line.split()[0]
        gcov, icov_f = line.split()[4:6]

        gcov_dict[sample] = int(gcov)
        insert_handles[sample] = pysam.Tabixfile(icov_f)

    rpc_classify(args.rpcluster, args.cluster_name, insert_handles,
                 gcov_dict, args.cluster_bedpe, args.event_bedpe)


if __name__ == '__main__':
    main()
