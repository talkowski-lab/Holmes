#!/usr/bin/env python


import argparse
import re
import sys
import pybedtools as pbt
# import pdb


def filter_sample_clusters(fin, size_cutoff=0, sample="", reformat=True,
                           excluded=pbt.BedTool("", from_string=True),
                           contigs=[], delim="___"):
    """
    Filter clusters for mapq=0 pairs and cluster size

    Args:
        fin (readable file) : readPairCluster output to filter
        size_cutoff (int) : Minimum cluster size required
        sample (string) : Sample ID to prepend to read IDs
        reformat (bool) : If true, output will be reformatted to
            bamstat output format. If false, output will remain in
            readPairCluster output format.
        excluded (BedTool) :
        contigs (list)

    Yields:
        list : List of reads in a cluster
    """

    empty_exp = re.compile(r'^\s*$')
    bamstat_format = [2, 3, 5, 11, 4, 6, 12, 9, 7, 10, 8, 13, 15, 14, 16]

    # Underscore only needed to separate sample ID if ID exists
    if sample:
        sample = "%s%s" % (sample, delim)

    # pdb.set_trace()

    cluster = []
    for line in fin:
        # Clusters are delimited by blank lines
        if not empty_exp.match(line):
            cl_data = line.rstrip().split()

            chrA, chrB = cl_data[3:5]
            posA, posB = [int(x) for x in cl_data[5:7]]
            mapqA, mapqB = [int(x) for x in cl_data[7:9]]
            lenA, lenB = [int(x) for x in cl_data[9:11]]

            # pdb.set_trace()

            cl_intA = pbt.Interval(chrA, posA, posA + lenA - 1)
            cl_intB = pbt.Interval(chrB, posB, posB + lenB - 1)

            # Skip
            if chrA not in contigs or chrB not in contigs:
                continue

            # Skip
            if excluded.any_hits(cl_intA) or excluded.any_hits(cl_intB):
                continue

            # Skip if read and mate map to identcal position
            if chrA == chrB and posA == posB:
                continue

            # Remove pairs where mapq=0 on both sides
            if mapqA == 0 and mapqB == 0:
                continue
                # data = line.rstrip().split()

            # output reads that pass all filters
            if reformat:
                read = "\t".join([cl_data[i] for i in bamstat_format])
                read = "%s%s" % (sample, read)
            else:
                read = "\t".join(cl_data[0:2] +
                                 ["%s%s" % (sample, cl_data[2])] +
                                 cl_data[4:])
            cluster.append(read)

        else:
            # Don't output clusters smaller than the size_cutoff
            if len(cluster) >= size_cutoff:
                yield cluster
            cluster = []


def main():
    parser = argparse.ArgumentParser(
        description="Filters readPairCluster output to remove "
        "clusters with mapq=0 on both sides and clusters with "
        "size below a minimum threshold")
    parser.add_argument('fin', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin,
                        help="rpCluster output file to filter "
                        "[default: stdin]")
    parser.add_argument('fout', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="filtered rpCluster output "
                        "[default: stdout]")
    parser.add_argument('sample', help="Sample ID to prepend to read IDs")
    parser.add_argument('-c', '--cutoff', type=int, default=4,
                        help="Minimum cluster size to allow, "
                        "[default: 4]")
    parser.add_argument('--no-reformat', dest='reformat', action='store_false',
                        help="Don't reformat output to bamstat format")
    parser.add_argument('--exclude',
                        help="Genomic regions to exclude from clustering "
                        "(BED format) [None]")
    parser.add_argument('--contigs', type=argparse.FileType('r'),
                        help="Contigs allowed [autosomes and sex chromosomes]")
    parser.set_defaults(reformat=True)
    args = parser.parse_args()

    if args.exclude is not None:
        # excluded = [line.strip() for line in args.exclude]
        excluded = pbt.BedTool(args.exclude)
    else:
        excluded = pbt.BedTool("", from_string=True)

    if args.contigs is not None:
        contigs = [line.strip() for line in args.contigs]
    # Allow chrN and N for autosomes, and chrX/chrY and X/Y for sex
    else:
        contigs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
                   'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                   'chr21', 'chr22', 'chrX', 'chrY']

        contigs = contigs + [contig[3:] for contig in contigs]

    for cluster in filter_sample_clusters(args.fin, args.cutoff,
                                          args.sample, args.reformat,
                                          excluded, contigs):
        args.fout.write("\n".join(cluster) + "\n")
        if not args.reformat:
            args.fout.write("\n")


if __name__ == '__main__':
    main()
