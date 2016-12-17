#!/usr/bin/env python


import argparse
import re


# Taken from http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics
HG19_NON_N_GENOME_SIZE = 2897310462

def scrape_bamstat(statsfile, threshold=7):
    proper_exp = re.compile(r'Proper pairs\s*(\d+)')
    insert_exp = re.compile(r'Actual FR median insert size:\s*(\d+)')
    dev_exp = re.compile(r'Actual FR median absolute deviation:\s*(\d+)')

    for line in statsfile:
        if proper_exp.match(line):
            proper = int(proper_exp.match(line).group(1))
        if insert_exp.match(line):
            insert = int(insert_exp.match(line).group(1))
        if dev_exp.match(line):
            dev = int(dev_exp.match(line).group(1))

    coverage = proper * insert / HG19_NON_N_GENOME_SIZE
    del_size = insert + threshold * dev

    return coverage, proper, insert, dev, del_size


def main():
    parser = argparse.ArgumentParser(
        description="Script that scrapes sample directories in Samples/ndd "
        "for bamstat stats files and writes a file of average library "
        "coverage, number of proper pairs in the library, median insert of "
        "the library, and the median absolute deviation from this median.")
    parser.add_argument('samples', type=argparse.FileType('r'),
                        help="Tab separated file containing sample names "
                        "and their subdirectories under "
                        "/data/talkowski/Samples/ndd")
    parser.add_argument('outfile', type=argparse.FileType('w'),
                        help="Output file to write to. File will be tab "
                        "separated and be in the format: Sample Coverage "
                        "Proper_Pair_Count Median_Insert Deviation")
    args = parser.parse_args()

    args.outfile.write('sample\tcoverage\tproper\tinsert\tdev\tdel_size\n')

    for line in args.samples:
        sample, bstat_dir = line.rstrip().split()[0:2]
        statsfile = open('%s/stats.file' % bstat_dir)

        coverage, proper, insert, dev, del_size = scrape_bamstat(statsfile)
        args.outfile.write("%s\t%d\t%d\t%d\t%d\t%d\n" %
                           (sample, coverage, proper, insert, dev, del_size))


if __name__ == '__main__':
    main()
