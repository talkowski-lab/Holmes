#!/usr/bin/env python


import argparse
import heapq
import filter_sample_clusters as fsc


class HeapIter():
    """
    Sortable class for filter_iterators
    """

    def __init__(self, read, iterator):
        self.read = read
        self.iterator = iterator

        self.contig = self.read.split()[1]
        self.pos = int(self.read.split()[2])
        self.contigB = self.read.split()[4]
        self.posB = int(self.read.split()[5])

    def __lt__(self, other):
        if self.contig == other.contig:
            if self.pos != other.pos:
                return self.pos < other.pos
            else:
                if self.contigB == other.contigB:
                    return self.posB < other.posB
                else:
                    return self.contigB < other.contigB
        else:
            return self.contig < other.contig

    def next(self):
        try:
            self.read = self.iterator.next()
            self.contig = self.read.split()[1]
            self.pos = int(self.read.split()[2])
        except StopIteration:
            self.read = None


def filter_iterator(cfile, cutoff, sample, tloc=False):
    """
    Wrapper for filter_sample_clusters to produce iterator over reads

    Args:
        cfile (file) : open rpcluster file
        cutoff (int) : minimum cluster size cutoff
        sample (str) : Sample ID to prepend to read IDs
        tloc (bool)  : cfile contains translocation cluster data. Reads in
            translocation files are not sorted relative to the file, and
            need to be handled differently.

    Yields:
        str : Next read from file in sorted order
    """

    # Translocation cluster files do not have reads in sorted order,
    # just sorted within each cluster. Unfortunately there is no
    # consistent ordering of the clusters, and there is no guarantee
    # that the maximum or minimum position in a cluster will be less than
    # the maximum or minimum of all following clusters. There is likely
    # a more elegant solution, but the following solution is still more
    # time and space efficient than sorting.
    # Tested on LOGIC translocation clusters:
    #   This implementation (filter + merge): 1m44s, 16MB mem, 262MB swap
    #   Unix sort of filtered reads: 13m57s, 869MB mem, 5.0GB swap
    #   Matlab filter and sort: 37m15s, 870MB mem, 5.1GB swap

    # For each contig, place all pairs with read A mapping to contig in a
    # priority queue, sorted by position. The input file maintains sort by
    # contig, so all pairs are consecutive in the file. Once all reads for the
    # contig have been read in, yield all, then reset for next contig
    if tloc:
        iterator = fsc.filter_sample_clusters(cfile, cutoff, sample)
        cluster = iterator.next()
        curr_contig = cluster[-1].split()[1]

        heap = []
        for read in cluster:
            heapq.heappush(heap, (curr_contig, int(read.split()[2]), read))

        for cluster in iterator:
            contig = cluster[0].split()[1]

            # Keep reads on current contig sorted in heap
            if contig == curr_contig:
                for read in cluster:
                    heapq.heappush(heap,
                                   (curr_contig, int(read.split()[2]), read))
            else:
                # Output all reads on current contig in order
                while len(heap) > 0:
                    contig, pos, read = heapq.heappop(heap)
                    yield read

                # Start the heap over for the new contig
                curr_contig = contig
                for read in cluster:
                    heapq.heappush(heap,
                                   (curr_contig, int(read.split()[2]), read))
            while len(heap) > 0:
                contig, pos, read = heapq.heappop(heap)
                yield read

    else:
        for cluster in fsc.filter_sample_clusters(cfile, cutoff, sample):
            for read in cluster:
                yield read


def merge_clusters(svtype, samples):
    """
    Merge filtered readpaircluster output for one sv type across all samples.

    Args:
        svtype (str)   : One of deletion, insertion, inversion, or transloc
        samples (list) : List of (sample, bamstat_dir, cluster_dist, cutoff)
                         tuples

    Yields:
        str: The next read in sorted (chr, pos) order across all files
    """

    # pdb.set_trace()
    # Generate list of (sample, rpcluster_output, size_cutoff) tuples
    cfiles = [(smpl[0],
               "%s/%s_clusters_d%s_q-1.txt" % (smpl[1], svtype, smpl[2]),
               smpl[3])
              for smpl in samples]

    # Create an iterator for every sample that will return reads
    # from filtered clusters in sorted order
    iters = [filter_iterator(open(cfile), int(cutoff), sample,
                             svtype == "transloc") for
             sample, cfile, cutoff in cfiles]

    # Instantiate the heap with the first read in each sample
    heap = []
    for cluster_iter in iters:
        read = cluster_iter.next()
        heapq.heappush(heap, HeapIter(read, cluster_iter))

    # Pull each read off the heap and add the next read from that sample
    # Heap remains sorted by chromosome (lexicographic) then position (numeric)
    while heap:
        heapiter = heapq.heappop(heap)
        yield heapiter.read

        heapiter.next()
        if not heapiter.read:
            continue

        heapq.heappush(heap, heapiter)


def main():
    parser = argparse.ArgumentParser(
        description="Filter and merge rpcluster output from multiple samples.")
    parser.add_argument(
        'sample_list', type=argparse.FileType('r'),
        help="List of samples to merge. Column format: sample_ID "
        "bamstat_directory cluster_distance cluster_size_cutoff")
    parser.add_argument(
        'prefix', help="Prefix of each merged readpair cluster output. "
        "Output will be written to [prefix]_[svtype].txt")
    parser.add_argument('--no-header', dest='header', action='store_false',
                        help="Sample list does not have a header "
                        "[default: True]")
    parser.add_argument('svtype', choices=['deletion', 'insertion',
                                           'inversion', 'transloc'],
                        help="Type of rpcluster output to merge. ")
    parser.set_defaults(header=True)
    args = parser.parse_args()

    if args.header:
        args.sample_list.next()

    samples = [line.rstrip().split() for line in args.sample_list]

    # Iterate over svtypes here or take svtype as argument?
    # svtypes = ['deletion', 'insertion', 'inversion', 'transloc']

    # svtype = 'deletion'
    fout = open('%s_%s.reads' % (args.prefix, args.svtype), 'w')
    for read in merge_clusters(args.svtype, samples):
        fout.write("%s\n" % read)


if __name__ == '__main__':
    main()
