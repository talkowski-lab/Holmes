

from __future__ import division
import sys
import numpy as np
import pdb


class Cluster(object):
    def __init__(self, data, reads=None, samples=None):
        self.chrA = data[0]
        self.startA = int(data[1])
        self.endA = int(data[2])
        self.chrB = data[3]
        self.startB = int(data[4])
        self.endB = int(data[5])
        self.name, self.score, self.strandA, self.strandB = data[6:10]
        self.size = int(data[10])
        self.qualA = float(data[11])
        self.qualB = float(data[12])
        self.uniqA = float(data[13])
        self.uniqB = float(data[14])
        self.reads = reads
        self.samples = samples

        if self.reads:
            self.qualA, self.qualB = self._mapq()
            self.uniqA, self.uniqB = self._uniq()
            self.innerA, self.innerB = self._innermost()

    def _mapq(self):
        mapqA = 0
        mapqB = 0
        for read in self.reads:
            A, B = [int(x) for x in read.split()[7:9]]
            mapqA += A
            mapqB += B

        return mapqA / self.size, mapqB / self.size

    def _uniq(self):
        """
        Compute weighted mean of 'uniqueness'

        Uniqueness is currently measured as the number of unique mapping
        positions within the cluster.
        """
        wmeanA = 0
        wmeanB = 0
        for sample in self.samples.keys():
            posA = []
            posB = []
            for read in self.samples[sample]:
                A, B = [int(x) for x in read.split()[5:7]]
                posA.append(A)
                posB.append(B)
            uniqA = len(set(posA)) / len(posA)
            uniqB = len(set(posB)) / len(posB)

            # Weights (fraction of total) sum to 1, so can simply add
            wmeanA += uniqA * (len(posA) / self.size)
            wmeanB += uniqB * (len(posB) / self.size)

        return wmeanA, wmeanB

    def global_cov(self, gcov_dict):
        """
        Compute global cov of cluster.
        Current computed as maximum global cov across samples

        Args:
            gcov_dict (dict) : {Sample_ID : global}
        """

        gcov = 0
        for sample in self.samples.keys():
            sample_gcov = len(self.samples[sample]) / gcov_dict[sample]
            if sample_gcov > gcov:
                gcov = sample_gcov

        return gcov

    def local_cov(self, insert_handles):
        """
        Args:
            insert_handles (dict) : Insert coverage files
                                    {sample_ID (str): icov_f (pysam.Tabixfile)}
        """

        localA = 0
        localB = 0

        # import pdb
        # pdb.set_trace()
        for sample in self.samples.keys():
            sample_size = len(self.samples[sample])

            posAs = [int(read.split()[5]) for read in self.samples[sample]]
            posBs = [int(read.split()[6]) for read in self.samples[sample]]
            icov_f = insert_handles[sample]

            innerA = max(posAs) if self.strandA == '+' else min(posAs)
            innerB = max(posBs) if self.strandB == '+' else min(posBs)

            try:
                regionA = icov_f.fetch(self.chrA, innerA, innerA + 1).next()
            except StopIteration:
                pdb.set_trace()
            smpl_localA = sample_size / (int(regionA.split()[3]) + sample_size)
            if smpl_localA > localA:
                localA = smpl_localA

            try:
                regionB = icov_f.fetch(self.chrB, innerB, innerB + 1).next()
            except StopIteration:
                pdb.set_trace()
            smpl_localB = sample_size / (int(regionB.split()[3]) + sample_size)
            if smpl_localB > localB:
                localB = smpl_localB

        return localA, localB

    def classify(self, gcov, localA, localB):
        # Constants for thresholds
        UNIQ_B1_T = 0.9
        GLOBAL_B1_T = 2.5
        MAPQ_B2_T = 25
        LOCAL_B3_T = 0.25
        GLOBAL_B4_T = 0.1
        LOCAL_B5_T = 0.25
        MAPQ_B5_T = 10

        # Branch 1 - Uniqueness > 90%
        if (self.uniqA >= UNIQ_B1_T and self.uniqB >= UNIQ_B1_T and
            gcov < GLOBAL_B1_T):

            # Branch 2 - MapQ > 25
            if self.qualA >= MAPQ_B2_T and self.qualB >= MAPQ_B2_T:
                # Branch 3 - Local coverage ratio >= 0.25
                if localA >= LOCAL_B3_T and localB >= LOCAL_B3_T:
                    return "Valid", 3
                else:
                    # Branch 4 - Global coverage ratio >= 0.1
                    if gcov >= GLOBAL_B4_T:
                        return "Valid", 4
                    else:
                        return "Invalid", 4
            else:
                # Branch 5 - Local coverage ratio >= 0.25, mapQ >= 10
                if (localA >= LOCAL_B5_T and localB >= LOCAL_B5_T and
                    self.qualA >= MAPQ_B5_T and self.qualB >= MAPQ_B5_T):

                    return "Valid", 5
                else:
                    return "Invalid", 5
        else:
            return "Invalid", 1

    def _innermost(self):
        """
        Return innermost positions

        Returns:
            (int, int) : innerA, innerB
        """

        innerA = self.endA if self.strandA == '+' else self.startA
        innerB = self.endB if self.strandB == '+' else self.startB

        return innerA, innerB

    @classmethod
    def init_rpc(cls, reads, name, sample_delim="___"):
        """
        Create a Cluster object from reads output by ReadPairCluster

        REQUIRES read IDs to be prepended with ID in format COHORT_SAMPLE
        """

        chrA = reads[0].split()[3]
        startA = sys.maxsize
        endA = 0

        chrB = reads[0].split()[4]
        startB = sys.maxsize
        endB = 0

        strandA = strand(int(reads[0].split()[11]))
        strandB = strand(int(reads[0].split()[12]))

        qualAs = []
        qualBs = []
        posAs = []
        posBs = []
        samples = {}

        for read in reads:
            [posA, posB] = [int(x) for x in read.split()[5:7]]

            posAs.append(posA)
            posBs.append(posB)

            if posA < startA:
                startA = posA
            if posA > endA:
                endA = posA
            if posB < startB:
                startB = posB
            if posB > endB:
                endB = posB

            [qualA, qualB] = [int(x) for x in read.split()[7:9]]

            qualAs.append(qualA)
            qualBs.append(qualB)

            # Allow sample level access to reads
            # Instead of assuming STUDY_SAMPLE ID format, 
            # require distinct delimiter
            # sample = "_".join(read.split()[2].split('_')[0:2])
            sample = read.split()[2].split(sample_delim)[0]
            if sample in samples.keys():
                samples[sample].append(read)
            else:
                samples[sample] = [read]

        size = len(reads)
        uniqA = len(set(posAs)) / size
        uniqB = len(set(posBs)) / size
        qualA = np.mean(qualAs)
        qualB = np.mean(qualBs)
        score = '.'

        return Cluster([chrA, startA, endA, chrB, startB, endB, name, score,
                        strandA, strandB, size, qualA, qualB, uniqA, uniqB],
                       reads, samples)

    def str_with_samples(self, samples):
        """
        Returns s
        """

        counts = {}
        for sample in samples:
            counts[sample] = 0

        for read in self.reads:
            name = read.split()[2]
            counts[name.split('__')[0]] += 1

        num_samples = len([x for x in samples if counts[x] > 0])

        probands = [x for x in samples if x.endswith('p1')]
        num_probands = len([x for x in probands if counts[x] > 0])

        count_str = '\t'.join(str(counts[sample]) for sample in samples)

        return '\t'.join([str(self), str(num_samples), str(num_probands), count_str])

    def __repr__(self):
        return ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t.\t%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t"
                "%.3f" %
                (self.chrA, self.startA, self.endA,
                 self.chrB, self.startB, self.endB,
                 self.name, self.strandA, self.strandB,
                 self.size, self.qualA, self.qualB, self.uniqA, self.uniqB))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)


def merge(clusters):
    if len(clusters) == 1:
        return clusters[0]

    reads = []
    for cluster in clusters:
        reads.extend(cluster.reads)

    # Any clusters being merged follow the .01, .02 naming convention
    name = clusters[0].name[:-3]

    merged = Cluster.init_rpc(reads, name)

    # Merged clusters have conflicting orientations
    merged.strandA = '.'
    merged.strandB = '.'

    return merged


def strand(flag):
    reverse = 16
    if flag & reverse:
        return '-'
    else:
        return '+'
