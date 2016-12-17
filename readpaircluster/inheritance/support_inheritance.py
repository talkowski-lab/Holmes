#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""

"""

import argparse
import pandas as pd
# import numpy as np


def read_sv(bed, size=500, freq=1.0, samples=None, bedpe=True):
    """
    Parameters
    ----------
    bed : str or handle
        a) Bedpe file in compress_clusters.py --samples format
        Column names: chrA, startA, endA, chrB, startB, endB, score,
        strandA, strandB, size, qualA, qualB, uniqA, uniqB, num_samples,
        num_probands, [sample for sample in samples]
        b) Bed file in deepclustering.sh format
        Column names: chr, start, end, svID, score, strand1, strand2, size,
        num_samples, num_probands, samples, num_progs, progs, uniqueIDs
    size : int
        Minimum size of variants to allow [500]
        (Calculated as median(startB, endB) - median(startA, endA) for bedpe)
    freq : float
        Maximum allele frequency (within probands) to allow [0.9]
    samples : list
        List of sample names
        If None, parsed from header (bedpe) or samples column (deepcluster)
    bedpe : bool
        bed in bedpe format [True]

    Returns
    -------
    bed: pd.DataFrame
        Index: svID
        Column names: chr, start, end

    sv : pd.DataFrame
        Index: svID
        Column names: [sample for sample in samples]
    """
    if bedpe:
        sv = pd.read_table(bed, index_col=6)
        samples = list(sv.columns.values[16:])
        sv['start'] = sv[['startA', 'endA']].median(axis=1).astype(int)
        sv['end'] = sv[['startB', 'endB']].median(axis=1).astype(int)
    else:
        colnames = ['#chrA', 'start', 'end', 'svID', 'score', 'strand1',
                    'strand2', 'size', 'num_samples', 'samples',
                    'num_progs', 'progs', 'uniqueIDs']

        sv = pd.read_table(bed, index_col=3,
                           header=None, names=colnames)

        # Extract bool vectors
        for sample in samples:
            sv[sample] = pd.Series(sv['samples'].str.contains(sample))

        if samples is None:
            samples = sorted(list(set(','.join(sv['samples']).split(','))))

    # import pdb
    # pdb.set_trace()
    # TODO: make robust against future sample IDs
    # quads = sorted(list(set([x[:-2] for x in samples])))

    # Size filter
    sv = sv[abs(sv['end'] - sv['start']) >= size].copy()

    # Compute allele frequency
    # frequency = compute_freq(sv)

    # Exclude variants that appears in >90% of parents
    # max_count = int(freq * len(quads) * 2)
    # sv = sv[frequency['parent'] <= max_count].copy()

    sv = sv.loc[:, samples[0]:samples[-1]].copy()

    # return sv, frequency
    return sv


def bool_sv(sv, null='N', keep=''):
    # Convert table of 'N', 'PASS', etc to boolean table
    supports = pd.Series(sv.values.ravel()).unique()
    replacements = {s: 1 if (s != null and keep in s) else 0 for s in supports}
    return sv.replace(replacements)


def compute_freq(sv):
    # def member(sample):
        # if 'SFARI_d' not in sample:
            # return 'not_sample'
        # else:
            # return sample[-2:]

    member = lambda x: x[-2:]

    freq = sv.groupby(member, axis=1).sum()
    # freq = freq.drop('not_sample', axis=1)

    freq['total'] = freq.sum(axis=1)
    freq['parent'] = freq['mo'] + freq['fa']

    return freq


def freq_filter(sv, freq=0.9):
    # frequency = compute_freq(sv)
    # frequency['parent'] = frequency['fa'] + frequency['mo']

    # Exclude variants that appears in >90% of parents
    max_count = int(freq * len(quads) * 2)

    return sv[frequency['parent'] <= max_count].copy()

def status(quad, member='p1'):
    q = quad.index.values[0][:-2]
    p1 = quad[q + 'p1']
    s1 = quad[q + 's1']
    fa = quad[q + 'fa']
    mo = quad[q + 'mo']

    if member == 'p1':
        target = p1
        alt = s1
    elif member == 's1':
        target = s1
        alt = p1
    else:
        raise Exception('Invalid member')

    if target & ~(alt | mo | fa):
        return 'denovo'
    if target & alt & ~(mo | fa):
        return 'false_denovo'
    if target & (mo ^ fa):
        return 'uniparental'
    if target & (mo & fa):
        return 'biparental'
    return 'absent'

def compute_inheritance(sv):
    """
    Parameters
    ----------
    sv : pd.DataFrame
        Index: svID
        Column names: [sample for sample in samples]

    Returns
    -------
    inherit : pd.DataFrame
        Matrix of (absent, denovo, false_denovo, uniparental, biparental)
        Index: svID
        Column names: [proband for proband in probands]
    counts : pd.DataFrame
        Index: quad names (corresponding to probands)
        Column names: absent, denovo, uniparental, biparental
    """

    def status(quad, member='p1'):
        q = quad.index.values[0][:-2]
        p1 = quad[q + 'p1']
        s1 = quad[q + 's1']
        fa = quad[q + 'fa']
        mo = quad[q + 'mo']

        if member == 'p1':
            target = p1
            alt = s1
        elif member == 's1':
            target = s1
            alt = p1
        else:
            raise Exception('Invalid member')

        if target & ~(alt | mo | fa):
            return 'denovo'
        if target & alt & ~(mo | fa):
            return 'false_denovo'
        if target & (mo ^ fa):
            return 'uniparental'
        if target & (mo & fa):
            return 'biparental'
        return 'absent'

    # Group by quadname and compute inheritance status
    family = lambda x: x[:-2]
    inherit = sv.groupby(family, axis=1) \
        .apply(lambda g: g.apply(status, axis=1))

    melted = pd.melt(inherit, var_name='sample', value_name='inherit')
    counts = pd.crosstab(index=[melted['sample']], columns=[melted['inherit']])

    return inherit, counts


def plot_inheritance(svfile, size=500, freqs=[0.1, 0.25, 0.5, 0.75, 0.9]):
    rates = pd.DataFrame({})
    for freq in freqs:
        bed, sv = read_sv(svfile, size=size, freq=freq)
        inherit, counts = compute_inheritance(sv)
        rates[freq] = counts['denovo'] / counts[['biparental', 'denovo', 'uniparental']].sum(axis=1)
    return rates


def main():
    parser = argparse.ArgumentParser(
        description="")

    args = parser.parse_args()

    compute_inheritance()


if __name__ == '__main__':
    main()
