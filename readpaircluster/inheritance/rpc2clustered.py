#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse
import sys


def rpc2clustered(bedpe):

    header = next(bedpe).strip().split()
    samples = header[17:]

    for line in bedpe:
        data = line.strip().split()
        chrom = data[0]
        start = data[1] if data[8] == '-' else data[2]
        end = data[4] if data[9] == '-' else data[5]
        name = data[6]
        score = data[7]
        strandA = data[8]
        strandB = data[9]
        size = int(end) - int(start)
        counts = [int(x) for x in data[17:]]
        num_samples = len([x for x in counts if x > 0])
        sample_names = ','.join([
            sample for sample, count in zip(samples, counts) if count > 0])
        num_progs = 1
        progs = 'rpc'
        uniqueIDs = name

        record = '\t'.join([str(x) for x in [
            chrom, start, end, name, score, strandA, strandB, size,
            num_samples, sample_names, num_progs, progs, name]])

        yield record

    # bedpe = pd.read_table(bedpe, index_col=False)
    # samples = list(bedpe.columns.values)[15:]
    # c = pd.DataFrame()
    # c['chr'] = bedpe['#chrA']
    # c['start'] = bedpe.loc[bedpe['strandA'] == '+', 'endA']
    # c['start'].fillna(bedpe.loc[bedpe['strandA'] == '-', 'startA'], inplace=True)
    # c['end'] = bedpe.loc[bedpe['strandB'] == '+', 'endB']
    # c['end'].fillna(bedpe.loc[bedpe['strandB'] == '-', 'startB'], inplace=True)
    # c['name'] = bedpe['name']
    # c['score'] = '.'
    # c['strandA'] = bedpe['strandA']
    # c['strandB'] = bedpe['strandB']
    # c['size'] = c['end'] - c['start']
    # c['num_samples'] = bedpe.loc[:, samples[0]:samples[-1]] \
        # .applymap(lambda x: 1 if x > 0 else 0) \
        # .sum(axis=1)

    # c['num_progs'] = 1
    # c['progs'] = 'rpc'
    # c['uniqueIDs'] = bedpe['name']

def main():
    parser = argparse.ArgumentParser(
        description="")
    parser.add_argument('bedpe', type=argparse.FileType('r'))
    parser.add_argument('fout', type=argparse.FileType('w'),
                        nargs='?', default=sys.stdout)

    args = parser.parse_args()

    for record in rpc2clustered(args.bedpe):
        args.fout.write(record + '\n')


if __name__ == '__main__':
    main()
