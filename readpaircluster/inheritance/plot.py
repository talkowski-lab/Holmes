#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 msto <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""


"""
# from sys import maxsize
import matplotlib.pyplot as plt
# from inheritance import plot_inheritance
import pandas as pd


flierprops = dict(marker='o', markerfacecolor='green', markersize=7,
                  linestyle='none')

# Plot counts 
# for sv in 'del dup inv tloc'.split():
    # for caller in 'lumpy delly rpc'.split():
        # calls = pd.DataFrame({})

        # for freq in [0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
            # counts = pd.read_table(
                # 'all.{0}.{1}.{2}.counts.txt'.format(caller, sv, freq),
                # index_col=0)
            # try:
                # calls[freq] = counts[['denovo', 'uniparental', 'biparental']].sum(axis=1)
            # except:
                # import pdb
                # pdb.set_trace()
        
        # if caller != 'rpc':
            # bp500 = pd.read_table('raw_counts/probands.raw_counts.500bp.{0}.{1}.txt'.format(caller, sv), header=None, names=('sample', 'count'), index_col=0)
            # calls['500bp'] = bp500['count']

            # raw = pd.read_table('raw_counts/probands.raw_counts.{0}.{1}.txt'.format(caller, sv), header=None, names=('sample', 'count'), index_col=0)
            # calls['raw'] = raw['count']

        # calls.plot(kind='box', title='{0} {1}'.format(caller, sv),
                   # flierprops=flierprops)
        # plt.savefig('{0}.{1}.box_counts.full.png'.format(caller, sv))

for sv in 'del dup inv tloc'.split():
    for caller in 'lumpy delly rpc'.split():
        rates = pd.DataFrame({})
        for freq in [0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
            counts = pd.read_table(
                'all.{0}.{1}.{2}.counts.txt'.format(caller, sv, freq),
                index_col=0)

            counts.drop('absent', axis=1, inplace=True)
            rates[freq] = counts['denovo'] / counts.sum(axis=1)

        rates.plot(kind='box', ylim=(0, 1), title='{0} {1}'.format(caller, sv),
                   flierprops=flierprops)
        plt.savefig('{0}.{1}.family_rates.png'.format(caller, sv))
