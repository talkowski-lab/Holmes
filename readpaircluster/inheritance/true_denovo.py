#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2015 msto <mstone5@mgh.harvard.edu>
#
# Distributed under terms of the MIT license.

"""

"""

import sys
import pandas as pd
import matplotlib.pyplot as plt


def true_denovo(caller, svtype):
    rates = pd.DataFrame({})
    dc = pd.DataFrame({})
    denovo_vars = {}
    for freq in [0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
        svfile = 'all.{0}.{1}.{2}.sv.txt'.format(caller, svtype, freq)
        sv = pd.read_table(svfile, index_col=0)
        sv = sv.astype(bool).astype(int)  # convert > 0 to 1
    
        # proband_totals = sv.filter(regex='p1').sum()
        variant_totals = sv.sum(axis=1)
    
        sample_counts = sv.groupby(lambda s: s[-2:], axis=1).sum()
        denovo = sv[sample_counts['p1'] == 1]  # Appear in only one proband
        denovo = denovo.loc[variant_totals == 1]  # and only that proband
        denovo_vars[freq] = denovo

        denovo_counts = denovo.sum().filter(regex='p1')
        denovo_counts = denovo_counts.rename(index=lambda x: x[:-2])
        dc[freq] = denovo_counts
    
        cfile = 'all.{0}.{1}.{2}.counts.txt'.format(caller, svtype, freq)
        counts = pd.read_table(cfile, index_col=0)
        rates[freq] = denovo_counts / counts[['biparental', 'denovo', 'uniparental']].sum(axis=1)

    return dc, rates, denovo_vars

if __name__ == '__main__':

    caller = sys.argv[1]
    svtype = sys.argv[2]

    dc, rates, dv = true_denovo(caller, svtype)

    flierprops = dict(marker='o', markerfacecolor='green', markersize=7, linestyle='none')
    rates.plot(kind='box', ylim=(0,1), title='{0} {1}'.format(caller, svtype), flierprops=flierprops)
    plt.savefig('{0}.{1}.true_denovo_rates.png'.format(caller, svtype))
    
    dc.plot(kind='box', title='{0} {1}'.format(caller, svtype), flierprops=flierprops)
    plt.savefig('{0}.{1}.true_denovo_counts.png'.format(caller, svtype))
