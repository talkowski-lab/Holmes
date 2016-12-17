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
import os
import argparse
from inheritance import read_sv, compute_inheritance


def save_inheritance(bedpe, svtype=None, prefix=None):
    # terrible hack don't do this
    if svtype is None:
        if 'tloc' in bedpe:
            svtype = 'tloc'
        else:
            svtype = 'other'

    if svtype == 'tloc':
        sv, freq = read_sv(bedpe, size=(0 - sys.maxsize), freq=1)
    else:
        sv, freq = read_sv(bedpe, freq=1)

    inherit, counts = compute_inheritance(sv)

    if not prefix:
        prefix = os.path.splitext(bedpe)[0]
    sv.to_csv('{0}.sv.csv'.format(prefix))
    freq.to_csv('{0}.freq.csv'.format(prefix))
    inherit.to_csv('{0}.inherit.csv'.format(prefix))
    counts.to_csv('{0}.counts.csv'.format(prefix))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bedpe')
    parser.add_argument('-s', '--svtype', default=None)
    args = parser.parse_args()

    save_inheritance(args.bedpe, args.svtype)

if __name__ == '__main__':
    main()
