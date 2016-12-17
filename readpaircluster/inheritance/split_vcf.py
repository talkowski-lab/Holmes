#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import argparse


def split_vcf(vcf, fouts):
    for line in vcf:
        if line.startswith('#'):
            continue
        
        data = line.strip().split()
        chrA = data[0]
        INFO = {k: v for k, v in [f.split('=') for f in data[7].split(';') if '=' in f]}
        svtype = INFO['SVTYPE']
        alt = data[4]

        if chrA.startswith('GL') or chrA.startswith('MT'):
            continue

        if svtype not in ['BND']:
            fouts[svtype].write(line)
        else:
            chrB = alt.strip('[]N').split(':')[0]
            if chrB.startswith('GL') or chrB.startswith('MT'):
                continue
            if chrA == chrB:
                fouts['INV'].write(line)
            else:
                fouts['TRA'].write(line)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument('vcf', type=argparse.FileType('r'))
    parser.add_argument('sample')
    parser.add_argument('prog')

    args = parser.parse_args()

    fouts = {
        'DEL': open('{0}.{1}.{2}.vcf'.format(args.sample, args.prog, 'del'), 'w'),
        'DUP': open('{0}.{1}.{2}.vcf'.format(args.sample, args.prog, 'dup'), 'w'),
        'INV': open('{0}.{1}.{2}.vcf'.format(args.sample, args.prog, 'inv'), 'w'),
        'TRA': open('{0}.{1}.{2}.vcf'.format(args.sample, args.prog, 'tloc'), 'w')
    }

    split_vcf(args.vcf, fouts)


if __name__ == '__main__':
    main()
