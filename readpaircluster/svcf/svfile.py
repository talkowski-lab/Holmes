#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

import os
from io import open, IOBase
from collections import OrderedDict
import numpy as np
import vcf
from pybedtools import BedTool
from .svcall import DellyCall, LumpyCall, RPCCall, BedCall, BedpeCall, VCFCall, PindelCall


class UnsupportedFiletypeError(Exception):
    """Unsupported or unexpected filetypes"""


class ImproperFormattingError(Exception):
    """Wrong number of fields or other formatting error"""


def BedpeParser(filename):
    if not filename.endswith('bedpe'):
        filetype = os.path.splitext(filename)[1]
        raise UnsupportedFiletypeError('%s (expected bedpe)' % filetype)

    bedpe = open(filename)
    header = next(bedpe)
    if not header.startswith('#'):
        header = None
        bedpe.seek(0)
    else:
        header = header.strip('#')

    for line in bedpe:
        data = line.strip().split()
        if len(data) < 10:
            msg = 'Too few fields (%d) (bedpe requires 10)' % len(data)
            raise ImproperFormattingError(msg)
        if header:
            INFO = OrderedDict(zip(header.split()[10:], data[10:]))
        else:
            INFO = {'INFO': data[10:]}

        yield BedpeRecord(*data[:10], INFO=INFO)


class BedpeRecord(object):
    def __init__(self, chrA, startA, endA, chrB, startB, endB, name, svtype,
                 strandA, strandB, INFO={}):
        self.chrA = chrA
        self.startA = int(startA)
        self.endA = int(endA)
        self.chrB = chrB
        self.startB = int(startB)
        self.endB = int(endB)
        self.name = name
        self.svtype = svtype
        self.strandA = strandA
        self.strandB = strandB
        self.INFO = INFO

    @property
    def posA(self):
        if self.strandA == '.':
            return int(np.median([self.startA, self.endA]))
        else:
            return self.startA if self.strandA == '-' else self.endA

    @property
    def posB(self):
        if self.strandB == '.':
            return int(np.median([self.startB, self.endB]))
        else:
            return self.startB if self.strandB == '-' else self.endB

    def __str__(self):
        if 'INFO' in self.INFO.keys():
            INFO = '\t'.join([str(x) for x in self.INFO['INFO']])
        else:
            INFO = '\t'.join([str(v) for k, v in self.INFO])

        vals = [str(x) for x in [self.chrA, self.startA, self.endA,
                                 self.chrB, self.startB, self.endB,
                                 self.name, self.svtype,
                                 self.strandA, self.strandB, INFO]]
        return '\t'.join(vals)


class SVWriter(object):
    def __init__(self, *svfiles):
        # Set: infos, metadata, formats, filters, alts, contigs, _column_headers, samples
        pass

class SVFile(object):
    def __init__(self, filename, filetype=None,
                 sample=None, source=None, svtype=None):

        #TODO: allow passing of open file handle for piping

        gzipped = filename.endswith('.gz')

        self.sample = sample
        self.source = source

        # TODO: add bstat, rpc?
        filetypes = 'vcf bed bedpe bstat gz'.split()
        if filetype is None:
            def _get_suffix(fname):
                return os.path.splitext(fname)[-1].strip('.')

            if gzipped:
                filetype = _get_suffix(filename[:-3])
            else:
                filetype = _get_suffix(filename)

            if filetype not in filetypes:
                raise UnsupportedFiletypeError(filetype)

        self.filetype = filetype

        if filetype == 'vcf':
            reader = vcf.Reader(filename=filename, compressed=gzipped)
        elif filetype == 'bedpe':
            reader = BedpeParser(filename)
        elif filetype == 'bed':
            reader = BedTool(filename)
        elif filetype == 'bstat':
            raise UnsupportedFiletypeError('bstat not yet supported')

        if source == 'delly':
            parser = (DellyCall(record) for record in reader)
        elif source == 'lumpy':
            parser = (LumpyCall(record) for record in reader)
        elif source == 'pindel':
            parser = (PindelCall(record) for record in reader)
        elif source == 'rpc':
            parser = (RPCCall(record) for record in reader)
        else:
            if filetype == 'bed':
                parser = (BedCall(record) for record in reader)
            elif filetype == 'bedpe':
                parser = (BedpeCall(record, sample=sample, source=source) for record in reader)
            elif filetype == 'vcf':
                parser = (VCFCall(record, sample=sample, source=source) for record in reader)
            else:
                raise UnsupportedFiletypeError('%s not yet supported' % source)
        # elif source is None:
            # if filetype == 'vcf':
                # try:
                    # parser = (BedpeCall(record) for record in reader)
                # except:
                    # raise UnknownFormatError(

        # TODO: use filter methods not Nonetype comparison
        self.parser = (call for call in parser if call is not None)

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        return next(self.parser)


class Reader(object):
    def __init__(self):
        pass


class Writer(object):
    def __init__(self):
        #  header = '\t'.join('chrA posA posB chrB name'.split())
        #  header = header + '\t'.join(samples)
        pass
