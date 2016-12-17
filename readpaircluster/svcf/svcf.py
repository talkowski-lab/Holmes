#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""


from collections import OrderedDict


# TODO: define dict of FORMAT fields and functions to merge them across Records

def _order_dict(d):
    """Convert dict to sorted OrderedDict"""
    if isinstance(d, dict):
        d = [(k, v) for k, v in d.items()]
        d = sorted(d, key=lambda x: x[0])
        return OrderedDict(d)
    elif isinstance(d, OrderedDict):
        return d
    else:
        raise Exception('Must be dict or OrderedDict')


class SVCFRecord(object):
    """
    An SVCF Record.

    Attributes
    ----------
    chrA: str
        Chromosome A (first breakpoint)
    posA: int
        Position A
    chrB: str
        Chromosome B (second breakpoint)
    posB: int
        Position B
    name: str
        Record identifier
    svtype: str
        SV type of the record (currently limited to del, dup, inv, tloc)
    SVINFO: dict
        Variant metrics and annotations
    samples: dict
        Sample-specific call info
    SVFORMAT: list

    FORMAT: list or str
        List or colon-delimited string of sample-specific fields

    """

    def __init__(self, chrA, posA, chrB, posB, svtype, SVINFO, samples,
                 SVFORMAT=None, FORMAT=None):
        self.chrA = chrA
        self.posA = posA
        self.chrB = chrB
        self.posB = posB
        self.svtype = svtype

        # TODO: check performance impact sorting here
        if SVFORMAT is None:
            self.SVINFO = _order_dict(SVINFO)
            self.SVFORMAT = list(self.SVINFO.keys())
        else:
            self.SVINFO = OrderedDict([(f, SVINFO[f]) for f in SVFORMAT])
            self.SVFORMAT = SVFORMAT

        # TODO: check for or fill blank fields
        self.samples = _order_dict(samples)
        if FORMAT is None:
            for s in self.samples:
                self.samples[s] = _order_dict(self.samples[s])
            self.FORMAT = list(self.samples[s].keys())
        else:
            if isinstance(FORMAT, str):
                FORMAT = FORMAT.split(':')
            for s in self.samples:
                ordered_kv = [(f, self.samples[s][f]) for f in FORMAT]
                self.samples[s] = OrderedDict(ordered_kv)
            self.FORMAT = FORMAT

    def __repr__(self):
        base = '\t'.join([self.chrA, self.posA, self.posB, self.chrB,
                          self.svtype])
        SVFORMAT = ':'.join(self.SVFORMAT)
        SVINFO = ':'.join(self.SVINFO.values())

        FORMAT = ':'.join(self.FORMAT)
        samples = [':'.join(v) for v in self.samples.values()]
        samples = '\t'.join(samples)

        return '\t'.join([base, SVFORMAT, SVINFO, FORMAT, samples])
