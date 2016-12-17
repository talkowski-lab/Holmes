#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

from collections import defaultdict, OrderedDict
import numpy as np
from .rpc import RPCNode
from .svcf import SVCFRecord


class JointVCFError(Exception):
    """Multi-sample VCFs not yet supported"""


class UnparsableSVTypeError(Exception):
    """Unable to translate caller-specific SV type to del, dup, inv, tloc"""


class NonstandardVCFError(Exception):
    """VCF not formatted according to 4.2 spec and no template defined"""


SOURCE_FIELDS = OrderedDict([
    ('delly', ['PE_delly',
               'SR_delly',
               'MAPQ_delly',
               'DR_delly',
               'DV_delly',
               'RR_delly',
               'RV_delly']),
    ('lumpy', ['PE_lumpy',
               'SR_lumpy',]),
    ('pindel', ['RR_pindel',
                'RV_pindel']),
    ('rpc', ['PE_rpc',
             'QualA_rpc',
             'QualB_rpc',
             'UniqA_rpc',
             'UniqB_rpc'])])


class SVCallCluster(object):
    def __init__(self, samples, sources, calls, name='.', merge_fn=np.median):
        """
        Parameters
        ----------
        samples: list of str
            List of sample IDs
        sources: list of str
            List of source program names
        calls: list of SVCall
            List of SVCalls represented by the Record
        name: str, optional
            Record identifier (default: '.')
        merge_fn: function, optional
            Function to obtain Record start/end from constituent calls.
            Currently takes a list of numeric types (default: np.median),
            taking min/max based on call strands a todo

        Attributes
        ----------
        samples: dict
        sources: dict of str, int
            Dict of source names, counts of SVCalls from source
            source
        name: str
            Record identifier
        chrA: str
            Chromosome A (first breakpoint)
        posA: int
            Position A
        chrB: str
            Chromosome B (second breakpoint)
        posB: int
            Position B
        svtype: str
            SV type of the record (currently limited to del, dup, inv, tloc)
        """

        self.calls = calls
        self.name = name
        self.chrA = calls[0].chrA
        self.posA = int(merge_fn([call.posA for call in calls]))
        self.chrB = calls[0].chrB
        self.posB = int(merge_fn([call.posB for call in calls]))
        self.svtype = calls[0].svtype

        self._sample_list = sorted(samples)
        self._source_list = sorted(sources)

        self.samples = defaultdict(list)
        for call in calls:
            self.samples[call.sample].append(call)

        #  self.samples = sorted(set([c.sample for c in calls]))

        #  self.samples = defaultdict(list)
        #  self.sources = defaultdict(int)
        #  self._sample_sources = defaultdict(list)

        #  for call in calls:
            #  self.samples[call.sample].append(call)
            #  # count number of samples with source support, not number of
            #  # calls from source
            #  if (call.sample, call.source) not in self._sample_sources:
                #  self.sources[call.source] += 1
            #  self._sample_sources[(call.sample, call.source)].append(call)

        # Set basic attributes

    @property
    def FORMAT(self):
        """Record FORMAT field"""
        #  sources = ':'.join(self._source_list)
        #  fields = ':'.join([key for source in self._source_list
                           #  for key in SOURCE_FIELDS[source]])

        fields = [key for source in self._source_list
                      for key in SOURCE_FIELDS[source]]

        return self._source_list + fields

    @property
    def freq(self):
        """Number of samples in the Record"""
        return len([s for s in self.samples if len(self.samples[s]) > 0])

    @property
    def svcf(self):
        """Return corresponding SVCF Record"""

        SVINFO = {'freq': str(self.freq)}

        def _compile_genotypes(calls):
            sources = defaultdict(list)
            for call in calls:
                sources[call.source].append(call)

            gt = OrderedDict()
            for source in self._source_list:
                if source in sources:
                    if len(sources[source]) == 1:
                        sgt = sources[source].svcf
                    else:
                        scalls = sources[source]
                        sggt = scalls[0].joint_svcf(scalls[1:])

                    for field in SOURCE_FIELDS[source]:
                        gt[field] = sgt[field]
                else:
                    for field in SOURCE_FIELDS[source]:
                        gt[field] = 0

            return gt

        FORMAT = [f for s in self._source_list for f in SOURCE_FIELDS[s]]
        null_gt = OrderedDict([(f, 0) for f in FORMAT])

        genotypes = OrderedDict()
        for sample in self._sample_list:
            if sample in self.samples:
                genotypes[sample] = _compile_genotypes(self.samples[sample])
            else:
                genotypes[sample] = null_gt

        return SVCFRecord(self.chrA, self.posA, self.chrB, self.posB,
                          self.svtype, SVINFO, genotypes, FORMAT=FORMAT)

    @property
    def links(self):
        """Link joint SVCF ID to constituent caller IDs"""

        links = [(self.name, c.svtype, c.sample, c.source, c.name,
                  c.caller_svtype) for c in self.calls]

        links = ['\t'.join(link) for link in links]
        links = '\n'.join(links)
        return links

    @property
    def inherit(self):
        #  self.samples = OrderedDict((sample, []) for sample in samples)
        #  self.sources = OrderedDict((source, 0) for source in sources)
        #  self._sample_sources = defaultdict(list)
        #  for call in calls:
            #  self.samples[call.sample].append(call)
            #  self.sources[call.source] += 1
            #  self._sample_sources[(call.sample, call.source)].append(call)

        probands = [s for s in self.samples if s.endswith('p1')]
        probands = sorted(probands)
        members = 'fa mo s1'.split()

        # Consider inheritance independent of source
        # TODO: allow inheritance only if sources concordant
        inherit = [[m for m in members if (p[:-2] + m) in self.samples] for p in probands]
        inherit = [','.join(status) for status in inherit]
        inherit = [i if len(i) > 0 else 'denovo' for i in inherit]

        source_lists = [[c.source for c in self.samples[p]] for p in probands]
        source_lists = [','.join(sorted(set(s))) for s in source_lists]

        freq = str(self.freq)
        stats = zip(probands, source_lists, inherit)

        statuses = ['\t'.join([str(self), freq, p, s, i]) for p, s, i in stats]

        return '\n'.join(statuses)

    def __str__(self):
        return ('{chrA}\t{posA}\t{posB}\t{chrB}\t{name}\t{svtype}'.format(
                **self.__dict__))


class SVCall(RPCNode):
    def __init__(self, chrA, posA, chrB, posB,
                 name='.', sample=None, source=None, svtype=None):
        """
        Common format for SV calls for intersection analyses.

        Includes methods for RPC-based clustering
        """
        super().__init__(chrA, posA, chrB, posB, name)

        self.name = name
        self.sample = sample
        self.source = source
        self.svtype = svtype

    def clusters_with(self, other, dist):
        return (super().clusters_with(other, dist) and
                self.svtype == other.svtype)

    @property
    def HQ(self):
        return True

    @property
    def secondary(self):
        return False

    @property
    def caller_svtype(self):
        return self.svtype

    @property
    def svsize(self):
        if self.chrA == self.chrB:
            return self.posB - self.posA
        else:
            return -1

    @property
    def svcf(self):
        return {}

    def joint_svcf(self, others):
        s1 = self.svcf

        for other in others:
            s2 = other.svcf
            for key in s1:
                s1[key] += s2[key]

        return s1

    # TODO: include samples in equality test?
    def __eq__(self, other):
        return (self.chrA == other.chrA and
                self.posA == other.posA and
                self.chrB == other.chrB and
                self.posB == other.posB and
                self.name == other.name)

    def __str__(self):
        return ('{chrA}\t{posA}\t{posB}\t{chrB}\t{callname}'.format(
                chrA=self.chrA, posA=self.posA, posB=self.posB, chrB=self.posB,
                callname=self.name))

STD_SVTYPES = {
    'DEL': 'del',
    'DUP': 'dup',
    'DUP:TANDEM': 'dup',
    'INV': 'inv',
    'TRA': 'tloc',
    'INS': 'dup',
    'RPL': 'del'}

class VCFCall(SVCall):
    def __init__(self, record, chrB=None, posB=None,
                 sample=None, source=None, svtype=None):
        """
        Wrapper for VCF-formatted callers.

        Computation of chromosome and position left to subclasses
        """
        if len(record.samples) > 1:
            raise JointVCFError('Multi-sample VCFs not yet supported')

        self.record = record
        chrA = str(record.CHROM)
        posA = record.POS
        name = record.ID
        if name is None:
            name = '.'
        if sample is None:
            sample = self.record.samples[0].sample

        # Detect svtype
        if svtype is None:
            svtype = record.INFO['SVTYPE']
            if svtype in STD_SVTYPES:
                svtype = STD_SVTYPES[svtype]
            else:
                allowed = ','.join(STD_SVTYPES)
                msg = 'SV type %s not in known types: %s\n' % (svtype, allowed)
                msg = msg + 'Sample: %s\tSource: %s\tRecord: %s' % (sample, source, name)
                raise UnparsableSVTypeError(msg)

        # Parse chrB based on VCF standard
        if chrB is None:
            if svtype == 'tloc':
                if hasattr(record.ALT[0], 'chr'):
                    chrB = str(record.ALT[0].str)
                    posB = record.ALT[0].pos
                else:
                    msg = 'Could not parse chrB from translocation record\n'
                    msg = msg + 'Sample: %s\tSource: %s\tRecord: %s' % (sample, source, name)
                    raise NonstandardVCFError(msg)
            else:
                chrB = chrA

        # Parse posB based on VCF standard
        if posB is None:
            if 'END' in record.INFO:
                posB = record.INFO['END']
            else:
                msg = 'Could not parse SV end position\n'
                msg = msg + 'Sample: %s\tSource: %s\tRecord: %s' % (sample, source, name)
                raise NonstandardVCFError(msg)

        super().__init__(chrA, posA, chrB, posB, name, sample, source, svtype)

    @property
    def caller_svtype(self):
        return self.record.INFO['SVTYPE']

    def __hash__(self):
        return id(self)


class PindelCall(VCFCall):
    def __init__(self, record):
        super().__init__(record, source='pindel')

    @property
    def svcf(self):
        return OrderedDict([
            ('RR_pindel', self.record.samples[0].data.AD[0]),
            ('RV_pindel', self.record.samples[0].data.AD[1])])


class LumpyCall(VCFCall):
    def __init__(self, record):
        # Parse B chrom/pos
        if record.INFO['SVTYPE'] == 'BND':
            chrB = str(record.ALT[0].chr)
            posB = record.ALT[0].pos
        else:
            chrB = str(record.CHROM)
            posB = record.INFO['END']

        # Convert svtype to common string
        svtype = record.INFO['SVTYPE']
        if svtype in STD_SVTYPES:
            svtype = STD_SVTYPES[svtype]
        elif svtype == 'CNV':
            svtype = 'del'
        elif svtype == 'BND':
            chrA = str(record.CHROM)
            if chrA == chrB:
                svtype = 'inv'
            else:
                svtype = 'tloc'
        else:
            raise UnparsableSVTypeError(svtype)

        super().__init__(record, chrB, posB=posB, source='lumpy', svtype=svtype)

    @property
    def secondary(self):
        return ('SECONDARY' in self.record.INFO and
                self.record.INFO['SECONDARY'])

    @property
    def HQ(self):
        """
        By default, variants with both PE and SR support are considered HQ
        """

        qual_filter = None
        if qual_filter is not None:
            return qual_filter(self.record)
        else:
            return (self.record.samples[0].data.PE > 0 and
                    self.record.samples[0].data.SR > 0)

    @property
    def svcf(self):
        """Track LUMPY PE and SR counts"""
        #  return ':'.join(str(lrec.samples[0].data.PE),
                        #  str(lrec.samples[0].data.SR))

        return OrderedDict([
            ('PE_lumpy', self.record.samples[0].data.PE),
            ('SR_lumpy', self.record.samples[0].data.SR)])

    def __hash__(self):
        return id(self)


class DellyCall(VCFCall):
    def __init__(self, record):
        """Delly has non standard chrB"""
        chrB = str(record.INFO['CHR2'])
        posB = int(record.INFO['END'])

        super().__init__(record, chrB=chrB, posB=posB, source='delly')

    @property
    def HQ(self):
        """
        By default, variants with both PE and SR support and a PASS filter
        (min 3 reads, mapq > 20) are considered HQ
        """

        qual_filter = None
        if qual_filter is not None:
            return qual_filter(self.record)
        else:
            record = self.record
            return (len(record.FILTER) == 0 and
                    ('PRECISE' in record.INFO and record.INFO['PRECISE']))

    @property
    def svsize(self):
        if self.caller_svtype == 'INS':
            return self.record.INFO['INSLEN']
        else:
            return super().svsize

    @property
    def svcf(self):
        #  PE (DR, DV), SR (RR, RV), MAPQ
        rec = self.record

        def _get_info(record, key):
            return record.INFO.get(key) or 0

        return OrderedDict([
            ('PE_delly', _get_info(rec, 'PE')),
            ('SR_delly', _get_info(rec, 'SR')),
            ('MAPQ_delly', _get_info(rec, 'MAPQ')),
            ('DR_delly', rec.samples[0].data.DR),
            ('DV_delly', rec.samples[0].data.DV),
            ('RR_delly', rec.samples[0].data.RR),
            ('RV_delly', rec.samples[0].data.RV)])

    def joint_svcf(self, others):
        """Merge genotype fields of overlapping calls"""
        # TODO: take max instead of sum?

        s1 = self.svcf

        mapq_key = 'MAPQ_delly'
        mapqs = [s1['MAPQ_delly']]
        for other in others:
            s2 = other.svcf
            for key in s1:
                if key == mapq_key:
                    mapqs.append(s2[key])
                else:
                    s1[key] += s2[key]

        mapq = np.mean(mapqs)
        s1[mapq_key] = mapq

        return s1

    def __hash__(self):
        return id(self)


class BedpeCall(SVCall):
    def __init__(self, record, sample=None, source=None, svtype=None):

        self.record = record
        chrA = str(record.chrA)
        posA = record.posA
        chrB = str(record.chrB)
        posB = record.posB
        name = record.name

        if svtype is None:
            svtype = record.svtype

        super().__init__(chrA, posA, chrB, posB, name, sample, source, svtype)

    def __hash__(self):
        return id(self)


class RPCCall(BedpeCall):
    def __init__(self, record, sample=None, svtype=None):

        svtype = record.svtype
        if svtype == '.':
            svtypes = 'del dup inv tloc'.split()
            svtypes = [s for s in svtypes if s in record.name]
            if len(svtypes) > 0:
                svtype = svtypes[0]
            else:
                msg = 'RPC file: %s, %s' % (sample, record.name)
                raise UnparsableSVTypeError(msg)

        if sample is None:
            sample = record.name.split('_%s_' % svtype)[0]

        super().__init__(record, sample=sample, source='rpc', svtype=svtype)

        self.cluster_size = int(record.INFO['INFO'][0])

    def __hash__(self):
        return id(self)

    @property
    def HQ(self):
        return self.cluster_size >= 5

    @property
    def svcf(self):
        return OrderedDict([
            ('PE_rpc', int(self.record.INFO['INFO'][0])),
            ('QualA_rpc', float(self.record.INFO['INFO'][1])),
            ('QualB_rpc', float(self.record.INFO['INFO'][2])),
            ('UniqA_rpc', float(self.record.INFO['INFO'][3])),
            ('UniqB_rpc', float(self.record.INFO['INFO'][4]))
        ])


class BedCall(SVCall):
    def __init__(self, interval, sample=None):
        """
        pybedtools.Interval
        """
        self.record = interval
        chrA = str(interval.chrom)
        posA = interval.start
        chrB = str(interval.chrom)
        posB = interval.end

        super().__init__(chrA, posA, chrB, posB, name='.')

    def __hash__(self):
        return id(self)
