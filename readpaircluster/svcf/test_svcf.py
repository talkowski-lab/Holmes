#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2016 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""

"""

from __future__ import print_function
import unittest
try:
    unittest.skip
except AttributeError:
    import unittest2 as unittest

import vcf
import pysam
import rpc
import svcall
import svfile


class TestHelpers(unittest.TestCase):
    def test_is_smaller_chrom(self):
        one = '1'
        two = '2'
        ten = '10'
        x = 'X'
        y = 'Y'

        self.assertTrue(rpc.is_smaller_chrom(one, two))
        self.assertTrue(rpc.is_smaller_chrom(one, one))
        self.assertFalse(rpc.is_smaller_chrom(one, one, le=False))
        self.assertFalse(rpc.is_smaller_chrom(two, one))
        self.assertTrue(rpc.is_smaller_chrom(two, ten))
        self.assertTrue(rpc.is_smaller_chrom(ten, x))
        self.assertTrue(rpc.is_smaller_chrom(x, y))
        self.assertFalse(rpc.is_smaller_chrom(y, x))
        self.assertTrue(rpc.is_smaller_chrom(x, x))
        self.assertFalse(rpc.is_smaller_chrom(x, x, le=False))


class TestSVCall(unittest.TestCase):
    def test_delly(self):
        reader = vcf.Reader(filename='example.delly.vcf')
        record = next(reader)
        call = svcall.DellyCall(record)

        self.assertEqual(call.chrA, '1')
        self.assertEqual(call.posA, 869478)
        self.assertEqual(call.chrB, '1')
        self.assertEqual(call.posB, 870222)
        self.assertEqual(call.name, 'DEL00000001')
        self.assertFalse(call.is_HQ())

        # Test HQ variant
        for record in reader:
            if record.ID == 'DUP00000003':
                break
        call = svcall.DellyCall(record)
        self.assertTrue(call.is_HQ())

    def test_lumpy(self):
        reader = vcf.Reader(filename='example.lumpy.vcf.gz', compressed=True)
        record = next(reader)
        call = svcall.LumpyCall(record)
        self.assertEqual(call.chrA, '1')
        self.assertEqual(call.posA, 869476)
        self.assertEqual(call.chrB, '1')
        self.assertEqual(call.posB, 870221)
        self.assertEqual(call.svtype, 'DEL')
        self.assertFalse(call.is_HQ())

        # Test HQ variant
        for record in reader:
            if record.ID == '101':
                break
        call = svcall.LumpyCall(record)
        self.assertEqual(call.svtype, 'INV')
        self.assertTrue(call.is_HQ())

        # Test secondary
        for record in reader:
            if record.ID.endswith('_2'):
                break
        call = svcall.LumpyCall(record)
        self.assertTrue(call.is_secondary)
        self.assertEqual(call.svtype, 'BND')
        self.assertFalse(call.is_HQ())

    def test_rpc(self):
        reader = svfile.BedpeParser('example.rpc.bedpe')

        # Test parsing
        record = next(reader)
        call = svcall.RPCCall(record)
        self.assertEqual(call.chrA, '1')
        self.assertEqual(call.posA, 869575)
        self.assertEqual(call.chrB, '1')
        self.assertEqual(call.posB, 870221)
        self.assertEqual(call.name, 'SFARI_d11194p1_del_1')
        self.assertEqual(call.sample, 'SFARI_d11194p1')

    def test_equality(self):
        c1 = svcall.SVCall('1', 1, '1', 10, 'c1')
        c1b = svcall.SVCall('1', 1, '1', 10, 'c1')
        c2 = svcall.SVCall('1', 1, '2', 10, 'c1')
        c3 = svcall.SVCall('1', 1, '1', 11, 'c1')

        self.assertEqual(c1, c1b)
        self.assertNotEqual(c1, c2)
        self.assertNotEqual(c1, c3)

    def test_inequality(self):
        c1 = svcall.SVCall('1', 1, '1', 10, 'c1')
        c2 = svcall.SVCall('1', 2, '1', 9, 'c2')
        c3 = svcall.SVCall('1', 1, '1', 8, 'c3')

        c4 = svcall.SVCall('10', 1, '10', 10, 'c4')
        c5 = svcall.SVCall('2', 1, '2', 10, 'c5')
        c6 = svcall.SVCall('X', 1, 'X', 10, 'c6')
        c7 = svcall.SVCall('X', 1, 'Y', 10, 'c7')

        self.assertTrue(c1 < c2)
        self.assertTrue(c1 <= c3)
        self.assertTrue(c1 < c4)
        self.assertTrue(c1 < c5)
        self.assertTrue(c5 < c4)
        self.assertTrue(c1 < c6)
        self.assertTrue(c5 < c7)

    def test_clusters(self):
        dist = 500
        c1 = svcall.SVCall('1', 1, '1', 1000, 'c1')

        c2 = svcall.SVCall('2', 1, '2', 1000, 'c2')
        self.assertFalse(c1.clusters_with(c2, dist))

        c3 = svcall.SVCall('1', 550, '1', 1001, 'c3')
        self.assertFalse(c1.clusters_with(c3, dist))

        c4 = svcall.SVCall('1', 450, '1', 1001, 'c4')
        self.assertTrue(c1.clusters_with(c4, dist))

        c6 = svcall.SVCall('1', 450, '1', 1505, 'c6')
        self.assertTrue(c1.is_clusterable_with(c6, dist))
        self.assertFalse(c1.clusters_with(c6, dist))

    def test_is_in(self):
        tabixfile = pysam.TabixFile('example.tabix.bed.gz')

        c = svcall.SVCall('1', 100, '1', 20000, 'c')
        self.assertTrue(c.is_in(tabixfile))

        c = svcall.SVCall('1', 4000, '1', 6000, 'c')
        self.assertTrue(c.is_in(tabixfile))

        c = svcall.SVCall('2', 4000, '1', 6000, 'c')
        self.assertTrue(c.is_in(tabixfile))

        c = svcall.SVCall('2', 4000, '2', 6000, 'c')
        self.assertFalse(c.is_in(tabixfile))

        c = svcall.SVCall('X', 4000, '2', 6000, 'c')
        self.assertFalse(c.is_in(tabixfile))

        c = svcall.SVCall('X', 400, '2', 6000, 'c')
        self.assertTrue(c.is_in(tabixfile))


class TestSVCallCluster(unittest.TestCase):
    def setUp(self):
        reader = vcf.Reader(filename='example.lumpy.vcf.gz', compressed=True)
        _ = next(reader)
        r1 = next(reader)
        r2 = next(reader)

        c1 = svcall.LumpyCall(r1)
        c2 = svcall.LumpyCall(r2)

        self.cluster = svcall.SVCallCluster([c1, c2], name='merged')

    def test_init(self):
        self.assertEqual(self.cluster.posA, 964450)
        self.assertEqual(self.cluster.posB, 964890)

        self.assertTrue('SFARI_d11194p1' in self.cluster.obs)

        obs_list = ['lumpy', 'lumpy']
        self.assertTrue(self.cluster.obs['SFARI_d11194p1'] == obs_list)

    def test_support_str(self):
        samples = 'SFARI_d11194p1 SFARI_d11194s1'.split()
        progs = 'delly lumpy rpc'.split()

        support_str = self.cluster.support_str(samples, progs)
        p1 = '0:2:0'
        s1 = '0:0:0'

        self.assertEqual(support_str.split()[-1], s1)
        self.assertEqual(support_str.split()[-2], p1)



# class TestSVFileParsing(unittest.TestCase):
    # """Test VCF parsing"""

    # def test_vcf(self):
        # sv_file = svfile.SVFile('example.delly.vcf')
        # call = next(sv_file)

        # self.assertEqual(call.chrA, '1')
        # self.assertEqual(call.posA, 869478)
        # self.assertEqual(call.chrB, '1')
        # self.assertEqual(call.posB, 870222)
        # self.assertEqual(call.name, 'DEL00000001')

    # def test_gzipped_vcf(self):
        # svfile = svfile.SVFile('example.lumpy.vcf.gz')
        # call = svfile.next()

        # self.assertEqual(call.chrA, '1')
        # self.assertEqual(call.posA, 869476)
        # self.assertEqual(call.chrB, '1')
        # self.assertEqual(call.posB, 870221)
        # self.assertEqual(call.name, '1')


if __name__ == '__main__':
    unittest.main(warnings='ignore')
