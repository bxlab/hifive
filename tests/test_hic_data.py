#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from operator import mul

import numpy
try:
    import pysam
except:
    pass

from hifive import hic_data
import h5py


class HiCData(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.data = h5py.File('%s/tests/data/test_import.hcd' % self.basedir, 'r')
        self.fend_fname = '%s/tests/data/test.fends' % self.basedir
        self.mat_fname = '%s/tests/data/test.mat' % self.basedir
        self.raw_fname = '%s/tests/data/test.raw' % self.basedir
        self.bam_fname1 = '%s/tests/data/test_hic_1.bam' % self.basedir
        self.bam_fname2 = '%s/tests/data/test_hic_2.bam' % self.basedir

    def test_hic_raw_data_creation(self):
        data = hic_data.HiCData('%s/tests/data/test_temp.hcd' % self.basedir, 'w', silent=True)
        data.load_data_from_raw(self.fend_fname, self.raw_fname, 1000)
        data.save()
        data = h5py.File('%s/tests/data/test_temp.hcd' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            self.assertTrue(name in data['/'].attrs,
                "%s missing from data attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == data['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in data,
                "%s missing from data arrays" % name)
            self.compare_arrays(self.data[name][...], data[name][...], name)

    def test_hic_mat_data_creation(self):
        data = hic_data.HiCData('%s/tests/data/test_temp.hcd' % self.basedir, 'w', silent=True)
        data.load_data_from_mat(self.fend_fname, self.mat_fname, 1000)
        data.save()
        data = h5py.File('%s/tests/data/test_temp.hcd' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            self.assertTrue(name in data['/'].attrs,
                "%s missing from data attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == data['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in data,
                "%s missing from data arrays" % name)
            self.compare_arrays(self.data[name][...], data[name][...], name)

    def test_hic_bam_data_creation(self):
        if 'pysam' not in sys.modules.keys():
            print >> sys.stderr, "pysam required for bam import"
            return None
        data = hic_data.HiCData('%s/tests/data/test_temp.hcd' % self.basedir, 'w', silent=True)
        data.load_data_from_bam(self.fend_fname, [self.bam_fname1, self.bam_fname2], 1000)
        data.save()
        data = h5py.File('%s/tests/data/test_temp.hcd' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            self.assertTrue(name in data['/'].attrs,
                "%s missing from data attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == data['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in data,
                "%s missing from data arrays" % name)
            self.compare_arrays(self.data[name][...], data[name][...], name)

    def test_hic_mat_export(self):
        data = hic_data.HiCData('%s/tests/data/test_temp.hcd' % self.basedir, 'w', silent=True)
        data.load_data_from_raw(self.fend_fname, self.raw_fname, 1000)
        data.export_to_mat('%s/tests/data/test_temp.mat' % self.basedir)
        data1 = []
        for line in open('%s/tests/data/test.mat' % self.basedir, 'r'):
            data1.append(line)
        data2 = []
        for line in open('%s/tests/data/test_temp.mat' % self.basedir, 'r'):
            data2.append(line)
        self.assertEqual(data1, data2,
                "generated mat file doesn't match original")

    def tearDown(self):
        subprocess.call('rm -f %s/tests/data/test_temp.hcd' % self.basedir, shell=True)
        subprocess.call('rm -f %s/tests/data/test_temp.mat' % self.basedir, shell=True)

    def compare_arrays(self, array1, array2, name):
        self.assertTrue(array1.shape == array2.shape,
            "%s shape doesn't match target value" % name)
        self.assertTrue(array1.dtype == array2.dtype,
            "%s dtype doesn't match target value" % name)
        if str(array1.dtype).count('S') + str(array1.dtype).count('a') > 0:
            self.assertTrue(numpy.sum(array1 == array2) == reduce(mul, array1.shape),
                "%s don't match target values." % name)
        else:
            self.assertTrue(numpy.allclose(array1, array2),
                "%s don't match target values" % name)
        return None


if __name__ == "__main__":
    unittest.main()