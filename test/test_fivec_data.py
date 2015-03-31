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

from hifive import fivec_data
import h5py


class FiveCData(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.data = h5py.File('%s/test/data/test_import.fcd' % self.basedir, 'r')
        self.frag_fname = '%s/test/data/test.frags' % self.basedir
        self.count_fname = '%s/test/data/test.counts' % self.basedir
        self.bam_fname1 = '%s/test/data/test_fivec_1.bam' % self.basedir
        self.bam_fname2 = '%s/test/data/test_fivec_2.bam' % self.basedir

    def test_fivec_counts_data_creation(self):
        data = fivec_data.FiveCData('%s/test/data/test_temp.fcd' % self.basedir, 'w', silent=True)
        data.load_data_from_counts(self.frag_fname, self.count_fname)
        data.save()
        data = h5py.File('%s/test/data/test_temp.fcd' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            if name == 'history':
                continue
            self.assertTrue(name in data['/'].attrs,
                "%s missing from data attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == data['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in data,
                "%s missing from data arrays" % name)
            self.compare_arrays(self.data[name][...], data[name][...], name)

    def test_fivec_bam_data_creation(self):
        if 'pysam' not in sys.modules.keys():
            print >> sys.stderr, "pysam required for bam import"
            return None
        data = fivec_data.FiveCData('%s/test/data/test_temp.fcd' % self.basedir, 'w', silent=True)
        data.load_data_from_bam(self.frag_fname, [self.bam_fname1, self.bam_fname2])
        data.save()
        data = h5py.File('%s/test/data/test_temp.fcd' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            if name == 'history':
                continue
            self.assertTrue(name in data['/'].attrs,
                "%s missing from data attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == data['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in data,
                "%s missing from data arrays" % name)
            self.compare_arrays(self.data[name][...], data[name][...], name)

    def tearDown(self):
        subprocess.call('rm -f %s/test/data/test_temp.fcd' % self.basedir, shell=True)

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
