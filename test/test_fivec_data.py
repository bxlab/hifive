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
        self.data = h5py.File('test/data/test_import.fcd', 'r')
        self.frag_fname = 'test/data/test.frags'
        self.count_fname = 'test/data/test.counts'
        self.bam_fname1 = 'test/data/test_fivec_1.bam'
        self.bam_fname2 = 'test/data/test_fivec_2.bam'

    def test_fivec_counts_data_creation(self):
        subprocess.call("./bin/hifive 5c-data -q -C %s %s test/data/test_temp.fcd" %
                        (self.count_fname, self.frag_fname), shell=True)
        data = h5py.File('test/data/test_temp.fcd', 'r')
        self.compare_hdf5_dicts(self.data, data, 'data')

    def test_fivec_bam_data_creation(self):
        if 'pysam' not in sys.modules.keys():
            print >> sys.stderr, "pysam required for bam import"
            return None
        subprocess.call("./bin/hifive 5c-data -q -B %s %s %s test/data/test_temp.fcd" %
                        (self.bam_fname1, self.bam_fname2, self.frag_fname), shell=True)
        data = h5py.File('test/data/test_temp.fcd', 'r')
        self.compare_hdf5_dicts(self.data, data, 'data')

    def tearDown(self):
        subprocess.call('rm -f test/data/test_temp.fcd', shell=True)

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

    def compare_hdf5_dicts(self, dict1, dict2, name):
        for key in dict1['/'].attrs.keys():
            if key == 'history':
                continue
            self.assertTrue(key in dict2['/'].attrs,
                "%s missing from %s attributes" % (key, name))
            value1 = dict1['/'].attrs[key]
            value2 = dict2['/'].attrs[key]
            if isinstance(value1, float):
                self.assertTrue(
                    numpy.allclose(numpy.array([value1], dtype=numpy.float32), numpy.array([value2], dtype=numpy.float32)),
                    "%s in %s doesn't match target value" % (key, name))
            elif isinstance(value1, int):
                self.assertTrue(
                    numpy.allclose(numpy.array([value1], dtype=numpy.int32), numpy.array([value2], dtype=numpy.int32)),
                    "%s in %s doesn't match target value" % (key, name))
            else:
                self.assertTrue(value1 == value2, "%s in %s doesn't match target value" % (key, name))
        for key in dict1.keys():
            self.assertTrue(key in dict2,
                "%s missing from heatmap arrays" % name)
            self.compare_arrays(dict1[key][...], dict2[key][...], name)
        return None


if __name__ == "__main__":
    unittest.main()
