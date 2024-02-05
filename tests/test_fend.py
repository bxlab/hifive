#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from functools import reduce
from operator import mul

import numpy

from hifive import fend
import h5py


class Fends(unittest.TestCase):
    def setUp(self):
        self.fends = h5py.File('tests/data/test.fends', 'r')
        self.binned_fends = h5py.File('tests/data/test_binned.fends', 'r')
        self.len_fends = h5py.File('tests/data/test_len.fends', 'r')
        self.bed_fname = 'tests/data/test_fend.bed'
        self.len_fname = 'tests/data/test.len'

    def test_fend_creation(self):
        subprocess.call("./bin/hifive fends -q -g mm9 -r HindIII -B %s tests/data/test_temp.fends" %
                        self.bed_fname, shell=True)
        fends = h5py.File('tests/data/test_temp.fends', 'r')
        self.compare_hdf5_dicts(self.fends, fends, 'fends')

    def test_binned_fend_creation(self):
        subprocess.call("./bin/hifive fends -q -g mm9 -B %s tests/data/test_temp.fends --binned 40000" %
                        self.bed_fname, shell=True)
        fends = h5py.File('tests/data/test_temp.fends', 'r')
        self.compare_hdf5_dicts(self.binned_fends, fends, 'RE binned fends')

    def test_binned_creation(self):
        subprocess.call("./bin/hifive fends -q -g mm9 -L %s tests/data/test_temp.fends --binned 40000" %
                        self.len_fname, shell=True)
        fends = h5py.File('tests/data/test_temp.fends', 'r')
        self.compare_hdf5_dicts(self.len_fends, fends, 'binned fends')

    def tearDown(self):
        subprocess.call('rm -f tests/data/test_temp.fends', shell=True)

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
                "%s missing from fend arrays" % key)
            self.compare_arrays(dict1[key][...], dict2[key][...], name)
        return None


if __name__ == "__main__":
    unittest.main()