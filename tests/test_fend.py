#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from operator import mul

import numpy

from hifive import fend
import h5py


class Fends(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.fends = h5py.File('%s/tests/data/test.fends' % self.basedir, 'r')
        self.bed_fname = '%s/tests/data/test_fend.bed' % self.basedir

    def test_fend_creation(self):
        fends = fend.Fend('%s/tests/data/test_temp.fends' % self.basedir, 'w', silent=True)
        fends.load_fends(self.bed_fname, genome_name="hg19", re_name='HindIII')
        fends.save()
        fends = h5py.File('%s/tests/data/test_temp.fends' % self.basedir, 'r')
        for name in self.fends['/'].attrs.keys():
            if name == 'history':
                continue
            self.assertTrue(name in fends['/'].attrs,
                "%s missing from fend attributes" % name)
            self.assertTrue(self.fends['/'].attrs[name] == fends['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.fends.keys():
            self.assertTrue(name in fends,
                "%s missing from fend arrays" % name)
            self.compare_arrays(self.fends[name][...], fends[name][...], name)

    def tearDown(self):
        subprocess.call('rm -f %s/tests/data/test_temp.fends' % self.basedir, shell=True)

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