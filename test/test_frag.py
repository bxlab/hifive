#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from operator import mul

import numpy

from hifive import fragment
import h5py


class Frags(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.frags = h5py.File('%s/tests/data/test.frags' % self.basedir, 'r')
        self.bed_fname = '%s/tests/data/test_frag.bed' % self.basedir

    def test_fragment_creation(self):
        frags = fragment.Fragment('%s/tests/data/test_temp.frags' % self.basedir, 'w', silent=True)
        frags.load_fragments(self.bed_fname, genome_name="mm9", re_name='HindIII')
        frags.save()
        frags = h5py.File('%s/tests/data/test_temp.frags' % self.basedir, 'r')
        for name in self.frags['/'].attrs.keys():
            if name == 'history':
                continue
            self.assertTrue(name in frags['/'].attrs,
                "%s missing from fragment attributes" % name)
            self.assertTrue(self.frags['/'].attrs[name] == frags['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.frags.keys():
            self.assertTrue(name in frags,
                "%s missing from fragment arrays" % name)
            self.compare_arrays(self.frags[name][...], frags[name][...], name)

    def tearDown(self):
        subprocess.call('rm -f %s/tests/data/test_temp.frags' % self.basedir, shell=True)

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