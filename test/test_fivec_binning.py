#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from operator import mul

import numpy

from hifive import fivec
import h5py


class FiveCBinning(unittest.TestCase):
    def setUp(self):
        self.project_fname = 'test/data/test_probability.fcp'
        self.binned_fname = 'test/data/test_binned.hdf5'
        self.heatmap_fname = 'test/data/test.fch'
        self.data = h5py.File(self.binned_fname, 'r')
        self.heatmap = h5py.File(self.heatmap_fname, 'r')
        self.project = fivec.FiveC(self.project_fname, 'r', silent=True)

    def test_cis_binning(self):
        cis, mapping = self.project.cis_heatmap(0, start=50000, stop=200000, binsize=25000, datatype='enrichment',
                                                arraytype='upper', skipfiltered=False, returnmapping=True,
                                                dynamically_binned=True, minobservations=10, searchdistance=0,
                                                expansion_binsize=0, removefailed=False)
        self.compare_arrays(self.data['fivec_cis'][...], cis, 'cis binned enrichments')
        self.compare_arrays(self.data['fivec_mapping'][...], mapping, 'cis mappings')

    def test_trans_binning(self):
        trans, mapping1, mapping2 = self.project.trans_heatmap(0, 1, binsize=50000, datatype='enrichment',
                                                               returnmapping=True, dynamically_binned=True,
                                                               minobservations=100, searchdistance=0,
                                                               expansion_binsize=50000, removefailed=False)
        self.compare_arrays(self.data['fivec_trans'][...], trans, 'trans binned enrichments')
        self.compare_arrays(self.data['fivec_mapping1'][...], mapping1, 'trans mappings')
        self.compare_arrays(self.data['fivec_mapping2'][...], mapping2, 'trans mappings')

    def test_generate_heatmap(self):
        subprocess.call("./bin/hifive 5c-heatmap -q -b 50000 -t -d fragment -a full %s test/data/test_temp.fch" %
                        self.project_fname, shell=True)
        heatmap = h5py.File("test/data/test_temp.fch")
        self.compare_hdf5_dicts(self.heatmap, heatmap, 'heatmap')

    def tearDown(self):
        subprocess.call('rm -f test/data/test_temp.fch', shell=True)

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