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
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.project_fname = '%s/tests/data/test_probability.fcp' % self.basedir
        self.binned_fname = '%s/tests/data/test_binned.hdf5' % self.basedir
        self.heatmap_fname = '%s/tests/data/test.fch' % self.basedir
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
        self.project.write_heatmap("%s/tests/data/test_temp.fch" % self.basedir, 50000, includetrans=True,
                                   datatype='fragment')
        heatmap = h5py.File("%s/tests/data/test_temp.fch" % self.basedir)
        for name in self.heatmap['/'].attrs.keys():
            if name == 'history':
                continue
            self.assertTrue(name in heatmap['/'].attrs,
                "%s missing from heatmap attributes" % name)
            self.assertTrue(self.heatmap['/'].attrs[name] == heatmap['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.heatmap.keys():
            self.assertTrue(name in heatmap,
                "%s missing from heatmap arrays" % name)
            self.compare_arrays(self.heatmap[name][...], heatmap[name][...], name)

    def tearDown(self):
        subprocess.call('rm -f %s/tests/data/test_temp.fch' % self.basedir, shell=True)

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