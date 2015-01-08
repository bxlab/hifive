#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from operator import mul

import numpy

from hifive import hic
import h5py


class HiCBinning(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.project_fname = '%s/tests/data/test_analyzed.hcp' % self.basedir
        self.binned_fname = '%s/tests/data/test_binned.hdf5' % self.basedir
        self.heatmap_fname = '%s/tests/data/test.hch' % self.basedir
        self.data = h5py.File(self.binned_fname, 'r')
        self.heatmap = h5py.File(self.heatmap_fname, 'r')
        self.project = hic.HiC(self.project_fname, 'r', silent=True)

    def test_cis_binning(self):
        cis, mapping = self.project.cis_heatmap('1', start=50000, stop=200000, binsize=25000, datatype='enrichment',
                                                arraytype='compact', maxdistance=0, skipfiltered=False,
                                                returnmapping=True, dynamically_binned=True, minobservations=10,
                                                searchdistance=0, expansion_binsize=0, removefailed=False)
        self.compare_arrays(self.data['hic_cis'][...], cis, 'cis binned enrichments')
        self.compare_arrays(self.data['hic_mapping'][...], mapping, 'cis mappings')

    def test_trans_binning(self):
        trans, mapping1, mapping2 = self.project.trans_heatmap('1', '2', binsize=50000, datatype='enrichment',
                                                               returnmapping=True, dynamically_binned=True,
                                                               minobservations=100, searchdistance=0,
                                                               expansion_binsize=50000, removefailed=False)
        self.compare_arrays(self.data['hic_trans'][...], trans, 'trans binned enrichments')
        self.compare_arrays(self.data['hic_mapping1'][...], mapping1, 'trans mappings')
        self.compare_arrays(self.data['hic_mapping2'][...], mapping2, 'trans mappings')

    def test_generate_heatmap(self):
        self.project.write_heatmap_dict("%s/tests/data/test_temp.hch" % self.basedir, 50000, includetrans=True,
                                        remove_distance=False)
        heatmap = h5py.File("%s/tests/data/test_temp.hch" % self.basedir)
        for name in self.heatmap['/'].attrs.keys():
            self.assertTrue(name in heatmap['/'].attrs,
                "%s missing from heatmap attributes" % name)
            self.assertTrue(self.heatmap['/'].attrs[name] == heatmap['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.heatmap.keys():
            self.assertTrue(name in heatmap,
                "%s missing from heatmap arrays" % name)
            self.compare_arrays(self.heatmap[name][...], heatmap[name][...], name)

    def tearDown(self):
        subprocess.call('rm -f %s/tests/data/test_temp.hch' % self.basedir, shell=True)

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