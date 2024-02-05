#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from functools import reduce
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
        self.mat_data = h5py.File('tests/data/test_import_mat.hcd', 'r')
        self.raw_data = h5py.File('tests/data/test_import_raw.hcd', 'r')
        self.bam_data = h5py.File('tests/data/test_import_bam.hcd', 'r')
        self.bin_mat_data = h5py.File('tests/data/test_import_bin_mat.hcd', 'r')
        self.bin_raw_data = h5py.File('tests/data/test_import_bin_raw.hcd', 'r')
        self.bin_bam_data = h5py.File('tests/data/test_import_bin_bam.hcd', 'r')
        self.matrix_data = h5py.File('tests/data/test_import_matrix.hcd', 'r')
        self.fend_fname = 'tests/data/test.fends'
        self.binned_fend_fname = 'tests/data/test_binned.fends'
        self.mat_fname = 'tests/data/test.mat'
        self.raw_fname = 'tests/data/test.raw'
        self.bam_fname1 = 'tests/data/test_hic_1.bam'
        self.bam_fname2 = 'tests/data/test_hic_2.bam'
        self.matrix_fname = 'tests/data/*.test.matrix'

    def test_hic_raw_data_creation(self):
        subprocess.call("./bin/hifive hic-data -q -R %s -i 500 %s tests/data/test_temp.hcd" %
                        (self.raw_fname, self.fend_fname), shell=True)
        data = h5py.File('tests/data/test_temp.hcd', 'r')
        self.compare_hdf5_dicts(self.raw_data, data, 'data')

    def test_hic_mat_data_creation(self):
        subprocess.call("./bin/hifive hic-data -q -M %s -i 500 %s tests/data/test_temp.hcd" %
                        (self.mat_fname, self.fend_fname), shell=True)
        data = h5py.File('tests/data/test_temp.hcd', 'r')
        self.compare_hdf5_dicts(self.mat_data, data, 'data')

    def test_hic_bam_data_creation(self):
        if 'pysam' not in sys.modules.keys():
            print >> sys.stderr, "pysam required for bam import"
            return None
        subprocess.call("./bin/hifive hic-data -q -S %s %s -i 500 %s tests/data/test_temp.hcd" %
                        (self.bam_fname1, self.bam_fname2, self.fend_fname), shell=True)
        data = h5py.File('tests/data/test_temp.hcd', 'r')
        self.compare_hdf5_dicts(self.bam_data, data, 'data')

    def test_hic_bin_raw_data_creation(self):
        subprocess.call("./bin/hifive hic-data -q -R %s -i 500 %s tests/data/test_temp.hcd" %
                        (self.raw_fname, self.binned_fend_fname), shell=True)
        data = h5py.File('tests/data/test_temp.hcd', 'r')
        self.compare_hdf5_dicts(self.bin_raw_data, data, 'data')

    def test_hic_bin_mat_data_creation(self):
        subprocess.call("./bin/hifive hic-data -q -M %s -i 500 %s tests/data/test_temp.hcd" %
                        (self.mat_fname, self.binned_fend_fname), shell=True)
        data = h5py.File('tests/data/test_temp.hcd', 'r')
        self.compare_hdf5_dicts(self.bin_mat_data, data, 'data')

    def test_hic_bin_bam_data_creation(self):
        if 'pysam' not in sys.modules.keys():
            print >> sys.stderr, "pysam required for bam import"
            return None
        subprocess.call("./bin/hifive hic-data -q -S %s %s -i 500 %s tests/data/test_temp.hcd" %
                        (self.bam_fname1, self.bam_fname2, self.binned_fend_fname), shell=True)
        data = h5py.File('tests/data/test_temp.hcd', 'r')
        self.compare_hdf5_dicts(self.bin_bam_data, data, 'data')

    def test_hic_matrix_data_creation(self):
        subprocess.call("./bin/hifive hic-data -q -X '%s' -i 500 %s tests/data/test_temp.hcd" %
                        (self.matrix_fname, self.binned_fend_fname), shell=True)
        data = h5py.File('tests/data/test_temp.hcd', 'r')
        self.compare_hdf5_dicts(self.matrix_data, data, 'data')

    def test_hic_mat_export(self):
        data = hic_data.HiCData('tests/data/test_import_raw.hcd', 'r', silent=True)
        data.export_to_mat('tests/data/test_temp.mat')
        data1 = []
        with open('tests/data/test.mat', 'r') as fp:
            for line in fp:
                data1.append(line)
        data2 = []
        with open('tests/data/test_temp.mat', 'r') as fp:
            for line in fp:
                data2.append(line)
        self.assertEqual(data1, data2,
                "generated mat file doesn't match original")

    def tearDown(self):
        subprocess.call('rm -f tests/data/test_temp.hcd', shell=True)
        subprocess.call('rm -f tests/data/test_temp.mat', shell=True)

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