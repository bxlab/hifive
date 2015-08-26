#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from operator import mul

import numpy

from hifive import hic
import h5py


class HiCProject(unittest.TestCase):
    def setUp(self):
        self.data_fname = 'test/data/test.hcd'
        self.project_fname = 'test/data/test.hcp'
        self.probbin_fname = 'test/data/test_probbin.hcp'
        self.probpois_fname = 'test/data/test_probpois.hcp'
        self.express_fname = 'test/data/test_express.hcp'
        self.binning_fname = 'test/data/test_binning.hcp'
        self.data = h5py.File(self.project_fname, 'r')
        self.probbin = hic.HiC(self.probbin_fname, 'r')
        self.probpois = hic.HiC(self.probpois_fname, 'r')
        self.express = hic.HiC(self.express_fname, 'r')
        self.binning = hic.HiC(self.binning_fname, 'r')

    def test_hic_project_preanalysis(self):
        subprocess.call("./bin/hifive hic-project -q -m 20000 -f 10 -j 30000 -n 5 %s test/data/test_temp.hcp" %
                        (self.data_fname), shell=True)
        project = h5py.File('test/data/test_temp.hcp', 'r')
        self.compare_hdf5_dicts(self.data, project, 'project')

    def test_hic_project_probability_binomial(self):
        subprocess.call("./bin/hifive hic-normalize probability -q -m 20000 -o test/data/test_temp.hcp -b 15 -l 0.4 -g 0.0015 -p %s" %
                        (self.project_fname), shell=True)
        project = hic.HiC("test/data/test_temp.hcp", 'r', silent=True)
        self.assertTrue(numpy.allclose(self.probbin.corrections, project.corrections, atol=1e-4),
            "learned correction values don't match target values")
        self.assertTrue(numpy.allclose(self.probbin.chromosome_means, project.chromosome_means, atol=1e-4),
            "chromosome means don't match target values")

    def test_hic_project_probability_poisson(self):
        subprocess.call("./bin/hifive hic-normalize probability -q -m 20000 -o test/data/test_temp.hcp -b 15 -l 0.4 -g 0.0015 -p -a poisson %s" %
                        (self.project_fname), shell=True)
        project = hic.HiC("test/data/test_temp.hcp", 'r', silent=True)
        self.assertTrue(numpy.allclose(self.probpois.corrections, project.corrections, atol=1e-1),
            "learned correction values don't match target values")
        self.assertTrue(numpy.allclose(self.probpois.chromosome_means, project.chromosome_means, atol=1e-4),
            "chromosome means don't match target values")

    def test_hic_project_express(self):
        subprocess.call("./bin/hifive hic-normalize express -q -m 20000 -o test/data/test_temp.hcp -e 100 -w cis -f 10 %s" %
                        (self.project_fname), shell=True)
        project = hic.HiC("test/data/test_temp.hcp", 'r', silent=True)
        self.assertTrue(numpy.allclose(self.express.corrections, project.corrections),
            "learned express correction values don't match target values")
        self.assertTrue(numpy.allclose(self.express.chromosome_means, project.chromosome_means),
            "chromosome means don't match target values")

    def test_hic_project_binning(self):
        subprocess.call("./bin/hifive hic-normalize binning -q -m 20000 -o test/data/test_temp.hcp -r 5 -y cis -t 1.0 -v len,distance -s 3,3 -u even,fixed-const %s" %
                        (self.project_fname), shell=True)
        project = hic.HiC("test/data/test_temp.hcp", 'r', silent=True)
        self.assertTrue(numpy.allclose(self.binning.binning_corrections, project.binning_corrections),
            "learned binning correction values don't match target values")

    def tearDown(self):
        subprocess.call('rm -f test/data/test_temp.hcp', shell=True)

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