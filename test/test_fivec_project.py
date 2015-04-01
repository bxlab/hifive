#!/usr/bin/env python

import os
import sys
import subprocess
import unittest
from operator import mul

import numpy

from hifive import fivec
import h5py


class FiveCProject(unittest.TestCase):
    def setUp(self):
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.data_fname = '%s/test/data/test.fcd' % self.basedir
        self.project_fname = '%s/test/data/test.fcp' % self.basedir
        self.probability_fname = '%s/test/data/test_probability.fcp' % self.basedir
        self.express_fname = '%s/test/data/test_express.fcp' % self.basedir
        self.binning_fname = '%s/test/data/test_binning.fcp' % self.basedir
        self.raw = h5py.File(self.project_fname, 'r')
        self.probability = fivec.FiveC(self.probability_fname, 'r')
        self.express = fivec.FiveC(self.express_fname, 'r')
        self.binning = fivec.FiveC(self.binning_fname, 'r')

    def test_fivec_project_preanalysis(self):
        subprocess.call("hifive 5c-project -q -f 20 %s %s/test/data/test_temp.fcp" %
                        (self.data_fname, self.basedir), shell=True)
        project = h5py.File('%s/test/data/test_temp.fcp' % self.basedir, 'r')
        self.compare_hdf5_dicts(self.raw, project, 'project')

    def test_fivec_project_probability(self):
        subprocess.call("hifive 5c-normalize probability -q -o %s/test/data/test_temp.fcp -b 100 -a 10 -l 0.01 -p %s" %
                        (self.basedir, self.project_fname), shell=True)
        project = fivec.FiveC("%s/test/data/test_temp.fcp" % self.basedir, 'r', silent=True)
        self.assertTrue(numpy.allclose(self.probability.corrections, project.corrections),
            "learned correction values don't match target values")
        self.assertTrue(numpy.allclose(self.probability.region_means, project.region_means),
                                       "region means don't match target values")

    def test_fivec_project_express(self):
        subprocess.call("hifive 5c-normalize express -q -o %s/test/data/test_temp.fcp -e 100 -d -w cis %s" %
                        (self.basedir, self.project_fname), shell=True)
        project = fivec.FiveC("%s/test/data/test_temp.fcp" % self.basedir, 'r', silent=True)
        self.assertTrue(numpy.allclose(self.express.corrections, project.corrections),
            "learned express correction values don't match target values")
        self.assertTrue(numpy.allclose(self.express.region_means, project.region_means),
                                       "region means don't match target values")

    def test_fivec_project_binning(self):
        subprocess.call("hifive 5c-normalize binning -q -o %s/test/data/test_temp.fcp -i 10 -t 1.0 -y cis -v len -n 5 -u even %s" %
                        (self.basedir, self.project_fname), shell=True)
        project = fivec.FiveC("%s/test/data/test_temp.fcp" % self.basedir, 'r', silent=True)
        self.assertTrue(numpy.allclose(self.binning.binning_corrections, project.binning_corrections),
            "learned binning correction values don't match target values")

    def tearDown(self):
        subprocess.call('rm -f %s/test/data/test_temp.fcp' % self.basedir, shell=True)

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