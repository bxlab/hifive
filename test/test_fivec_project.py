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
        project = fivec.FiveC('%s/test/data/test_temp.fcp' % self.basedir, 'w', silent=True)
        project.load_data(self.data_fname)
        project.filter_fragments(mininteractions=20)
        project.find_distance_parameters()
        project.save()
        project = h5py.File('%s/test/data/test_temp.fcp' % self.basedir, 'r')
        for name in self.raw['/'].attrs.keys():
            if name == 'history':
                continue
            self.assertTrue(name in project['/'].attrs,
                "%s missing from project attributes" % name)
            self.assertTrue(self.raw['/'].attrs[name] == project['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.raw.keys():
            self.assertTrue(name in project,
                "%s missing from project arrays" % name)
            self.compare_arrays(self.raw[name][...], project[name][...], name)

    def test_fivec_project_probability(self):
        project = fivec.FiveC(self.project_fname, 'r', silent=True)
        project.find_probability_fragment_corrections(maxdistance=0, burnin_iterations=100, annealing_iterations=10,
                                                      learningrate=0.01)
        self.assertTrue(numpy.allclose(self.probability.corrections, project.corrections),
            "learned correction values don't match target values")
        self.assertTrue(numpy.allclose(self.probability.region_means, project.region_means),
                                       "region means don't match target values")

    def test_fivec_project_express(self):
        project = fivec.FiveC(self.project_fname, 'r', silent=True)
        project.find_express_fragment_corrections(iterations=100, remove_distance=True)
        self.assertTrue(numpy.allclose(self.express.corrections, project.corrections),
            "learned express correction values don't match target values")
        self.assertTrue(numpy.allclose(self.express.region_means, project.region_means),
                                       "region means don't match target values")

    def test_fivec_project_binning(self):
        project = fivec.FiveC(self.project_fname, 'r', silent=True)
        project.find_binning_fragment_corrections(model=['len'], num_bins=[5], parameters=['even'], usereads='cis',
                                                     max_iterations=0)
        print list(self.binning.binning_corrections)
        print list(project.binning_corrections)
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


if __name__ == "__main__":
    unittest.main()