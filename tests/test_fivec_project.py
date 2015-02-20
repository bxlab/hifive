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
        self.data_fname = '%s/tests/data/test.fcd' % self.basedir
        self.project_fname = '%s/tests/data/test.fcp' % self.basedir
        self.analyzed_fname = '%s/tests/data/test_analyzed.fcp' % self.basedir
        self.express_fname = '%s/tests/data/test_express.fcp' % self.basedir
        self.data = h5py.File(self.project_fname, 'r')
        self.analyzed = fivec.FiveC(self.analyzed_fname, 'r')
        self.express = fivec.FiveC(self.express_fname, 'r')

    def test_fivec_project_creation(self):
        project = fivec.FiveC('%s/tests/data/test_temp.fcp' % self.basedir, 'w', silent=True)
        project.load_data(self.data_fname)
        project.save()
        project = h5py.File('%s/tests/data/test_temp.fcp' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            self.assertTrue(name in project['/'].attrs,
                "%s missing from project attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == project['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in project,
                "%s missing from project arrays" % name)
            self.compare_arrays(self.data[name][...], project[name][...], name)

    def test_fivec_project_preanalysis(self):
        project = fivec.FiveC(self.project_fname, 'r', silent=True)
        project.filter_fragments(mininteractions=20)
        self.assertTrue(numpy.allclose(self.analyzed.filter, project.filter),
            "filtered fragments don't match target values")
        project.find_distance_parameters()
        print self.analyzed.mu, project.mu
        print self.analyzed.gamma, project.gamma
        print self.analyzed.sigma, project.sigma
        self.assertTrue(self.analyzed.mu == project.mu and self.analyzed.gamma == project.gamma and
                        self.analyzed.sigma == project.sigma,
                        "distance parameters don't match target values")

    def test_fivec_project_analysis(self):
        project = fivec.FiveC(self.analyzed_fname, 'r', silent=True)
        project.corrections.fill(0.0)
        project.find_fragment_corrections(maxdistance=0, burnin_iterations=100, annealing_iterations=10,
                                          learningrate=0.1, display=0)
        self.assertTrue(numpy.allclose(self.analyzed.corrections, project.corrections),
            "learned correction values don't match target values")

    def test_fivec_project_express(self):
        project = fivec.FiveC(self.express_fname, 'r', silent=True)
        project.corrections.fill(0.0)
        project.find_express_fragment_corrections(iterations=100, remove_distance=True)
        self.assertTrue(numpy.allclose(self.express.corrections, project.corrections),
            "learned express correction values don't match target values")

    def tearDown(self):
        subprocess.call('rm -f %s/tests/data/test_temp.fcp' % self.basedir, shell=True)

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