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
        self.basedir = os.path.abspath(os.path.dirname(sys.argv[0]))
        self.data_fname = '%s/tests/data/test.hcd' % self.basedir
        self.project_fname = '%s/tests/data/test.hcp' % self.basedir
        self.analyzed_fname = '%s/tests/data/test_analyzed.hcp' % self.basedir
        self.express_fname = '%s/tests/data/test_express.hcp' % self.basedir
        self.data = h5py.File(self.project_fname, 'r')
        self.analyzed = hic.HiC(self.analyzed_fname, 'r')
        self.express = hic.HiC(self.express_fname, 'r')

    def test_hic_project_creation(self):
        project = hic.HiC('%s/tests/data/test_temp.hcp' % self.basedir, 'w', silent=True)
        project.load_data(self.data_fname)
        project.save()
        project = h5py.File('%s/tests/data/test_temp.hcp' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            self.assertTrue(name in project['/'].attrs,
                "%s missing from project attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == project['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in project,
                "%s missing from project arrays" % name)
            self.compare_arrays(self.data[name][...], project[name][...], name)

    def test_hic_project_preanalysis(self):
        project = hic.HiC(self.project_fname, 'r', silent=True)
        project.filter_fends(mininteractions=10, mindistance=20000, maxdistance=1000000)
        self.assertTrue(numpy.allclose(self.analyzed.filter, project.filter),
            "filtered fends don't match target values")
        project.find_distance_means(numbins=5, minsize=30000, maxsize=1000000)
        self.assertTrue(numpy.allclose(self.analyzed.distance_parameters, project.distance_parameters),
            "distance parameters don't match target values")
        self.assertTrue(numpy.allclose(self.analyzed.chromosome_means, project.chromosome_means),
            "chromosome means don't match target values")

    def test_hic_project_analysis(self):
        project = hic.HiC(self.analyzed_fname, 'r', silent=True)
        project.corrections.fill(1.0)
        project.find_fend_corrections(mindistance=10000, maxdistance=1000000, minchange=0.0015, burnin_iterations=100,
                                      annealing_iterations=10, learningrate=0.4, display=0, maxiterations=100)
        self.assertTrue(numpy.allclose(self.analyzed.corrections, project.corrections),
            "learned correction values don't match target values")

    def test_hic_project_express(self):
        project = hic.HiC(self.express_fname, 'r', silent=True)
        project.corrections.fill(1.0)
        project.find_express_fend_corrections(mindistance=10000, iterations=100)
        self.assertTrue(numpy.allclose(self.express.corrections, project.corrections),
            "learned express correction values don't match target values")

    def tearDown(self):
        subprocess.call('rm -f %s/tests/data/test_temp.hcp' % self.basedir, shell=True)

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