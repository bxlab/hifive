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
        self.data_fname = '%s/test/data/test.hcd' % self.basedir
        self.project_fname = '%s/test/data/test.hcp' % self.basedir
        self.probability_fname = '%s/test/data/test_probability.hcp' % self.basedir
        self.express_fname = '%s/test/data/test_express.hcp' % self.basedir
        self.binning_fname = '%s/test/data/test_binning.hcp' % self.basedir
        self.data = h5py.File(self.project_fname, 'r')
        self.probability = hic.HiC(self.probability_fname, 'r')
        self.express = hic.HiC(self.express_fname, 'r')
        self.binning = hic.HiC(self.binning_fname, 'r')

    def test_hic_project_preanalysis(self):
        project = hic.HiC('%s/test/data/test_temp.hcp' % self.basedir, 'w', silent=True)
        project.load_data(self.data_fname)
        project.filter_fends(mininteractions=10, mindistance=20000, maxdistance=1000000)
        project.find_distance_parameters(numbins=5, minsize=30000, maxsize=1000000)
        project.save()
        project = h5py.File('%s/test/data/test_temp.hcp' % self.basedir, 'r')
        for name in self.data['/'].attrs.keys():
            if name == 'history':
                continue
            self.assertTrue(name in project['/'].attrs,
                "%s missing from project attributes" % name)
            self.assertTrue(self.data['/'].attrs[name] == project['/'].attrs[name],
                "%s doesn't match target value" % name)
        for name in self.data.keys():
            self.assertTrue(name in project,
                "%s missing from project arrays" % name)
            self.compare_arrays(self.data[name][...], project[name][...], name)

    def test_hic_project_probability(self):
        project = hic.HiC(self.project_fname, 'r', silent=True)
        project.find_probability_fend_corrections(mindistance=10000, maxdistance=1000000, minchange=0.0015,
                                                  burnin_iterations=100, annealing_iterations=10, learningrate=0.4,
                                                  display=0)
        self.assertTrue(numpy.allclose(self.probability.corrections, project.corrections),
            "learned correction values don't match target values")
        self.assertTrue(numpy.allclose(self.probability.chromosome_means, project.chromosome_means),
            "chromosome means don't match target values")

    def test_hic_project_express(self):
        project = hic.HiC(self.project_fname, 'r', silent=True)
        project.find_express_fend_corrections(mindistance=10000, iterations=100)
        for i in range(project.corrections.shape[0]):
            if project.filter[i] == 1:
                print self.express.corrections[i], project.corrections[i]
        self.assertTrue(numpy.allclose(self.express.corrections, project.corrections),
            "learned express correction values don't match target values")
        self.assertTrue(numpy.allclose(self.express.chromosome_means, project.chromosome_means),
            "chromosome means don't match target values")

    def test_hic_project_binning(self):
        project = hic.HiC(self.project_fname, 'r', silent=True)
        project.find_binning_fend_corrections(model=['len','distance'], num_bins=[3, 3], usereads='cis',
                                                 max_iterations=5)
        self.assertTrue(numpy.allclose(self.binning.binning_corrections, project.binning_corrections),
            "learned binning correction values don't match target values")

    def tearDown(self):
        subprocess.call('rm -f %s/test/data/test_temp.hcp' % self.basedir, shell=True)

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