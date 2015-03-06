#!/usr/bin/env python

import sys
import multiprocessing

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
from distutils.extension import Extension
import numpy


if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

def main():
    setup(name = "hifive",
          version = "2.2.0",
          description = 'Python library for normalizing and analyzing HiC and 5C data',
          zip_safe = False,
          include_package_data = True,
          package_dir = {'':'lib'},
          packages = find_packages(exclude=['examples', 'tests', 'ez_setup.py'], where='lib/'),
          install_requires = ['numpy', 'scipy', 'h5py'],
          setup_requires = ['setuptools_cython'],
          ext_modules = get_extension_modules(),
          test_suite = 'nose.collector',
          tests_require = 'nose',
          #extras_require = {'pyx':[], 'PIL':[], 'mlpy':[]},
          author = "Michael Sauria",
          author_email = "mike.sauria@gmail.com",
          url='https://bitbucket.org/bxlab/hifive')

# ---- Extension Modules ----------------------------------------------------

def get_extension_modules():
    extensions = []
    # Distance functions
    extensions.append(Extension("hifive.libraries._hic_distance", ["lib/hifive/libraries/_hic_distance.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    extensions.append(Extension("hifive.libraries._fivec_distance", ["lib/hifive/libraries/_fivec_distance.pyx",
                                "lib/hifive/libraries/_normal.cpp"],
                                include_dirs=[numpy.get_include()], language="c++"))
    # Binning functions
    extensions.append(Extension("hifive.libraries._hic_binning", ["lib/hifive/libraries/_hic_binning.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    extensions.append(Extension("hifive.libraries._fivec_binning", ["lib/hifive/libraries/_fivec_binning.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    # Interaction functions
    extensions.append(Extension("hifive.libraries._hic_interactions", ["lib/hifive/libraries/_hic_interactions.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    # Optimization functions
    extensions.append(Extension("hifive.libraries._hic_optimize", ["lib/hifive/libraries/_hic_optimize.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    # Boundary Index functions
    extensions.append(Extension("hifive.libraries._bi", ["lib/hifive/libraries/_bi.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    return extensions

if __name__ == "__main__":
    main()