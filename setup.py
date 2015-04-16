#!/usr/bin/env python

import sys
import multiprocessing

import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
from distutils.extension import Extension
import numpy

MAJOR = 1
MINOR = 0
MICRO = 2
ISRELEASED = False
VERSION = '%d.%d' % (MAJOR, MINOR)
FULLVERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

def write_version_py(filename='hifive/version.py'):
    print "Writing hifive/version.py"
    cnt = """
# THIS FILE IS GENERATED FROM HIFIVE SETUP.PY
version = '%(version)s'
full_version = '%(full_version)s'
release = %(isrelease)s
if not release:
    version = full_version
    """
    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
        'full_version' : FULLVERSION,
        'isrelease': str(ISRELEASED)})
    finally:
        a.close()

def main():
    write_version_py()

    build_requires = ['setuptools_cython']

    metadata = dict(
        name = "hifive",
        version = FULLVERSION,
        description = 'Python library for normalizing and analyzing HiC and 5C data',
        zip_safe = False,
        include_package_data = True,
        package_dir = {'':'./'},
        packages = find_packages(exclude=['examples', 'test', 'ez_setup.py'], where='./'),
        install_requires = ['numpy', 'scipy', 'h5py'],
        setup_requires = build_requires,
        ext_modules = get_extension_modules(),
        scripts = ['bin/hifive',],
        test_suite = 'nose.collector',
        tests_require = 'nose',
        #extras_require = {'pyx':[], 'PIL':[], 'mlpy':[]},
        author = "Michael Sauria",
        author_email = "mike.sauria@jhu.edu",
        url='https://github.com/bxlab/hifive',
        download_url='https://github.com/bxlab/hifive/tarball/%s' % FULLVERSION)
    setup(**metadata)

# ---- Extension Modules ----------------------------------------------------

def get_extension_modules():
    extensions = []
    # Distance functions
    extensions.append(Extension("hifive.libraries._hic_distance", ["hifive/libraries/_hic_distance.pyx"],
                                include_dirs=[numpy.get_include()], language="c++",
                                extra_compile_args=[]))
    # Binning functions
    extensions.append(Extension("hifive.libraries._hic_binning", ["hifive/libraries/_hic_binning.pyx"],
                                include_dirs=[numpy.get_include()], language="c++",
                                extra_compile_args=[]))
    extensions.append(Extension("hifive.libraries._fivec_binning", ["hifive/libraries/_fivec_binning.pyx"],
                                include_dirs=[numpy.get_include()], language="c++",
                                extra_compile_args=[]))
    # Interaction functions
    extensions.append(Extension("hifive.libraries._hic_interactions", ["hifive/libraries/_hic_interactions.pyx"],
                                include_dirs=[numpy.get_include()], language="c++",
                                extra_compile_args=[]))
    # Optimization functions
    extensions.append(Extension("hifive.libraries._hic_optimize", ["hifive/libraries/_hic_optimize.pyx"],
                                include_dirs=[numpy.get_include()], language="c++",
                                extra_compile_args=[]))
    extensions.append(Extension("hifive.libraries._fivec_optimize", ["hifive/libraries/_fivec_optimize.pyx",
                                "hifive/libraries/_normal.cpp"],
                                include_dirs=[numpy.get_include()], language="c++",
                                extra_compile_args=[]))
    return extensions

if __name__ == "__main__":
    main()
