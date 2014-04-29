#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


def main():
    setup(name = "hifive",
          version = "1.0.0",
          description = 'Python library for normalizing and analyzing HiC and 5C data',
          package_dir = {'':'lib'},
          packages = ['hifive','hifive.hic', 'hifive.fivec'],
          ext_modules = get_extension_modules(),
          cmdclass = {'build_ext': build_ext},
          requires = ['numpy', 'scipy', 'h5py'],
          author = "Michael Sauria",
          author_email = "mike.sauria@gmail.com",
          url='https://bitbucket.org/bxlab/hifive')

# ---- Extension Modules ----------------------------------------------------

def get_extension_modules():
    extensions = []
    # Distance functions
    extensions.append(Extension("hifive.hic._distance", ["lib/hifive/hic/_distance.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    extensions.append(Extension("hifive.fivec._distance", ["lib/hifive/fivec/_distance.pyx",
                                "lib/hifive/fivec/_normal.cpp"],
                                include_dirs=[numpy.get_include()], language="c++"))
    # Binning functions
    extensions.append(Extension("hifive.hic._binning", ["lib/hifive/hic/_binning.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    extensions.append(Extension("hifive.fivec._binning", ["lib/hifive/fivec/_binning.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    # Boundary Index functions
    extensions.append(Extension("hifive._bi", ["lib/hifive/_bi.pyx"],
                                include_dirs=[numpy.get_include()], language="c++"))
    return extensions

if __name__ == "__main__":
    main()