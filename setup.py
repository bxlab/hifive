#!/usr/bin/env python

from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler import Options
import numpy


def get_extension_modules(include_dirs):
    extensions = []
    # Distance functions
    extensions.append(Extension("hifive.libraries._hic_distance", ["hifive/libraries/_hic_distance.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    # Binning functions
    extensions.append(Extension("hifive.libraries._hic_binning", ["hifive/libraries/_hic_binning.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    extensions.append(Extension("hifive.libraries._fivec_binning", ["hifive/libraries/_fivec_binning.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    # Interaction functions
    extensions.append(Extension("hifive.libraries._hic_interactions", ["hifive/libraries/_hic_interactions.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    # Optimization functions
    extensions.append(Extension("hifive.libraries._hic_optimize", ["hifive/libraries/_hic_optimize.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    extensions.append(Extension("hifive.libraries._fivec_optimize", ["hifive/libraries/_fivec_optimize.pyx",
                                "hifive/libraries/_normal.cpp"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    # Domain functions
    extensions.append(Extension("hifive.libraries._hic_domains", ["hifive/libraries/_hic_domains.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    extensions.append(Extension("hifive.libraries._hmm", ["hifive/libraries/_hmm.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    # Quasar functions
    extensions.append(Extension("hifive.libraries._quasar", ["hifive/libraries/_quasar.pyx"],
                                include_dirs=include_dirs, language="c++",
                                extra_compile_args=[]))
    return extensions

extensions = get_extension_modules([numpy.get_include()])

setup(
    name='hifive',
    ext_modules = cythonize(extensions, compiler_directives={"language_level": 3, "profile": False}),
)
