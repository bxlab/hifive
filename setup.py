#!/usr/bin/env python

import sys
import os
import subprocess
#import multiprocessing
if sys.version_info[:2] < (2, 7) or sys.version_info[0:2] >= (3, 0):
    raise RuntimeError("Python version 2.7+ is required.")

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
from setuptools import setup, find_packages
import ez_setup
ez_setup.use_setuptools()
from distutils.extension import Extension

MAJOR = 1
MINOR = 4
PATCH = None
ISRELEASED = True
VERSION = '%d.%d' % (MAJOR, MINOR)
if not PATCH is None:
    VERSION += '.%d' % (PATCH)

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

def get_version_info():
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('hifive/version.py'):
        # must be a source distribution, use existing version file
        # load it as a separate module to not load scipy/__init__.py
        import imp
        version = imp.load_source('hifive.version', 'hifive/version.py')
        GIT_REVISION = version.git_revision
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev0+' + GIT_REVISION[:7]

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='hifive/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM HIFIVE SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()

def setup_package():
    write_version_py()

    build_requires = ['setuptools_cython']
    try:
        import numpy
    except:
        build_requires.append('numpy')


    FULLVERSION, GIT_REVISION = get_version_info()

    metadata = dict(
        name = "hifive",
        version = FULLVERSION,
        description = 'Python library for normalizing and analyzing HiC and 5C data',
        zip_safe = False,
        include_package_data = True,
        package_dir = {'':'./'},
        install_requires = ['numpy', 'scipy', 'h5py'],
        setup_requires = build_requires,
        scripts = ['bin/hifive', 'bin/fetch_mrh_data'],
        test_suite = 'nose.collector',
        tests_require = 'nose',
        packages = [],
        ext_modules = [],
        #extras_require = {'pyx':[], 'PIL':[], 'mlpy':[]},
        author = "Michael Sauria",
        author_email = "mike.sauria@jhu.edu",
        url='https://github.com/bxlab/hifive',
        download_url='https://github.com/bxlab/hifive/tarball/%s' % FULLVERSION)

    if len(sys.argv) >= 2 and ('--help' in sys.argv[1:] or
            sys.argv[1] in ('--help-commands', 'egg_info', '--version',
                            'clean')):
        # For these actions, NumPy is not required.
        #
        # They are required to succeed without Numpy for example when
        # pip is used to install HiFive when Numpy is not yet present in
        # the system.
        pass
    else:
        include_dirs = [numpy.get_include()]
        metadata['packages'] = find_packages(exclude=['examples', 'test', 'ez_setup.py'], where='./')
        print metadata['packages']
        metadata['ext_modules'] = get_extension_modules(include_dirs)

    setup(**metadata)

# ---- Extension Modules ----------------------------------------------------

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
    return extensions

if __name__ == "__main__":
    setup_package()
