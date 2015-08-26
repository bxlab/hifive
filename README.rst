.. image:: https://readthedocs.org/projects/bxlab-hifive/badge/?version=latest
  :target: http://bxlab-hifive.readthedocs.org/en/latest/
  :alt: Latest Documentation
  
.. image:: https://travis-ci.org/bxlab/hifive.svg?branch=master
  :target: https://travis-ci.org/bxlab/hifive
  :alt: Latest Build

This repository contains code for the hifive package, a set of tools for
handling HiC and 5C data. This includes managing data from mapped reads, either
in bam, mat, or raw formats. All stages use hdf5 dictionaries for fast access
and minimal memory and storage usage.

This package includes methods for normalizing data from either HiC or 5C
experiments at the fragment-end, or fragment level resolution, respectively.
Once normalized, data can be used for plotting, binning, or other statistical
tests within the package very quickly.

This package makes extensive use of the h5py and numpy packages. In addition, if
mpi4py is installed, several methods are compatible with parallelization,
including HiC normalization and heatmap generation.

Documentation can be found `here <http://hifive.docs.taylorlab.org/en/latest/>`_.

The HiFive binary requires Python 2.7. However, the library can be used by Python 2.6 and 2.7.

Installing HiFive
=============================

The easiest way to get HiFive is using pip::

  > pip install hifive

HiFive can also be obtained from `github <https://github.com/bxlab/hifive/>`_ using the following command::

  > git clone https://github.com/bxlab/hifive.git
or alternatively, download a snapshot of the repository using the following commands::

  > wget https://github.com/bxlab/hifive/tarball/1.0
  > tar -xf hifive_v1.0.tar

HiFive depends on a few packages and has several others that extend its functionality.

Required Packages
-----------------
  * `Scipy <http://www.scipy.org>`_
  * `Numpy <http://www.numpy.org>`_
  * `Cython <http://www.cython.org>`_
  * `h5py <http://www.h5py.org>`_

Recommended Packages
--------------------
  * `Pysam <http://code.google.com/p/pysam/>`_
  * `Pyx <http://pyx.sourceforge.net/>`_
  * `PIL <http://www.pythonware.com/products/pil/>`_
  * `mpi4py <http://mpi4py.scipy.org>`_
  * `mlpy <http://mlpy.sourceforge.net>`_

To install HiFive, simply enter the directory that the repository was cloned or downloaded to and use the following command::

  > python setup.py install

If you wish to install HiFive in a location other than the default, you can use the prefix option::

  > python setup.py install --prefix=/your/desired/path

.. _installing_docs:

Installing Documentation
================================

In order to build HiFive's documentation locally, you need to execute the following command::

  > cd doc
  > make html

This will create the documentation suitable for viewing with a web browser. In order to create a pdf version of the documentation, simply run the following::

  > cd doc
  > make latexpdf
