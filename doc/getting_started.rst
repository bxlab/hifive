.. _getting_started:


***************
Getting started
***************

.. _installing_HiFive:

Installing HiFive
=============================

:mod:`HiFive` can be obtained from `HiFive <https://bitbucket.org/bxlab/hifive/>`_ using the following command::

  > hg clone https://bitbucket.org/bxlab/hifive/

or alternatively, download a snapshot of the repository using the following commands::

  > wget https://bitbucket.org/bxlab/hifive/downloads/hifive_v2.2.tar.bz2
  > tar -xjf hifive_v2.2.tar.bz2

:mod:`HiFive` depends on a few packages and has several others that extend its functionality.

Required Packages
-----------------
  * `Scipy <http://www.scipy.ord>`_
  * `Numpy <http://www.numpy.org>`_
  * `h5py <http://www.h5py.org>`_

Recommended Packages
--------------------
  * `Sphinx <https://pypi.python.org/pypi/Sphinx>`_ for generating local documentation
  * `Pysam <http://code.google.com/p/pysam/>`_ for reading BAM files
  * `Pyx <http://pyx.sourceforge.net/>`_ for generating PDF images
  * `PIL <http://www.pythonware.com/products/pil/>`_ for generating bitmap images
  * `mpi4py <http://mpi4py.scipy.org>`_ for utilizing MPI capabilities of several HiC functions
  * `mlpy <http://mlpy.sourceforge.net>`_ for modeling 3D structure

To install :mod:`HiFive`, simply enter the directory that the repository was cloned or downloaded to and use the following command::

  > python setup.py install

If you wish to install :mod:`HiFive` in a location other than the default, you can use the prefix option::

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