.. _getting_started:


***************
Getting started
***************

The :mod:`HiFive` binary requires Python 2.7. However, the library can be used by Python 2.6 and 2.7. To check your version, run::

  > python -V

.. _installing_HiFive:

Installing HiFive
=============================

The easiest way to get :mod:`HiFive` is using pip::

  > pip install hifive

:mod:`HiFive` can also be obtained from `HiFive <https://github.com/bxlab/hifive/>`_ using the following command::

  > git clone https://github.com/bxlab/hifive.git

or alternatively, download a snapshot of the repository using the following commands::

  > wget https://github.com/bxlab/hifive/tarball/1.0
  > tar -xf hifive_v1.0.tar

:mod:`HiFive` depends on a few packages and has several others that extend its functionality.

Required Packages
-----------------
  * `Scipy <http://www.scipy.org>`_
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

For a more traditional Python installation simply setup a virtualevn for :mod:`HiFive` (this example creates one in .venv in your home directory) and then install the required packages and HiFive.

::

  > virtualenv ~/.venv; source ~/.venv/bin/activate

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