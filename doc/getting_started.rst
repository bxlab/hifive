.. _getting_started:


***************
Getting started
***************

.. _installing_HiFive:

Installing HiFive
=============================

HiFive can bo obtained from `HiFive <https://bitbucket.org/bxlab/hifive/>`_ using the following command::

  > hg clone https://bitbucket.org/bxlab/hifive/

or alternatively, download a snapshot of the repository using the following commands::

  > wget https://bitbucket.org/bxlab/hifive/downloads/hifive_v2.0.tar.bz2
  > tar -xjf hifive_v2.0.tar.bz2

HiFive depends on a few packages and has several others that extend its functionality.

Required Packages
-----------------
  * `Scipy <http://www.scipy.ord>`_
  * `Numpy <http://www.numpy.org>`_
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

This will create the documentation suitale for viewing with a web browser. In order to create a pdf version of the documentation, simply run the following::

  > cd doc
  > make latexpdf