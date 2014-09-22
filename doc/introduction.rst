.. _introduction:


************
Introduction
************

HiFive is a Python package for normalization and analysis of chromatin structural data produced using either the 5C of HiC assay. This library contains tools for handling all steps after mapping of reads.

.. _organization_of_the_hifive_package:

Organization of the HiFive package
==================================

The HiFive package is split into several modules, each serving a specific purpose.

================================  =====================
Functionality                     Module               
================================  =====================
Restriction enzyme information    Fragment / Fend
Read counts and orientations      FiveCData / HiCData
Model parameters and filtering    FiveC / HiC
Boundary index                    BI
Plotting functions                plotting
================================  =====================

The classes Fragment, Fend, FiveCData, HiCData, FiveC, HiC, and BI are all available from the top level of the HiFive namespace and can be imported using::

  from hifive import *

at the beginning of the Python program. However, in order to prevent namespace pollution, you may also simply use import hifive.

HiFive is organized into a hierarchy of data structures. Each structure represents a set of data that may be shared with any number of structures higher in the hierarchy, thus eliminating redundency of information. For example, a Fragment object which contains information about the fragments being interogated in a 5C experiment can be used for all replicates and conditions that use the same primer scheme. Likewise, a HiCDataset object which contains all of the mapped read information for a specific HiC experiment can be used for multiple analyses with different parameter values. This helps reduce the space these data occupy as well as reduce the time to run mutiple analyses since each object need only be created once.

The organization of structures is as follows:

.. image:: _static/flowchart.png