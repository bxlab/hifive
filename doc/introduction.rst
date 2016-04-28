.. _introduction:


************
Introduction
************

:mod:`HiFive` is a Python package for normalization and analysis of chromatin structural data produced using either the 5C of HiC assay. The aim of HiFive is to create an efficient framework for handling the large amounts of data associated chromatin interaction studies and provide an easy interface for working with and normalizing these data. There are two ways of utilizing HiFive, either through its binary executable or by using its set of classes and functions for more flexibility in analysis. In order to use HiFive, only two types of information are required. First, you need a record of restriction enzyme fragments in an appropriate format (see :class:`hifive.fragment.Fragment` or :class:`hifive.fend.Fend`) for either the subset of targeted fragments or the whole genome for 5C and HiC, respectively. Second, you need 5C or HiC chromatin interaction data, either in a mapped BAM format or in an appropriate post-processed format (see :class:`hifive.fragment.Fragment` or :class:`hifive.fend.Fend`). HiFive will handle everything else. In addition, because of the large memory and time requirements, HiFive uses modularity to avoid repeating processes such as loading data. This allows easy exploration of analysis parameters at a faster rate.

.. _overview_hifive:

HiFive Overview
======================

::

  hifive [-h] [--version]
              {fragments,5c-data,5c-project,5c-normalize,5c-complete,
              5c-heatmap,5c-interval,combine-5c-replicates,
              fends,hic-data,hic-project,hic-normalize,hic-complete,
              hic-heatmap,hic-interval,combine-hic-replicates}

There are eight major functions for both 5C and HiC data types available as subcommands in HiFive.

5C
---

:fragments:             Create a fragment file from a BED file containing targeted RE fragment data.
:5c-data:               Create a data file from mapped BAM or fragment-pair counts files.
:5c-project:            Create a project file, filter fragments, and estimate distance-dependence.
:5c-normalize:          Find correction parameter values using one of the available algorithms (see normalization_).
:5c-complete:           Perform all of the steps of the subcommands fragments, 5c-data, 5c-project, and 5c-normalization in one command.
:5c-heatmap:            Using an already created 5C project, generate an HDF5-formatted heatmap file and optional image.
:5c-interval:           Using an already created 5C project, generate a tabular genomic-interval file for a specified region and optional image.
:5c-combine-replicates: Combine multiple 5C data files into a single file without needing to reload the data.

.. _normalization: normalization.html

HiC
----

:fends:                   Create a fend file from either a BED or HiCPipe-style fend file containing RE fragment data or create an arbitrarily-binned interval file from chromosome length file.Create a data file from mapped BAM, MAT, or paired coordinate text (RAW) files or from binned matrix files.
:hic-project:             Create a project file, filter fends, and estimate distance-dependence.
:hic-normalize:           Find correction parameter values using one of the available algorithms (see normalization_).
:hic-complete:            Perform all of the steps of the subcommands fends, hic-data, hic-project, and hic-normalization in one command.
:hic-heatmap:             Using an already created HiC project, generate an HDF5-formatted heatmap file and optional image.
:hic-interval:            Using an already created HiC project, generate a tabular genomic-interval file for a specified region and optional image.
:hic-combine-replicates:  Combine multiple HiC data files into a single file without needing to reload the data.

A more detailed explanation of the subcommands and all of the command options is here_.

.. _here: hifive.html

.. _organization_of_the_hifive_package:

Organization of the HiFive package
=========================================

The :mod:`HiFive` package is split into several modules, each serving a specific purpose.

================================  ====================================
Functionality                     Module               
================================  ====================================
Restriction enzyme information    :class:`Fragment <hifive.fragment.Fragment>` / :class:`Fend <hifive.fend.Fend>`
Read counts and orientations      :class:`FiveCData <hifive.fivec_data.FiveCData>` / :class:`HiCData <hifive.hic_data.HiCData>`
Model parameters and filtering    :class:`FiveC <hifive.fivec.FiveC>` / :class:`HiC <hifive.hic.HiC>`
Plotting functions                :mod:`plotting <hifive.plotting>`
================================  ====================================

The classes :class:`Fragment <hifive.fragment.Fragment>`, :class:`Fend <hifive.fend.Fend>`, :class:`FiveCData <hifive.fivec_data.FiveCData>`, :class:`HiCData <hifive.hic_data.HiCData>`, :class:`FiveC <hifive.fivec.FiveC>`, and :class:`HiC <hifive.hic.HiC>` are all available from the top level of the :mod:`HiFive` namespace and can be imported using::

  from hifive import *

at the beginning of the Python program. However, in order to prevent namespace pollution, you may also simply use::

  import hifive

HiFive is organized into a hierarchy of data structures. Each structure represents a set of data that may be shared with any number of structures higher in the hierarchy, thus eliminating redundancy of information. For example, a :class:`Fragment <hifive.fragment.Fragment>` object which contains information about the fragments being interrogated in a 5C experiment can be used for all replicates and conditions that use the same primer scheme. Likewise, a :class:`HiCData <hifive.hic_data.HiCData>` object which contains all of the mapped read information for a specific HiC experiment can be used for multiple analyses with different parameter values. This helps reduce the space these data occupy as well as reduce the time to run multiple analyses since each object need only be created once.

Release History
================

- :ref:`1_3`
- :ref:`1_2_2`
- :ref:`1_2`
- :ref:`1_1`
- :ref:`1_0_3`
- :ref:`1_0_2`
