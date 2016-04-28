.. _Release_History:

*******************
Release History
*******************

.. _1_3:

1.3 (2016-04-28)
-------------------

- Added support for creating non-restriction enzyme-based HiC genome partitions
- Added ability to load RE HiC data into binned genome partition
- Added ability to directly load pre-binned matrix files
- Added support for normalization and plotting of binned HiC data
- Added matrix output option in hic-interval subcommand
- Fixed error in sub-binning cis heatmap arrays
- Compiled domain calling approaches into new hic_domain module
- Renamed _hic_tads to _hic_domains in libraries
- Removed ununsed functions in _hic_domains
- Added explicit variable casting for compatibility with cythoning in Windows 

.. _1_2_2:

1.2.2 (2015-12-15)
-------------------

- Added modified arrowhead transformation domain calling
- Re-implemented modified version of boundary-index domain finding
- Added cis sub-region binning function (no longer needs to extend to the diagonal)

.. _1_2:

1.2 (2015-10-21)
-----------------

- Added detailed statistics for loaded HiC data
- Added testing of various data import functions
- Added option to skip PCR duplicate filtering
- Cleaned up HiC cis binning
- Various bug fixes

.. _1_1:

1.1 (2015-08-26)
-----------------

- Added multi-resolution heatmaps (MRH).
- Added stand-alone script for plotting MRH data.
- Fixed requirements resolution for pip loading without requirements already installed.
- Re-implemented Poisson model for HiC probability algorithm.
- Lowered memory requirement for HiC data loading.
- Added dockerfile for building docker image.

.. _1_0_3:

1.0.3 (2015-04-23)
-------------------

- Expanded functionality of Express to include support to Knight-Ruiz algorithm to both 5C and HiC algorthims.
- Added binary counts option to Express in HiC.
- Added non-logged counts functionality to Express in 5C.
- Added pseudo-count option to Binning in HiC.
- Dramatically sped up heatmap creation for 'fend'-corrected heatmaps using probability or express corrections.
- Added binary distance-dependence calculation alongside standard.

.. _1_0_2:

1.0.2 (2015-04-16)
-------------------

- Fixed Express algorithm stopping criteria.
- Updated Probability algorithm and documentation.
- Fixed bug introduced into HiC binning algorithm.
- Added binning documentation.
- Added Travis automated testing.