.. _Release_History:

*******************
Release History
*******************

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