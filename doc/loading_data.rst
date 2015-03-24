.. _loading_data:

*************************************
Loading Data
*************************************

HiFive data is handled using the :class:`FiveCData <hifive.fivec_data.FiveCData>` and :class:`HiCData <hfiive.hic_data.HiCData>` classes.

.. _fivec_data_loading:

Loading 5C data
===============

HiFive can load 5C data from one of two source file types.

BAM Files
---------

When loading 5C data from BAM files, they should always come in pairs, one for each end of the paired-end reads. HiFive can load any number of pairs of BAM files, such as when multiple sequencing lanes have been run for a single replicate. These files do not need to be indexed or sorted. All sequence names that these files were mapped against should exactly match the primer names in the BED file used to construct the Fragment object.

Count Files
------------

Counts files are tabular text files containing pairs of primer names and a count of the number of observed occurrences of that pairing.

::

  5c_for_primer1   5c_rev_primer2    10
  5c_for_primer1   5c_rev_primer4    3
  5c_for_primer3   5c_rev_primer4    18

.. _hic_data_loading:

Loading HiC Data
================

HiFive can load HiC data from three different types of source files.

BAM Files
---------

When loading HiC data from BAM files, they should always come in pairs, one for each end of the paired-end reads. HiFive can load any number of pairs of BAM files, such as when multiple sequencing lanes have been run for a single replicate. These files do not need to be indexed or sorted. For faster loading, especially with very large numbers of reads, it is helpful to parse out single-mapped reads to reduce the number of reads that HiFive needs to traverse in reading the BAM files.

RAW Files
---------

RAW files are tabular text files containing pairs of read coordinates from mapped reads containing the chromosome, coordinate, and strand for each read end. HiFive can load any number of RAW files into a single HiC Data object.

::

  chr1    30002023    +    chr3    4020235    -
  chr5    9326220     -    chr1    3576222    +
  chr8    1295363     +    chr6    11040321   +

MAT Files
---------

MAT files are in a tabular text format previously defined for `HiCPipe <http://www.wisdom.weizmann.ac.il/~eitany/hicpipe/>`_. This format consists of a pair of fend indices and a count of observed occurrences of that pairing. These indices must match those associated with the Fend object used when loading the data. Thus it is wise when using this format to also create the Fend object from a HiCPipe-style fend file to ensure accurate fend-count association.

::

  fend1    fend2    count
  1        4        10
  1        10       5
  1        13       1

.. note::
    In order to maintain compatibility with HiCPipe, both tabular fend files and MAT files are 1-indexed, rather than the standard 0-indexed used everywhere else with HiFive.
