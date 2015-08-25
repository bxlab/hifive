.. _multiresolution heatmap:

*****************************
Multi-Resolution Heatmap
*****************************

HiFive has a unique filetype called a 'multi-resolution heatmap' or MRH, which is a binary file holding information necessary for plotting any part of the genome across a wide range of resolutions using intelligent binning based on observation density. Files are indexed in a recursive fashion making retrieval incredibly fast. This makes exploring HiC data incredibly easy, either through the included plotting script or when combined with the interactive HiC browser for multi-resolution heatmaps support in `Galaxy <https://usegalaxy.org/>`_.

===============
MRH Approach
===============

The idea behind the multi-resolution heatmap is to have a file that is as small as possible while containing data  in the highest resolution reasonable given the density of information present. Each chromosome (or chromosome pair for inter-chromosomal interactions) is binned at a low resolution. Each bin is designated as valid or invalid based on an observation cutoff (the number of reads that must be seen in a given bin). Each valid bin is then subdivided into four equal-sized bins. This respresents the next level of resolution. These new bins are also subjected to the same read cutoff. Lower resolution bins also have associated with them an index specifying where their subdivided bins are located in the file and a shape indicator specifying how many and which sub-bins are valid. This has two advantages. First, there are no invalid bins included in the file past the lowest level of resolution, and second, retrieving data from any pair of coordinate ranges is easy. This is because to locate any 2D position in the heatmap, one simple identifies the correct low-resolution bin and then recursively selects the appropriate subdivision. The final advantage is that for positions that have insufficient data at a higher level of resolution, the value of a lower resolution bin containing the target position can be returned. This makes visual exploration much easier as there are effectively no missing values.

.. _mrh_program:

=======================
Standalone MRH program
=======================

HiFive includes a command line program, 'fetch_mrh_data', that allows data to be retrieved from an MRH file and either plotted as an image or written to a genomic interval text file. 

::

  > fetch_mrh_data [-h] -c CHROM [-C CHROM2] [-s START] [-S START2] [-e END]
  		[-E END2] [-R MAXRES] [-r MINRES] [-t] heatmap output

Arguments:

:project:  The HiFive MRH file to create pull data from.
:output: A filename to write the image or interval file to.

Options:

-h, --help            Display the help message, options, and arguments and exit.
-c, --chrom str       The first region chromosome to pull data from.
-C, --chrom2 str      The second region chromosome to pull data from. If no value is passed, this will be set to the same value as 'chrom'.
-s, --start int       The first coordinate of the chromosome to pull data from. If no value is passed, this will be set to the first bin position in the MRH file for 'chrom'. [None]
-S, --start2 int       The first coordinate of the chromosome to pull data from. If no value is passed, this will be set to the first bin position in the MRH file for 'chrom2'. [None]
-e, --end int          The last coordinate + 1 of the chromosome to pull data from. If no value is passed, this will be set to the last bin position in the MRH file for 'chrom'. [None]
-E, --end2 int         The last coordinate + 1 of the chromosome to pull data from. If no value is passed, this will be set to the last bin position in the MRH file for 'chrom2'. [None]
-R, --max-resolution   The maximum resolution bound for returned data. If no value is passed, this will be set to the highest resolution available in the heatmap for the chromosome(s). [default: None]
-r, --min-resolution   The minimum resolution bound for returned data. If no value is passed, this will be set to the lowest resolution available in the heatmap for the chromosome(s). [default: None]
-t, --text             Write a genomic interval text file instead of an image.

==========================
File specs
==========================

The MRH file is composed of a header and series of inline 1D arrays of data. The header contains information about what chromosomes are present, if inter-chromosomal interactions are present, and coordinate and score information for each chromosome/chromosome-pair. For each chromosome and chromosome-pair, there are three arrays of information. the First contains the enrichment data, the second contains bin index values, and the third contains information about the validity of bins for each resolution level. All values in the file are big-endian.

----------------
The File Header
----------------

The header contains the following information:

:magic number:                           A filetype identifier of 4 bytes, '0x42054205'
:interchromosomal flag:                  A 32-bit integer of 1 or 0, indicating if inter-chromosomal interactions are included or not.
:number of chromosomes:                  A 32-bit integer specifying how many chromosomes are in the files index.
:chromosome name sizes:                  An N * 32-bit integer array indicating the length in characters of each of N chromosome names.
:chromosome names:                       A sum(name sizes) * 1-byte character array containing each chromosome's name.
:chromosome bit index:                   An index of the starting bit for each chromosome/chromosome-pair. If inter-chromosomal data is not included, this is an N * 32-bit integer array. Otherwise this is a N * (N - 1) / 2 * 32-bit integer array with chromosomes ordered 1x1, 1x2... 1xN, 2x2, 2x3... NxN. All subsequent arrays for chromosomes/chromosome-pairs are ordered the same way. 
:chromosome partitions:                  An N * 32-bit integer array specifying the number of genomic paritions (bins along 1 axis) each chromosome is partitioned into.
:trans chromosome partitions:            If trans interactions are included, this is an N * 32-bit integer array specifying the number of genomic paritions (bins along 1 axis) each chromosome is partitioned into for inter-chromosomal interactions.
:data array sizes:                       A # chroms/pairs * 32-bit integer array specifying how many 32-bit floating point values each chromosome/pair's data array contains.
:index array sizes:                      A # chroms/pairs * 32-bit integer array specifying how many 32-bit integer values each chromosome/pair's index array contains.
:chromosome starting coordinates:        An N * 32-bit integer array specifying the first coordinate of each chromosome's smallest intra-chromosomal bin.
:trans chromosome starting coordinates:  If trans interactions are included, an N * 32-bit integer array specifying the first coordinate of each chromosome's smallest inter-chromosomal bin.
:chromosome stopping coordinates:        An N * 32-bit integer array specifying the last coordinate of each chromosome's largest intra-chromosomal bin.
:trans chromosome stopping coordinates:  If trans interactions are included, an N * 32-bit integer array specifying the last coordinate of each chromosome's largest inter-chromosomal bin.
:smallest enrichment scores:             A # chroms/pairs * 32-bit floating point array of values with the smallest enrichment value for each chromosome/pair.
:largest enrichment scores:              A # chroms/pairs * 32-bit floating point array of values with the largest enrichment value for each chromosome/pair.
:maximum bin size:                       A 32-bit integer specifying the largest bin size (lowest resolution) for intra-chromosomal interactions.
:maximum trans bin size:                 If trans interactions are included, a 32-bit integer specifying the largest bin size (lowest resolution) for inter-chromosomal interactions.
:minimum bin size:                       A 32-bit integer specifying the smallest bin size (highest resolution) for intra-chromosomal interactions.
:minimum trans bin size:                 If trans interactions are included, a 32-bit integer specifying the smallest bin size (highest resolution) for inter-chromosomal interactions.
:minimum observation cutoff:             A 32-bit integer specifying the minimum number of observations needed for a bin to be included in the MRH file.

-------------
Data Arrays
-------------

Data for each chromosome and chromosome-pair consists of three arrays, the interaction, index, and shape arrays. The first contains all of the actual enrichment values. The interaction array contains data as 32-bit floating point values grouped by resolution level, going from lowest (largest bin size) to highest. The lowest resolution data contains all possible bins within the coordinate bounds, either a flattened rectangular array (row-major) for inter-chromosomal data or a flattened upper-triangle array for intra-chromosomal data. Subsequent resolution data is presented in groups of 4 or fewer bins.

For each interaction bin, excepting the lowest resolution bins, there is a corresponding index bin and shape bin in the same order. The index bin specifies the index in the data array containing the next level of resolution up for data within the current bin. If no valid higher resolution bins are contained within the current bin's boundaries, then the index bin has a value of -1. The shape bin contains a 16-bit integer that uses the first 4 bits to indicate which of the four possible bins contain valid data in a row-major order. For example a value of 9 would indicate the top left and bottom right corners of the 2 by 2 grid would have valid values. This is also used to indicate how many bins at that higgher resolution correspond to the current bin.
