.. _5C_tutorial:


********************
A Brief 5C Tutorial
********************

In order to perform analysis with HiFive, you must construct a series of objects, each one relying on the previous one. This tutorial will walk you through the process for doing a basic analysis of 5C data.

.. _creating_a_fragment_object:

Creating a Fragment object
================================

In order to use HiFive to process 5C data, the first thing that you need to do is create a Fragment object. The Fragment object is an set of information contained in an HDF5 file. The object contains information about fragment sizes as determined by a series of boundary coordinates and orientations given in a BED formatted file, as well as an index for converting between string-based region names and numerical identifiers. The user may also include information about the genome and restriction enzyme, although these are optional and do not affect the functionality.

To create the Fragment object, you must provide a bed file in which each line contains the up- and downstream boundaries of a fragment and the strand to which this fragment belongs. Fragment names are also stored from the BED file but are not required. Finally, groups of fragments can be separated into distinct regions which will be handled separately in downstream normalization and analysis. Fragments occuring on different chromosomes are always assigned to different regions. In addition, the user has the option of either explicitly specifying a list of regions or giving a minimum separation distance between fragments that is used to partition fragments into distinct regions.

To create a basic Fragment object, use the following command::

  import hifive
  fragment = hifive.Fragment(out_filename, mode='w')
  fragment.load_fragments(bed_filename, genome_name='MM9', re_name='HindIII')
  fragment.save()

In this case, the 'out_filename' specifies the location to save the Fragment object to and should end with the '.hdf5' extension. The 'bed_filename' contains the fragment boundaries, primer names, and strand information in BED format. The 'genome_name' and 're_name' are option strings that may be passed. In this case, the regions would be automatically determined using the default minimum distance between regions of 1 Mb. This could set by specifying a value with the keyword 'minregionspacing'.

In order to specifiy how to partion fragments into regions, we could use the following commands::

  import hifive
  fragment = hifive.Fragment(out_filename, mode='w')
  fragment.load_fragments(bed_filename, regions=[['chr1', start1, stop1], ['chr2', start2, stop2]])
  fragment.save()

This would create two regions that would include all of the fragments between coordinate values start1 - stop1 and start2 - stop2. When regions are explicitly given, the 'minregionspacing' variable is ignored.

.. note:
  The Fragment object can now be used by any experiment that relies on the same set of probes and does not need to be created separately for different experiments or analyses.

.. _creating_a_5C_dataset:

Creating a 5C dataset
================================

In order to create a 5C dataset, you first need to have created an appropriate Fragment object. You can create the dataset object either from a file containing primer pairs and their observed counts or directly from mapped data in BAM format. If using counts data to create the dataset, the format should be::

  primer1   primer2 #_observed

where values are separated by tabs and '#_observed' is an integer. With either format the primer names must exacly match those in the BED file used to create the Fragment object.

To create the 5C dataset, you can run the following commands::

  import hifive
  data = hifive.FiveCDataset(out_filename, mode='w')
  data.load_data_from_counts(fragment_filename, [counts1.txt, counts2.txt])
  data.save()

In this case, 'out_filename' specifies the location to save the FiveCDataset object to and should end with the '.hdf5' extension. The 'fragment_filename' value is the location of the appropriate Fragment object. Multiple files containing counts data may be passed to the function as a list or, if only a single counts file is needed, it may be passed as a string. In order to load data from a set of BAM files, a similar procedure is used::

  import hifive
  data = hifive.FiveCDataset(out_filename, mode='w')
  data.load_data_from_bam(fragment_filename, [bam_prefix1, bam_prefix2])
  data.close()

In this case, the only difference is that prefices are passed instead of complete file names. The prefices should be the entire file path up until the strand specifier such that the two file names are created by appending either '_1.bam' or '_2.bam' to the prefix. Like the function for counts data, if only a single prefix is needed it may be passed as a string.

  Note: The FiveC data object can now be used by multiple analyses of this sample and does not need to be created separately for each one.

.. _creating_a_5C_analysis_object:

Creating a 5C analysis object
================================

The 5C analysis object, FiveC, contains links to a FiveCData and Fragment object, information about which fragments to include in the analysis, model parameters, and learned model values. This is the standard way of working with 5C data in HiFive and this object will be used for learning the model, extracting portions of data, plotting, and downstream analysis.

To create a FiveC object, you can use the following commands::

  import hifive
  fivec = hifive.FiveC(out_filename, 'w')
  fivec.load_data(data_filename)
  fivec.save()

In this case, 'out_filename' specifies the location to save the FiveC object to and should end with the '.hdf5' extension. The 'data_filename' value is the location of the appropriate data object.

.. warning:: Becauase data and fragment data are stored in their own objects, each object keeps track of the location of its dependents through relative file names. This means that links between them will break if the relative pathway is changed.

.. _filter_5C_fragments:

Filter 5C fragments
=====================

Prior to modeling the data, you need to filter out fragments that have few valid reads mapped to them. HiFive uses an iterative filtering approach such that only when all fragments satisfy a user-defined minimum number of valid interactions does the filtering process cease.

To filter fragments, you can use the following commands::

  import hifive
  fivec = hifive.FiveC(fivec_filename)   
  fivec.filter_fragments(mininteractions=25)
  fivec.save()

In this case, 'fivec_filename' is a previously saved FiveC analysis object. No value was passed to mode, since it defaults to 'r' for read. This loads the data from a previously created FiveCData object. In order for changes to be kept to a FiveC object, it must be written to file using the save command.

.. _find_5C_distance_function:

Find 5C distance function
============================

HiFive approximates the distance-signal relationship using a power-law regression such that the log of the distance between the midpoints of two fragments and the log of their observed interactions. To do an initial estimate of this function, you can use the following command::

 fivec.find_distance_parameters()

.. _learn_5C_normalization_parameters:

Learn 5C normalization parameters
=================================

In order to learn the correction model for 5C data, HiFive uses two rounds of gradient descent, one with constant learning rate (the 'burn-in' phase) and the second with a linearly decreasing learning rate (the 'annealing' phase). In addition, HiFive can recalculate the distance function parameters periodically using the correction-adjusted interaction values. Finally, HiFive limits which interactions it uses to learn the model parameters to those that fall within a user-specified maximum interaction distance.

To learn 5C corrections using the modeling approach, you can use the following command::

  fivec.find_fragment_corrections(display=100,
                                  maxdistance=0,
                                  learningrate=0.01,
                                  burnin_iterations=5000,
                                  annealing_iterations=10000,
                                  recalculate_distance=100)

In the above call, 'maxdistance' is set to zero, indicating that there is no upper limit on interaction distance to be used for learning model parameters. The 'recalculate_distance' parameters specifies how many iterations to wait before recalculating the distance parameters. The 'learningrate' specifies what percentage of the gradient to apply towards value updates. One last value passed to the function in 'display', which specifies how many iterations should pass before updating the display (via STDERR). This can also be set to zero to not display the progress.

.. _approximate_5C_normalization_parameters:

Approximate 5C normalization parameters
=======================================

HiFive also offers an approximation approach for learning correction values. The primary differences to the correction model from the user's perspective are a single learning phase and a lack of learning rate. The approximation learning approach can still recalculate the distance function parameters periodically.

To learn 5C corrections using the approximation approach, you can use the following command::

  fivec.find_fragment_corrections(iterations=1000,
                                  recalculate_distance=100,
                                  remove_distance=True)

In the above call, the 'remove_distance' argument specifies whether to remove the distance-dependent portion of the signal prior to approximating correction values. For best results, this should set to true (its default value).

.. _generating_a_fivec_heatmap:

Generating a heatmap
====================

In order to immediately make use of data, HiFive allows you to pull data from a regions and create a heatmap. The data can be returned unbinned, binned using a fixed-width bin size, or binned using boundaries passed by the user. There are  several options for the format the data can be passed back in. Please refer to the :meth:`hifive.fivec_binning.bin_cis_signal` function for more details. There are also several options for transformations to the data. These are used to remove the distance-dependence signal, fragment bias, both, or to return only the predicted signal. In this example, we'll get a set of data from an entire region binned into 10 Kb bins as follows::

  heatmap = hifive.fivec_binning(fivec,
                                 region=1,
                                 binsize=10000,
                                 arraytype='compact',
                                 datatype='enrichment')

In the above call, 'enrichment' specifies to find the observed counts and expected counts, which includes the distance-dependence and fragment bias values. The observed counts are in the first index of the last dimension of the returned array, the expected counts are in the second index of the last dimension. 'compact' specifies a rectangular array where the first axis is the forward primers and the second axis is the reverse primers. 'region' refers to the region index given by HiFive. To find out more details about that region, we could do the following::

  fivec.frags['regions'][1]

This returns the region's chromosome, starting fragment, stopping fragment (first fragment outside the region), starting coordinate and stopping coordinate.

.. _plotting_a_fivec_heatmap:

Plotting a heatmap
==================

In order to visualize the heatmap we just produced, HiFive has several plotting functions that take different shaped arrays. The function called needs to match the array produced. However, in this case, the 5C compact array is compatible with the :meth:`hifive.plotting.plot_full_array` function, so we'll use that as follows::

  img = hifive.plotting.plot_full_array(heatmap, symmetric_scaling=True)
  img.save(out_fname)

In calling the function, we pass the heatmap and that would be sufficient. There are, however, additional options. For example, 'symmetric_scaling' specifies whether the color scale should be expanded to run from the minimum value to the maximum (False) or so that the maximum absolute value determine both upper and lower color bounds. The image returnd is a PIL image of type 'png'.

.. note:: The next thing on the todo list is write wrappers within the :class:`FiveC` and :class:`HiC` classes for running binning and plotting through the analysis objects themselves.
