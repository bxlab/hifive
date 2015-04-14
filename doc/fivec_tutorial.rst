.. _5C_tutorial:


*****************************
Tutorial for 5C Classes
*****************************

In order to perform analysis using the classes defined in :mod:`HiFive`, you must construct a series of objects, each one relying on the previous one. This tutorial will walk you through the process for doing a basic analysis of 5C data directly using the library of classes and functions.

.. _creating_a_fragment_object:

Creating a :class:`Fragment <hifive.fragment.Fragment>`  object
=================================================================

In order to use :mod:`HiFive` to process 5C data, the first thing that you need to do is create a :class:`Fragment <hifive.fragment.Fragment>` object. The :class:`Fragment <hifive.fragment.Fragment>`  object is an set of information contained in an HDF5 file. The object contains information about fragment sizes as determined by a series of boundary coordinates and orientations given in a BED formatted file, as well as an index for converting between string-based region names and numerical identifiers and any addition fragment features included in that bed file. The user may also include information about the genome and restriction enzyme, although these are optional and do not affect the functionality.

To create the :class:`Fragment <hifive.fragment.Fragment>`  object, you must provide a bed file in which each line contains the up- and downstream boundaries of a fragment and the strand to which this fragment belongs. Fragment names are also stored from the BED file and are used to connect mapped reads to fragments. Finally, groups of fragments can be separated into distinct regions which will be handled separately in downstream normalization and analysis. Fragments occurring on different chromosomes are always assigned to different regions. In addition, the user has the option of either explicitly specifying a list of regions or giving a minimum separation distance between fragments that is used to partition fragments into distinct regions.

To create a basic :class:`Fragment <hifive.fragment.Fragment>`  object, use the following command::

  import hifive
  fragment = hifive.Fragment(out_filename, mode='w')
  fragment.load_fragments(bed_filename, genome_name='MM9', re_name='HindIII')
  fragment.save()

In this case, the 'out_filename' specifies the location to save the :class:`Fragment <hifive.fragment.Fragment>`  object to. The 'bed_filename' contains the fragment boundaries, primer names, and strand information in BED format. The 'genome_name' and 're_name' are option strings that may be passed. In this case, the regions would be automatically determined using the default minimum distance between regions of 1 Mb. This could set by specifying a value with the keyword 'minregionspacing'.

In order to specify how to partition fragments into regions, we could use the following commands::

  import hifive
  fragment = hifive.Fragment(out_filename, mode='w')
  fragment.load_fragments(bed_filename, regions=[['chr1', start1, stop1], ['chr2', start2, stop2]])
  fragment.save()

This would create two regions that would include all of the fragments between coordinate values start1 - stop1 and start2 - stop2. When regions are explicitly given, the 'minregionspacing' variable is ignored.

.. note:
  The :class:`Fragment <hifive.fragment.Fragment>`  object can now be used by any experiment that relies on the same set of probes and does not need to be created separately for different experiments or analyses.

.. _creating_a_5C_dataset:

Creating a :class:`FiveCData <hifive.fivec_data.FiveCData>` object
===================================================================

In order to create a 5C dataset, you first need to have created an appropriate :class:`Fragment <hifive.fragment.Fragment>`  object. You can create the :class:`FiveCData <hifive.fivec_data.FiveCData>` object either from a file containing primer pairs and their observed counts or directly from mapped data in BAM format. If using counts data to create the dataset, the format should be::

  primer1   primer2 #_observed

where values are separated by tabs and '#_observed' is an integer. With either format the primer names must exactly match those in the BED file used to create the :class:`Fragment <hifive.fragment.Fragment>`  object.

To create the 5C dataset, you can run the following commands::

  import hifive
  data = hifive.FiveCData(out_filename, mode='w')
  data.load_data_from_counts(fragment_filename, [counts1.txt, counts2.txt])
  data.save()

In this case, 'out_filename' specifies the location to save the :class:`FiveCData <hifive.fivec_data.FiveCData>` object to. The 'fragment_filename' value is the location of the appropriate :class:`Fragment <hifive.fragment.Fragment>`  object. Multiple files containing counts data may be passed to the function as a list or, if only a single counts file is needed, it may be passed as a string. In order to load data from a set of BAM files, a similar procedure is used::

  import hifive
  data = hifive.FiveCDataset(out_filename, mode='w')
  data.load_data_from_bam(fragment_filename, [[bam_file1, bam_file2]])
  data.save()

In this case, the only difference is that pairs of file names corresponding to the two mapped read ends are passed as lists. Like the function for counts data, if only a single pair of files is needed, it may be passed as a list (not nested).

  Note: The :class:`FiveCData <hifive.fivec_data.FiveCData>` object can now be used by multiple analyses of this sample and does not need to be created separately for each one.

.. _creating_a_5C_analysis_object:

Creating a :class:`FiveC <hifive.fivec.FiveC>` project object
================================================================

The 5C project object, :class:`FiveC <hifive.fivec.FiveC>`, contains links to a :class:`FiveCData <hifive.fivec_data.FiveCData>` and :class:`Fragment <hifive.fragment.Fragment>`  object, information about which fragments to include in the analysis, model parameters, and learned model values. This is the standard way of working with 5C data in HiFive and this object will be used for learning the model, extracting portions of data, plotting, and downstream analysis.

To create a :class:`FiveC <hifive.fivec.FiveC>` object, you can use the following commands::

  import hifive
  fivec = hifive.FiveC(out_filename, 'w')
  fivec.load_data(data_filename)
  fivec.save()

In this case, 'out_filename' specifies the location to save the :class:`FiveC <hifive.fivec.FiveC>` object to. The 'data_filename' value is the location of the appropriate :class:`FiveCData <hifive.fivec_data.FiveCData>` object.

.. warning:: Because data and fragment data are stored in their own objects, each object keeps track of the location of its dependents through relative file names. This means that links between them will break if the relative pathway is changed.

.. _filter_5C_fragments:

Filter 5C fragments
=====================

Prior to modeling the data, you need to filter out fragments that have few valid reads mapped to them. :mod:`HiFive` uses an iterative filtering approach such that only when all fragments satisfy a user-defined minimum number of valid interactions does the filtering process cease.

To filter fragments, you can use the following commands::

  import hifive
  fivec = hifive.FiveC(fivec_filename)   
  fivec.filter_fragments(mininteractions=25)
  fivec.save()

In this case, 'fivec_filename' is a previously saved :class:`FiveC <hifive.fivec.FiveC>` analysis object. No value was passed to mode, since it defaults to 'r' for read. This loads the data from a previously created :class:`FiveCData <hifive.fivec_data.FiveCData>` object. In order for changes to be kept to a :class:`FiveC <hifive.fivec.FiveC>` object, it must be written to file using the save command.

.. _find_5C_distance_function:

Find 5C distance function
============================

:mod:`HiFive` approximates the distance-signal relationship using a power-law regression such that the log of the distance between the midpoints of two fragments and the log of their observed interactions. To do an initial estimate of this function, you can use the following command::

 fivec.find_distance_parameters()

.. _learn_5C_normalization_parameters:

Learn 5C normalization parameters
=================================

Using the probability algorithm
+++++++++++++++++++++++++++++++

In order to learn the correction model for 5C data using the probability algorithm, :mod:`HiFive` uses a bactracking line gradient descent. In addition, :mod:`HiFive` limits which interactions it uses to learn the model parameters to those that fall within a user-specified maximum interaction distance.

To learn 5C corrections using the probability approach, you can use the following command::

  fivec.find_probability_fragment_corrections(mindistance=50000,
                                              maxdistance=0,
                                              learningstep=0.5,
                                              max_iterations=1000,
                                              minchange=0.005,
                                              regions=[0, 1, 2])

In the above call, 'mindistance' is set to 50 kb, indicating that interactions shorter than that distance are no useed in calculating correction values. maxdistance' is set to zero, indicating that there is no upper limit on interaction distance to be used for learning model parameters. TThe 'learningstep' specifies how quickly to scale down the step value if the current try doesn't meet the arjimo learning criterion. The 'max_iterations' specifies a limit for how long to run the learning process for. 'minchange' is the stopping threshold such that if all absolute gradient values are below this the learning terminates early. Finally, the 'regions' parameter specifies that we only want to learn corrections for regions 0 - 3. Not specifying a value for this parameter would default to including all regions.


Using the express algorithm
+++++++++++++++++++++++++++++++

:mod:`HiFive` also offers a matrix-balancing approach for learning correction values. The primary differences to the probability model from the user's perspective are a single learning phase and a lack of learning rate.

To learn 5C corrections using the approximation approach, you can use the following command::

  fivec.find_express_fragment_corrections(iterations=1000,
                                          mindistance=50000,
                                          maxdistance=0,
                                          remove_distance=True)

In the above call, the 'remove_distance' argument specifies whether to remove the distance-dependent portion of the signal prior to approximating correction values. For best results, this should set to true (its default value).


Using the binning algorithm
+++++++++++++++++++++++++++++++

:mod:`HiFive` also offers a fragment characteristic-based approach adapted from the learning model used by `HiCPipe <http://www.wisdom.weizmann.ac.il/~eitany/hicpipe/>`_. This algorithm takes a list of features to be partitioned and a number of bins to partition them into and learns correction values associated with each partition based on a log-normal distribution of non-zero interactions corrected for distance-dependence.

To learn 5C corrections using the binning approach, you can use the following command::

  fivec.find_binning_fragment_corrections(max_iterations=1000,
                                          mindistance=50000,
                                          maxdistance=0,
                                          num_bins=[10, 10],
                                          model=['len', 'gc'],
                                          parameters=['even', 'fixed-const'],
                                          usereads='cis',
                                          learning_threshold=1.0)

Unlike the other two learning algorithms, this approach caps the learning iterations using 'max_iterations' and provides a means of early termination. This is done with the 'learning_threhold' parameter, which specifies that if the change in log-likelihood drops below 1.0, then cease iterating. The 'model', 'num_bins', and 'parameters' values should consist of equal-length lists and describe the correction values that are to be learned. Here, we told HiFive to use the length and gc content (specified in our BED file) for each fragment. Each feature was partitioned into a number of bins specified in 'num_bins'. The partitioning of length was done to create bins containing equal numbers of fragments while the gc content was divided such that each bin spanned an equal portion of the characteristic's range. Finally, the '-const' suffix told HiFive not to optimize the values for gc content. The 'usereads' value 'cis' specified that only within-region interactions should be used to learn these correction values.

Chaining learning algorithms
++++++++++++++++++++++++++++++

Because they work in very different ways, :mod:`HiFive` allows the binning algorithm to be chained with either the probability or express algorithm. The learning algorithms can be run in either order and the corrections from the first algorithm are applied prior to learning corrections for the second algorithm. This can be done by using the 'precorrect' option as follows::

  fivec.find_express_fragment_corrections(iterations=1000,
                                          mindistance=50000,
                                          maxdistance=0,
                                          remove_distance=True)
  fivec.find_binning_fragment_corrections(max_iterations=1000,
                                          mindistance=50000,
                                          maxdistance=0,
                                          num_bins=[10],
                                          model=['len'],
                                          parameters=['even'],
                                          usereads='cis',
                                          learning_threshold=1.0,
                                          precorrect=True)

.. _generating_a_fivec_heatmap:

Generating a heatmap
====================

In order to immediately make use of data, :mod:`HiFive` allows you to pull data from a region and create a heatmap. The data can be returned unbinned, binned using a fixed-width bin size, or binned using boundaries passed by the user. There are  several options for the format the data can be passed back in. Please refer to the :meth:`cis_heatmap <hifive.fivec.FiveC.cis_heatmap>` function for more details. There are also several options for transformations to the data. These are used to remove the distance-dependence signal, fragment bias, both, or to return only the predicted signal. In this example, we'll get a set of data from an entire region binned into 10 Kb bins as follows::

  heatmap = fivec.cis_heatmap(region=1,
                              binsize=10000,
                              arraytype='compact',
                              datatype='enrichment')

In the above call, 'enrichment' specifies to find the observed counts and expected counts, which includes the distance-dependence and fragment bias values. The observed counts are in the first index of the last dimension of the returned array, the expected counts are in the second index of the last dimension. 'compact' specifies a rectangular array where the first axis is the forward primers and the second axis is the reverse primers. 'region' refers to the region index given by :mod:`HiFive`. To find out more details about that region, we could do the following::

  print fivec.frags['regions'][1]

This returns the region's index, chromosome, starting fragment, stopping fragment (first fragment outside the region), starting coordinate and stopping coordinate.

Accessing heatmap data
======================

When a heatmap is generated, data are stored in an HDF5 dictionary, a binary file format that allows easy access through python. In order to access data from your heatmap, you can load it as follows::

  import h5py
  import numpy
  print heatmap.keys()
  heatmap = h5py.File(heatmap_file, 'r')
  counts = heatmap['0.counts'][...]
  expected = heatmap['0.expected'][...]
  enrichment = numpy.zeros((counts.shape[0], 2), dtype=numpy.float32)
  where = numpy.where(counts > 0)[0]
  enrichment[where, 0] = numpy.log(counts[where] / expected[where])
  enrichment[where, 1] = 1

Note that we used the 'r' option when opening the file with h5py. This ensures that we are in 'read' mode. You could also use 'a' for 'append' mode, which is the default. First we printed out the available dataset names in our heatmap file. These are all of the arrays that are accessible to us by calling them like any other key value in a dictionary. Next, in order to load data from the heatmap into memory rather than access it from the disk every time we refer to it, we use the '[...]' indexing call after pass the heatmap filestream the name of the data we want. In this case, we asked for the counts and expected values for region zero. In order to look at the enrichments, we took the log of the ratio of observed to expected values for each bin. However, there are likely bins that contain no observed counts which would give us a divide by zero error in the log function. So, we can use numpy's 'where' function to get a index list of places that match our criterion, in this case non-zero counts. Finally, we have made the enrichment array 2D so we can keep track of which bins are valid (nonzero counts). If we were looking at trans data, we would need one more dimension as the counts and expected arrays would be two-dimensional instead of one.

.. _plotting_a_fivec_heatmap:

Plotting a heatmap
==================

In order to visualize the heatmap we just produced, :mod:`HiFive` has several plotting functions that take different shaped arrays. The function called needs to match the array produced. However, in this case, the 5C compact array is compatible with the :meth:`plot_full_array <hifive.plotting.plot_full_array>` function, so we'll use that as follows::

  img = hifive.plotting.plot_full_array(heatmap, symmetric_scaling=True)
  img.save(out_fname)

In calling the function, we pass the heatmap and that would be sufficient. There are, however, additional options. For example, 'symmetric_scaling' specifies whether the color scale should be expanded to run from the minimum value to the maximum (False) or so that the maximum absolute value determine both upper and lower color bounds. The image returned is a PIL image of type 'png'.
