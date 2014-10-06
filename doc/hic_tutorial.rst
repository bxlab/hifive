.. _HiC_tutorial:


********************
A Brief HiC Tutorial
********************

In order to perform analysis with :mod:`HiFive`, you must construct a series of objects, each one relying on the previous one. This tutorial will walk you through the process for doing a basic analysis of HiC data.

.. _creating_a_fend_object:

Creating a :class:`Fend <hifive.fend.Fend>` object
===================================================

In order to use :mod:`HiFive` to process HiC data, the first thing that you need to do is create a fragment-end (Fend) object. The :class:`Fend <hifive.fend.Fend>` object is an set of information contained in an HDF5 file. The object contains information about the fragments created by digestion of a genome by a specific restriction enzyme (RE) as well as an index for converting between string-based chromosome names and numerical identifiers. The user may also include information about the genome and restriction enzyme, although these are optional and do not affect the functionality.

To create a :class:`Fend <hifive.fend.Fend>` object, you must provide the location of RE fragments. This information is supplied in the form of either a 'fend' file with a format compatible with `HiCPipe <http://www.wisdom.weizmann.ac.il/~eitany/hicpipe/>`_ or a BED-formatted file containing fragment boundaries produced by RE digest or locations of the RE cut sites.

To create a basic :class:`Fend <hifive.fend.Fend>` object, use the following command::

  import hifive
  fend = hifive.Fend(out_filename, mode='w')
  fend.load_fends(RE_data_filename, genome_name='MM9', re_name='HindIII')
  fend.save()

In this case, the 'out_filename' specifies the location to save the :class:`Fend <hifive.fend.Fend>` object to and should end with the '.hdf5' extension. The 'RE_data_filename' contains the RE fragment boundaries in one of the formats described above. The 'genome_name' and 're_name' are option strings that may be passed.

.. note:: The :class:`Fend <hifive.fend.Fend>` object can now be used by any experiment that relies on the same genome / restriction enzyme combination and does not need to be created separately for different experiments or analyses.

.. _creating_a_HiC_dataset:

Creating a :class:`HiCData <hifive.hic_data.>` object
======================================================

In order to create a HiC dataset, you first need to have created an appropriate :class:`Fend <hifive.fend.Fend>` object. You can create the :class:`HiCData <hifive.hic_data.>` object from a file containing fend indices and their observed counts, the format used with `HiCPipe <http://www.wisdom.weizmann.ac.il/~eitany/hicpipe/>`_ data (the 'mat' format), a text file containing pairs of chromosomes, coordinates, and strands (referred to from now on as a 'raw' file), or directly from mapped data in BAM format. If using raw coordinate data to create the dataset, the format should be::

  chr1  coord1   strand1    chr2    coord2  strand2

where values are separated by tabs and strands are denoted by '+' or '-'. In addition to mapped reads, you need to provide a maximum insert size, i.e. the total acceptable length of a sequenced fragment as determined by the sum of each mapping location to its downstream RE cutsite.

To create the 5C dataset from a HiCPipe-compatible format, you can run the following commands::

  import hifive
  data = hifive.HiCData(out_filename, mode='w')
  data.load_data_from_mat(fend_filename, data_filename)
  data.save()

In this case, 'out_filename' specifies the location to save the :class:`HiCData <hifive.hic_data.HiCData>` object to and should end with the '.hdf5' extension. The 'fend_filename' value is the location of the appropriate :class:`Fend <hifive.fend.Fend>` object. To maintain compatibility with HiCPipe-formatted 'mat' files, :mod:`HiFive` expects that fend and fragment numbering begin at index 1, not 0.

To create the 5C dataset from raw coordinate data, you can run the following commands::

  import hifive
  data = hifive.HiCData(out_filename, mode='w')
  data.load_data_from_raw(fend_filename, [raw1.txt, raw2.txt])
  data.save()

In this case, 'out_filename' specifies the location to save the :class:`HiCData <hifive.hic_data.>` object to and should end with the '.hdf5' extension. The 'fend_filename' value is the location of the appropriate :class:`Fend <hifive.fend.Fend>` object. Multiple files containing paired-end mapped coordinates may be passed to the function as a list or, if only a single file is needed, it may be passed as a string.

In order to load data from a set of BAM files, a similar procedure is used::

  import hifive
  data = hifive.HiCData(out_filename, mode='w')
  data.load_data_from_bam(fragment_filename, [bam_prefix1, bam_prefix2])
  data.save()

In this case, the only difference is that prefices are passed instead of complete file names. The prefices should be the entire file path up until the strand specifier such that the two file names are created by appending either '_1.bam*' or '_2.bam*' to the prefix. The asterix denotes that all files with that name but different suffices are loaded. This allows files that are generated by multiple rounds of alignment that differ only by a suffix can be loaded without additional manipulation. Like the function for raw data, if only a single prefix is needed it may be passed as a string.

.. note:: The :class:`HiCData <hifive.hic_data.>` object can now be used by multiple analyses of this sample and does not need to be created separately for each one.

.. _creating_a_HiC_analysis_object:

Creating a :class:`HiC <hifive.hic.HiC>` analysis object
=========================================================

The HiC analysis object, :class:`HiC <hifive.hic.HiC>`, contains links to a HiCData and :class:`Fend <hifive.fend.Fend>` object, information about which fends to include in the analysis, model parameters, and learned model values. This is the standard way of working with HiC data in :mod:`HiFive` and this object will be used for learning the model, extracting portions of data, plotting, and downstream analysis.

To create a :class:`HiC <hifive.hic.HiC>` object, you can use the following commands::

  import hifive
  hic = hifive.HiC(out_filename, 'w')
  hic.load_data(data_filename)
  hic.save()

In this case, 'out_filename' specifies the location to save the :class:`HiC <hifive.hic.HiC>` object to and should end with the '.hdf5' extension. The 'data_filename' value is the location of the appropriate data object.

.. warning:: Becauase data and fragment data are stored in their own objects, each object keeps track of the location of its dependents through relative file names. This means that links between them will break if the relative pathway is changed.

.. _filter_HiC_fends:

Filter HiC fends
=====================

Prior to modeling the data, you need to filter out fends that have few valid reads mapped to them. :mod:`HiFive` uses an iterative filtering approach such that only when all fends satisfy a user-defined minimum number of valid interactions does the filtering process cease.

To filter fends, you can use the following commands::

  import hifive
  hic = hifive.HiC(hic_filename)   
  hic.filter_fends(mininteractions=25, maxdistance=5000000)
  hic.save()

In this case, 'hic_filename' is a previously saved :class:`HiC <hifive.hic.HiC>` analysis object. No value was passed to mode, since it defaults to 'r' for read. This loads the data from a previously created HiCData object. In order for changes to be kept to a FiveC object, it must be written to file using the save command. The 'maxdistance' argument specifies that only reads associated with interactions spanning that distance or less are counted for purposes of filtering fends.

.. _find_HiC_distance_function:

Find HiC distance function
============================

:mod:`HiFive` approximates the distance-signal relationship using a series of linear transitions between bin means of mean log interaction counts. Spanning from a user-defined minimum interaction distance up to the genome maximum interaction distance, the range is divided into equal-sized log distance bins. Values falling between bin midpoints are interpolated based on a linear transition between bins. To do an initial estimate of this function, you can use the following command::

  hic.find_distance_means(numbins=90,
                          minsize=200, 
                          maxsize=0,
                          smoothed=2,
                          corrected=False)

In this function call, the range of interaction sizes is being broken into 90 bins, 1 bin covering interactions <= 200 bp, and the other 89 spanning up to the maximum interaction distance with breaks evenly spaced in log space. The maximum of this range is set by 'maxsize', which can either be zero, as in this call, setting the maximum size equal to the longest interaction distance, or a positive integer value which would exclude any interaction distances greater than 'maxsize'. The 'smoothed' keyword specifies whether to apply a triangular smoothing to bin means after finding them. A value of zero would remain unsmoothed, one extends to using adjacent values, two would include values up to two steps away, and so on. The final argument, 'corrected', specifies whether to apply fend corrections to the interaction counts prior to calculating bin means.

Because this function involves scanning large amounts of data, it has been made to utilize MPI. To do so, you can use a scripts such as the following::

  import hifive
  from mpi4py import MPI

  rank = MPI.COMM_WORLD.Get_rank()
  hic = hifive.HiC(hic_filename)
  hic.find_distance_means(numbins=90,
                          minsize=200, 
                          maxsize=0,
                          smoothed=2,
                          corrected=False)
  if rank == 0:
    hic.save()

.. _learn_HiC_normalization_parameters:

Learn HiC normalization parameters
===================================

In order to learn the correction model for HiC data, :mod:`HiFive` uses two rounds of gradient descent, one with constant learning rate (the 'burn-in' phase) and the second with a linearly decreasing learning rate (the 'annealing' phase). In addition, :mod:`HiFive` can recalculate the distance function parameters periodically using the correction-adjusted interaction values. Finally, :mod:`HiFive` limits which interactions it uses to learn the model parameters to those that fall within a user-specified maximum interaction distance.

To learn HiC corrections using the modeling approach, you can use the following command::

  hic.find_fend_corrections(display=100,
                            maxdistance=5000000,
                            learningrate=0.01,
                            burnin_iterations=5000,
                            annealing_iterations=10000,
                            recalculate_distance=2500)

In the above call, 'maxdistance' indicates that interctions spanning more than 5 Mb are excluded from calculations. Setting this to zero would include all unfiltered cis interactions. The 'recalculate_distance' parameters specifies how many iterations to wait before recalculating the distance parameters. The 'learningrate' specifies what percentage of the gradient to apply towards value updates. One last value passed to the function in 'display', which specifies how many iterations should pass before updating the display (via STDERR). This can also be set to zero to not display the progress.

Because of the large numbers of calculations involved in this function, it has been made to utilize MPI. To do so, you can use a scripts such as the following::

  import hifive
  from mpi4py import MPI

  rank = MPI.COMM_WORLD.Get_rank()
  hic = hifive.HiC(hic_filename)
  hic.find_fend_corrections(display=100,
                            maxdistance=5000000,
                            learningrate=0.01,
                            burnin_iterations=5000,
                            annealing_iterations=10000,
                            recalculate_distance=2500)
  if rank == 0:
    hic.save()

.. _approximate_HiC_normalization_parameters:

Approximate HiC normalization parameters
=========================================

:mod:`HiFive` also offers an approximation approach for learning correction values. The primary differences to the correction model from the user's perspective are a single learning phase and a lack of learning rate. The approximation learning approach can still recalculate the distance function parameters periodically.

To learn HiC corrections using the approximation approach, you can use the following command::

  hic.find_express_fend_corrections(iterations=1000,
                                    mindistance=0,
                                    usereads='cis',
                                    remove_distance=True,
                                    recalculate_distance=200)

In the above call, 'mindistance' is used to exclude interaction distances shorter that the passed value. If this results in the exclusion of any reads, fends are refiltered using either the value passed under the keyword 'mininteractions' or, if that is not specified, the value passed the last time fends were filtered. The 'usereads' argument allows the user to base the correction value approximation on 'cis' interactions, 'trans' interactions, or 'all'. Selecting 'trans' interactions will also result in a refiltering of fends to ensure that all of them are involved in sufficient interactions as described previously. The 'remove_distance' argument specifies whether to remove the distance-dependent portion of the signal prior to approximating correction values. For best results, this should set to true (its default value).

Although this function is much more computationally efficient, the recalculation of the distance function can take time and so has been made to utilize the MPI environment when available as follows::

  import hifive
  from mpi4py import MPI

  rank = MPI.COMM_WORLD.Get_rank()
  hic = hifive.HiC(hic_filename)
  hic.find_express_fend_corrections(iterations=1000,
                                    mindistance=0,
                                    usereads='cis',
                                    remove_distance=True,
                                    recalculate_distance=200)
  if rank == 0:
    hic.save()

.. _generating_a_hic_heatmap:

Generating a heatmap
====================

In order to immediately make use of data, :mod:`HiFive` allows you to pull data from a region and create a heatmap. The data can be returned unbinned, binned using a fixed-width bin size, or binned using boundaries passed by the user. There are  several options for the format the data can be passed back in. Please refer to the :meth:`hifive.hic_binning.bin_cis_signal` function for more details. There are also several options for transformations to the data. These are used to remove the distance-dependence signal, fend bias, both, or to return only the predicted signal. In this example, we'll get a portion of chromosome 1 binned into 10 Kb bins as follows::

  heatmap = hifive.hic_binning(hic,
                              chrom='1',
                              start=1000000
                              stop=3000000
                              binsize=10000,
                              arraytype='upper',
                              datatype='enrichment')

In the above call, All valid possible interactions were queried from chromosome 1 between 1000000 and 3000000. For valid interactions that had no observation, an expected value was still added to the bin. 'enrichment' specifies to find the observed counts and expected counts, which includes the distance-dependence and fend bias values. The observed counts are in the first index of the last dimension of the returned array, the expected counts are in the second index of the last dimension. 'Upper' specifies a row-major upper triangle array (all values above the matrix diagonal flattened).

.. _plotting_a_hic_heatmap:

Plotting a heatmap
==================

In order to visualize the heatmap we just produced, :mod:`HiFive` has several plotting functions that take different shaped arrays. The function called needs to match the array produced. In this case, we produced an upper array which is compatible with the :meth:`hifive.plotting.plot_upper_array` function, so we'll use that as follows::

  img = hifive.plotting.plot_upper_array(heatmap, symmetric_scaling=True)
  img.save(out_fname)

In calling the function, we pass the heatmap and that would be sufficient. There are, however, additional options. For example, 'symmetric_scaling' specifies whether the color scale should be expanded to run from the minimum value to the maximum (False) or so that the maximum absolute value determine both upper and lower color bounds. The image returnd is a :mod:`PIL` image of type 'png'.

.. note:: The next thing on the todo list is write wrappers within the :class:`FiveC <hifive.fivec.FiveC>` and :class:`HiC <hifive.hic.HiC>` classes for running binning and plotting through the analysis objects themselves.