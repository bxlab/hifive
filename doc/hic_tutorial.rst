.. _HiC_tutorial:


*****************************
Tutorial for HiC Classes
*****************************

In order to perform analysis using the classes defined in :mod:`HiFive`, you must construct a series of objects, each one relying on the previous one. This tutorial will walk you through the process for doing a basic analysis of HiC data directly using the library of classes and functions.

.. _creating_a_fend_object:

Creating a :class:`Fend <hifive.fend.Fend>` object
===================================================

In order to use :mod:`HiFive` to process HiC data, the first thing that you need to do is create a fragment-end (Fend) object. The :class:`Fend <hifive.fend.Fend>` object is an set of information contained in an HDF5 file. The object contains information about the fragments created by digestion of a genome by a specific restriction enzyme (RE) as well as an index for converting between string-based chromosome names and numerical identifiers. The user may also include information about the genome and restriction enzyme, although these are optional and do not affect the functionality.

To create a :class:`Fend <hifive.fend.Fend>` object, you must provide the location of RE fragments. This information is supplied in the form of either a 'fend' file with a format compatible with `HiCPipe <http://www.wisdom.weizmann.ac.il/~eitany/hicpipe/>`_ or a BED-formatted file containing fragment boundaries produced by RE digest or locations of the RE cut sites.

To create a basic :class:`Fend <hifive.fend.Fend>` object, use the following command::

  import hifive
  fend = hifive.Fend(out_filename, mode='w')
  fend.load_fends(RE_data_filename, genome_name='MM9', re_name='HindIII', format='bed')
  fend.save()

In this case, the 'out_filename' specifies the location to save the :class:`Fend <hifive.fend.Fend>` object to. The 'RE_data_filename' contains the RE fragment boundaries in one of the formats described above. The 'genome_name' and 're_name' are option strings that may be passed. The 'format' argument specifies that we are passing a BED file containing the fend data to load.

To create a RE-based binned :class:`Fend <hifive.fend.Fend>` object, use the following command::

  import hifive
  fend = hifive.Fend(out_filename, mode='w', binned=40000)
  fend.load_fends(RE_data_filename, genome_name='MM9', re_name='HindIII', format='bed')
  fend.save()

This would create a 40Kb parition of the genome in addition to the RE-based parition. Now data associated with the Fend file would be filtered as usual for RE fragments but then assigned to the appropriate 40Kb bin instead of maintained as fend pairs. You could also skip the restriction enzyme parition information if it won't be needed as follows::

  import hifive
  fend = hifive.Fend(out_filename, mode='w', binned=40000)
  fend.load_bins(chrom_length_filename, genome_name='MM9', format='len')
  fend.save()

This creates a parition of uniform-sized bins starting at coordinate zero for each chromosome.

.. note:: The :class:`Fend <hifive.fend.Fend>` object can now be used by any experiment that relies on the same genome / restriction enzyme combination and does not need to be created separately for different experiments or analyses.

.. _creating_a_HiC_dataset:

Creating a :class:`HiCData <hifive.hic_data.>` object
======================================================

In order to create a HiC dataset, you first need to have created an appropriate :class:`Fend <hifive.fend.Fend>` object. You can create the :class:`HiCData <hifive.hic_data.>` object from a file containing fend indices and their observed counts, the format used with `HiCPipe <http://www.wisdom.weizmann.ac.il/~eitany/hicpipe/>`_ data (the 'mat' format), a text file containing pairs of chromosomes, coordinates, and strands (referred to from now on as a 'raw' file), or directly from mapped data in BAM format. If using raw coordinate data to create the dataset, the format should be::

  chr1  coord1   strand1    chr2    coord2  strand2

where values are separated by tabs and strands are denoted by '+' or '-'. In addition to mapped reads, you need to provide a maximum insert size, i.e. the total acceptable length of a sequenced fragment as determined by the sum of each mapping location to its downstream RE cutsite.

To create the HiC dataset from a HiCPipe-compatible format, you can run the following commands::

  import hifive
  data = hifive.HiCData(out_filename, mode='w')
  data.load_data_from_mat(fend_filename, data_filename)
  data.save()

In this case, 'out_filename' specifies the location to save the :class:`HiCData <hifive.hic_data.HiCData>` object to. The 'fend_filename' value is the location of the appropriate :class:`Fend <hifive.fend.Fend>` object. To maintain compatibility with HiCPipe-formatted 'mat' files, :mod:`HiFive` expects that fend and fragment numbering begin at index 1, not 0.

To create the HiC dataset from raw coordinate data, you can run the following commands::

  import hifive
  data = hifive.HiCData(out_filename, mode='w')
  data.load_data_from_raw(fend_filename, [raw1.txt, raw2.txt], maxinsert=500)
  data.save()

In this case, 'out_filename' specifies the location to save the :class:`HiCData <hifive.hic_data.>` object to. The 'fend_filename' value is the location of the appropriate :class:`Fend <hifive.fend.Fend>` object. Multiple files containing paired-end mapped coordinates may be passed to the function as a list or, if only a single file is needed, it may be passed as a string.

In order to load data from a set of BAM files, a similar procedure is used::

  import hifive
  data = hifive.HiCData(out_filename, mode='w')
  data.load_data_from_bam(fragment_filename,
    [[bam_file1_1, bam_prefix1_2], [bam_file2_1, bam_file2_2]],
    maxinsert=500)
  data.save()

In this case, the only difference is that pairs of file names corresponding to the two mapped read ends are passed as lists. Like the function for counts data, if only a single pair of files is needed, it may be passed as a list (not nested).

If your Fend file is 'binned', then you can also load data directly from a set of tab-delimited matrix files. These files can contain labels indicating bin positions (see :ref:`matrix_files`). If no labels are present, each column and row is expected to match the paritioning in the Fend file and start with the first bin of the chromosome(s). This is done using the command::

  data.load_data_from_matrices(fragment_filename,
    ['chr1.matrix', 'chr2.matrix', 'chr1_by_chr2.matrix'])

If only matrix files are to be loaded, a Fend file created using chromosome lengths is the best option as it does not contain fend data and gaurentees that bins start with the zero coordinate which is how most publicly available matrix files are organized.

.. note:: The :class:`HiCData <hifive.hic_data.>` object can now be used by multiple analyses of this sample and does not need to be created separately for each one.

.. _creating_a_HiC_project_object:

Creating a :class:`HiC <hifive.hic.HiC>` project object
=========================================================

The HiC project object, :class:`HiC <hifive.hic.HiC>`, contains links to a HiCData and :class:`Fend <hifive.fend.Fend>` object, information about which fends to include in the analysis, model parameters, and learned model values. This is the standard way of working with HiC data in :mod:`HiFive` and this object will be used for learning the model, extracting portions of data, plotting, and downstream analysis.

To create a :class:`HiC <hifive.hic.HiC>` object, you can use the following commands::

  import hifive
  hic = hifive.HiC(out_filename, 'w')
  hic.load_data(data_filename)
  hic.save()

In this case, 'out_filename' specifies the location to save the :class:`HiC <hifive.hic.HiC>` object to. The 'data_filename' value is the location of the appropriate data object.

.. warning:: Because data and fragment data are stored in their own objects, each object keeps track of the location of its dependents through relative file names. This means that links between them will break if the relative pathway is changed.

.. _filter_HiC_fends:

Filter HiC fends
=====================

Prior to modeling the data, you need to filter out fends that have few valid reads mapped to them. :mod:`HiFive` uses an iterative filtering approach such that only when all fends satisfy a user-defined minimum number of valid interactions does the filtering process cease.

To filter fends, you can use the following commands::

  import hifive
  hic = hifive.HiC(hic_filename)   
  hic.filter_fends(mininteractions=25, mindistance=5000000)
  hic.save()

In this case, 'hic_filename' is a previously saved :class:`HiC <hifive.hic.HiC>` analysis object. No value was passed to mode, since it defaults to 'r' for read. This loads the data from a previously created HiCData object. In order for changes to be kept to a FiveC object, it must be written to file using the save command. The 'mindistance' argument specifies that only reads associated with interactions spanning that distance or more are counted for purposes of filtering fends.

.. _find_HiC_distance_function:

Find HiC distance function
============================

:mod:`HiFive` approximates the distance-signal relationship using a series of linear transitions between bin means of mean log interaction counts. Spanning from a user-defined minimum interaction distance up to the genome maximum interaction distance, the range is divided into equal-sized log distance bins. Values falling between bin midpoints are interpolated based on a linear transition between bins. To estimate this function, you can use the following command::

  hic.find_distance_parameters(numbins=90,
                               minsize=200, 
                               maxsize=0)

In this function call, the range of interaction sizes is being broken into 90 bins, 1 bin covering interactions <= 200 bp, and the other 89 spanning up to the maximum interaction distance with breaks evenly spaced in log space. The maximum of this range is set by 'maxsize', which can either be zero, as in this call, setting the maximum size equal to the longest interaction distance, or a positive integer value which would exclude any interaction distances greater than 'maxsize'.

Because this function involves scanning large amounts of data, it has been made to utilize MPI. To do so, you can use a scripts such as the following::

  import hifive
  from mpi4py import MPI

  rank = MPI.COMM_WORLD.Get_rank()
  hic = hifive.HiC(hic_filename)
  hic.find_distance_means(numbins=90,
                          minsize=200, 
                          maxsize=0)
  if rank == 0:
    hic.save()

.. _learn_HiC_normalization_parameters:

Learn HiC normalization parameters
===================================

Using the probability algorithm
+++++++++++++++++++++++++++++++

In order to learn the correction model for HiC data using the probability algorithm, :mod:`HiFive` uses a bactracking line gradient descent. In addition, :mod:`HiFive` limits which interactions it uses to learn the model parameters to those that fall within a user-specified maximum interaction distance.

To learn HiC corrections using the modeling approach, you can use the following command::

  hic.find_probability_fend_corrections(mindistance=5000000,
                                        learningstep=0.5,
                                        max_iterations=1000,
                                        minchange=0.0005)

In the above call, 'mindistance' indicates that interactions spanning less than 5 Mb are excluded from calculations. Setting this to zero would include all unfiltered cis interactions. The 'learningstep' specifies how quickly to scale down the step value if the current try doesn't meet the arjimo learning criterion. The 'max_iterations' specifies a limit for how long to run the learning process for. Finally, 'minchange' is the stopping threshold such that if all absolute gradient values are below this the learning terminates early.

Because of the large numbers of calculations involved in this function, it has been made to utilize MPI. To do so, you can use a scripts such as the following::

  import hifive
  from mpi4py import MPI

  rank = MPI.COMM_WORLD.Get_rank()
  hic = hifive.HiC(hic_filename)
  hic.find_fend_corrections(mindistance=5000000,
                            learningstep=0.5,
                            max_iterations=1000,
                            minchange=0.0005)
  if rank == 0:
    hic.save()

Using the express algorithm
+++++++++++++++++++++++++++++++

:mod:`HiFive` also offers an express algorithm based on a matrix-balancing approach for learning correction values. The primary differences to the correction model from the user's perspective are a single learning phase and a lack of learning rate. The approximation learning approach can still recalculate the distance function parameters periodically.

To learn HiC corrections using the approximation approach, you can use the following command::

  hic.find_express_fend_corrections(iterations=1000,
                                    mindistance=0,
                                    usereads='cis',
                                    remove_distance=True)

In the above call, 'mindistance' is used to exclude interaction distances shorter that the passed value. If this results in the exclusion of any reads, fends are refiltered using either the value passed under the keyword 'mininteractions' or, if that is not specified, the value passed the last time fends were filtered. The 'usereads' argument allows the user to base the correction value approximation on 'cis' interactions, 'trans' interactions, or 'all'. Selecting 'trans' interactions will also result in a refiltering of fends to ensure that all of them are involved in sufficient interactions as described previously. The 'remove_distance' argument specifies whether to remove the distance-dependent portion of the signal prior to approximating correction values. For best results, this should set to true (its default value).

Although this function is much more computationally efficient, the calculation of the distance-dependence signal estimates can take time and so has been made to utilize the MPI environment when available as follows::

  import hifive
  from mpi4py import MPI

  rank = MPI.COMM_WORLD.Get_rank()
  hic = hifive.HiC(hic_filename)
  hic.find_express_fend_corrections(iterations=1000,
                                    mindistance=0,
                                    usereads='cis',
                                    remove_distance=True)
  if rank == 0:
    hic.save()

Using the binning algorithm
+++++++++++++++++++++++++++++++

:mod:`HiFive` also offers a fend characteristic-based approach adapted from the learning model used by `HiCPipe <http://www.wisdom.weizmann.ac.il/~eitany/hicpipe/>`_. This algorithm takes a list of features to be partitioned and a number of bins to partition them into and learns correction values associated with each partition based on a binomial distribution of binary data (observed / not observed).

To learn HiC corrections using the binning approach, you can use the following command::

  hic.find_binning_fragment_corrections(max_iterations=1000,
                                        mindistance=5000000,
                                        maxdistance=0,
                                        num_bins=[10, 10],
                                        model=['len', 'gc'],
                                        parameters=['even', 'fixed-const'],
                                        usereads='cis',
                                        learning_threshold=1.0)

Unlike the other two learning algorithms, this approach caps the learning iterations using 'max_iterations' and provides a means of early termination. This is done with the 'learning_threhold' parameter, which specifies that if the change in log-likelihood drops below 1.0, then cease iterating. The 'model', 'num_bins', and 'parameters' values should consist of equal-length lists and describe the correction values that are to be learned. Here, we told HiFive to use the length and gc content (specified in our BED file) for each fend. Each feature was partitioned into a number of bins specified in 'num_bins'. The partitioning of length was done to create bins containing equal numbers of fends while the gc content was divided such that each bin spanned an equal portion of the characteristic's range. Finally, the '-const' suffix told HiFive not to optimize the values for gc content. The 'usereads' value 'cis' specified that only intra-chromosomal interactions should be used to learn these correction values.

Although learning the correction values with this algorithm is much more computationally efficient, calculating the number of observations and possible observations per bin must span a large number of fend combinations so this function has been made to utilize the MPI environment when available as follows::

  import hifive
  from mpi4py import MPI

  rank = MPI.COMM_WORLD.Get_rank()
  hic = hifive.HiC(hic_filename)
  hic.find_binning_fragment_corrections(max_iterations=1000,
                                        mindistance=5000000,
                                        maxdistance=0,
                                        num_bins=[10, 10],
                                        model=['len', 'gc'],
                                        parameters=['even', 'fixed-const'],
                                        usereads='cis',
                                        learning_threshold=1.0)
  if rank == 0:
    hic.save()

Chaining learning algorithms
++++++++++++++++++++++++++++++

Because they work in very different ways, :mod:`HiFive` allows the binning algorithm to be chained with either the probability or express algorithm. The binning learning algorithm must be performed first and the corrections are applied prior to learning corrections with the second algorithm. This can be done by using the 'precorrect' option as follows::

  fivec.find_binning_fragment_corrections(max_iterations=1000,
                                          mindistance=5000000,
                                          maxdistance=0,
                                          num_bins=[20],
                                          model=['len','gc'],
                                          parameters=['even','even'],
                                          usereads='cis',
                                          learning_threshold=1.0)
  fivec.find_express_fragment_corrections(iterations=1000,
                                          mindistance=5000000,
                                          maxdistance=0,
                                          remove_distance=True,
                                          precorrect=True)

.. _generating_a_hic_heatmap:

Generating a heatmap
====================

In order to immediately make use of data, :mod:`HiFive` allows you to pull data from a region and create a heatmap. The data can be returned unbinned, binned using a fixed-width bin size, or binned using boundaries passed by the user. There are  several options for the format the data can be passed back in. Please refer to the :meth:`cis_heatmap <hifive.hic.HiC.cis_heatmap>` function for more details. There are also several options for transformations to the data. These are used to remove the distance-dependence signal, fend bias, both, or to return only the predicted signal. In this example, we'll get a portion of chromosome 1 binned into 10 Kb bins as follows::

  heatmap = hic.cis_heatmap(chrom='1',
                            start=1000000
                            stop=3000000
                            binsize=10000,
                            arraytype='upper',
                            datatype='enrichment')

In the above call, All valid possible interactions were queried from chromosome 1 between 1000000 and 3000000. For valid interactions that had no observation, an expected value was still added to the bin. 'enrichment' specifies to find the observed counts and expected counts, which includes the distance-dependence and fend bias values. The observed counts are in the first index of the last dimension of the returned array, the expected counts are in the second index of the last dimension. 'Upper' specifies a row-major upper triangle array (all values above the matrix diagonal flattened).

Accessing heatmap data
======================

When a heatmap is generated, data are stored in an HDF5 dictionary, a binary file format that allows easy access through python. In order to access data from your heatmap, you can load it as follows::

  import h5py
  import numpy
  heatmap = h5py.File(heatmap_file, 'r')
  print heatmap.keys()
  counts = heatmap['0.counts'][...]
  expected = heatmap['0.expected'][...]
  enrichment = numpy.zeros((counts.shape[0], 2), dtype=numpy.float32)
  where = numpy.where(counts > 0)[0]
  enrichment[where, 0] = numpy.log(counts[where] / expected[where])
  enrichment[where, 1] = 1

Note that we used the 'r' option when opening the file with h5py. This ensures that we are in 'read' mode. You could also use 'a' for 'append' mode, which is the default. First we printed out the available dataset names in our heatmap file. These are all of the arrays that are accessible to us by calling them like any other key value in a dictionary. Next, in order to load data from the heatmap into memory rather than access it from the disk every time we refer to it, we use the '[...]' indexing call after pass the heatmap filestream the name of the data we want. In this case, we asked for the counts and expected values for region zero. In order to look at the enrichments, we took the log of the ratio of observed to expected values for each bin. However, there are likely bins that contain no observed counts which would give us a divide by zero error in the log function. So, we can use numpy's 'where' function to get a index list of places that match our criterion, in this case non-zero counts. Finally, we have made the enrichment array 2D so we can keep track of which bins are valid (nonzero counts). If we were looking at trans data, we would need one more dimension as the counts and expected arrays would be two-dimensional instead of one.

.. _plotting_a_hic_heatmap:

Plotting a heatmap
==================

In order to visualize the heatmap we just produced, :mod:`HiFive` has several plotting functions that take different shaped arrays. The function called needs to match the array produced. In this case, we produced an upper array which is compatible with the :func:`plot_upper_array<hifive.plotting.plot_upper_array>` function, so we'll use that as follows::

  img = hifive.plotting.plot_upper_array(heatmap, symmetric_scaling=True)
  img.save(out_fname)

In calling the function, we pass the heatmap and that would be sufficient. There are, however, additional options. For example, 'symmetric_scaling' specifies whether the color scale should be expanded to run from the minimum value to the maximum (False) or so that the maximum absolute value determine both upper and lower color bounds. The image returned is a :mod:`PIL` image of type 'png'.

Note that if we were plotting data from a 'binned' HiC dataset, we would have to pass the 'diagonal_included' option as True for either :func:`plot_upper_array<hifive.plotting.plot_upper_array>` or :func:`plot_compact_array<hifive.plotting.plot_compact_array>`.

.. _making_an_mrh_file:

Creating a multi-resolution heatmap
====================================

An Alternative to the standard heatmap is a HiFive-specific filetype called a :ref:`multiresolution heatmap` or MRH. In order to create this compact heatmap file, we can use the built-in :class:`HiC <hifive.hic.HiC>` function as follows::

  hic.load( project_filename )
  hic.write_multiresolution_heatmap(out_fname, datatype='fend', maxbinsize=1280000, minbinsize=5000, includetrans=True, minobservations=5)

This call would create an MRH file under the name specified in 'out_fname'. Data would cover all intra-chromosomal interactions and, because we passed 'True' to the 'includetrans' argument, all of the inter-chromosomal interactions as well. Only bins with at least 5 reads would be included in the heatmaps because of the value passed to 'minobservations'. All resolutions from binsizes of 1.28Mb to 5Kb would be heatmapped in steps of 2X (i.e. 5Kb, 10Kb, 20Kb, etc). This imposes a limitation such that minbinsize and maxbinsize must differ from each other by an integer power of two.

In order to make use of this MRH file, we could either visualize it in `Galaxy <https://usegalaxy.org/>`_ or use the stand-alone program 'fetch_mrh_data' included with HiFive. For example if we wanted to get data from chromosome 1 between 10Mb and 20Mb including all resolution data down to 10Kb, we could call the program as follows::

  > fetch_mrh_data -c 1 -s 10000000 -e 20000000 -R 10000 mrh_fname img_fname

This would pull data from the MRH file 'mrh_fname' and plot it as an image in the file 'img_fname'. For specifics of the 'fetch_mrh_data' options, see :ref:`mrh_program`.
