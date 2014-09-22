.. _additional_scripts:


===================
Additional Scripts
===================

To make using :mod:`HiFive` easier, several  scripts are included to perform all of the standard functions.

--------------

create_fragmentset.py
---------------------

This script produces a :class:`Fragment` h5dict from a BED file containing restriction fragment boundaries and primer names targeting them. This script can be called as follows::

  > python create_fragmentset.py BED_FILE OUT_FILE [GENOME RE]

**Arguments:**

* **BED_FILE** (*str.*) - File containing 5C targeted fragments and primer names in BED format.
* **OUT_FILE** (*str.*) - File name to write :class:`Fragment` h5dict to.
* **GENOME** (*str.*) - Name of the genome.
* **RE** (*str.*) - Name of the restriction enzyme.

**Returns:**  None

--------------

create_fendset.py
-----------------

This script produces a :class:`Fend` h5dict from a HiCPipe-compatible fend file or BED file containing restriction fragment boundaries or sites. This script can be called as follows::

  > python create_fendset.py FEND_FILE OUT_FILE [GENOME RE]

**Arguments:**

* **FEND_FILE** (*str.*) - File containing restriction fragment data in HiCPipe-compatible or BED format.
* **OUT_FILE** (*str.*) - File name to write :class:`Fend` h5dict to.
* **GENOME** (*str.*) - Name of the genome.
* **RE** (*str.*) - Name of the restriction enzyme.

**Returns:**  None

--------------

create_fivec_dataset.py
-----------------------

This script produces a :class:`FiveCData` h5dict from a set of counts files or BAM files containing mapped paired-end reads. This script can be called as follows::

  > Usage python create_fivec_set.py DATA_FILE OUT_FILE MIN_INTERACTIONS

**Arguments:**

* **DATA_FILE** (*str.*) - File name of :class:`FiveCData` h5dict to link with analysis.
* **OUT_FILE** (*str.*) - File name to write :class:`FiveC` h5dict to.
* **MIN_INTERACTIONS** (*int.*) - Minimum number of interactions needed for valid fragment.

**Returns:**  None

--------------

create_hic_dataset.py
---------------------

This script produces a :class:`HiCData` h5dict from a set of BAM files containing mapped paired-end reads, HiCPipe-compatible MAT-formatted data files or raw paired-coordinate text files. This script can be called as follows::

  > python create_hic_dataset.py FEND_FILE DATA_FILE_1[,...,DATA_FILE_N] OUT_FILE MAX_INSERT

**Arguments:**

* **FEND_FILE** (*str.*) - File name of :class:`Fend` h5dict to link with data.
* **DATA_FILE_1[,...,DATA_FILE_N]** (*str.*) - A comma-separated list of either BAM file prefices, raw coordinate read pairs or HiCPipe-compatible MAT files.
* **OUT_FILE** (*str.*) - File name to write :class:`HiCData` h5dict to.
* **MAX_INSERT** (*int.*) - Integer specifying the maximum distance sum from each mapped end to restriction site.

**Returns:**  None

--------------

combine_replicates.py
---------------------

This script combines reads from multiple replicate :class:`HiCData` h5dicts and creates a new h5dict. This script can be called as follows::

  > python combine_replicates.py REP_FILE_1,REP_FILE_2[,...,REP_FILE_N] OUT_FILE

**Arguments:**

* **REP_FILE_1,REP_FILE2** (*str.*) - A comma-separated list of :class:`HiCData` h5dict files.
* **OUT_FILE** (*str.*) - File name to write :class:`HiCData` h5dict to.

**Returns:**  None

--------------

data2mat.py
-----------

This script exports read data from a :class:`HiCData` h5dict into a HiCPipe-compatible MAT-formatted text file. This script can be called as follows::

  > python data2mat.py DATA_FILE OUT_FILE

**Arguments:**

* **DATA_FILE** (*str.*) - File name of :class:`HiCData` h5dict.
* **OUT_FILE** (*str.*) - File name to write HiCPipe-compatible MAT-formatted data to.

**Returns:**  None

--------------

create_fivec_set.py
-------------------

This script creates a :class:`FiveC` h5dict analysis object, filters fragments, and calculates the distance dependence function. This script can be called as follows::

  > python create_fivec_set.py DATA_FILE OUT_FILE MIN_INTERACTIONS

**Arguments:**

* **DATA_FILE** (*str.*) - File name of :class:`FiveCData` h5dict to link with analysis.
* **OUT_FILE** (*str.*) - File name to write :class:`FiveC` h5dict to.
* **MIN_INTERACTIONS** (*int.*) - Minimum number of interactions needed for valid fragment.

**Returns:**  None

--------------

create_hic_set.py
-----------------

This script creates a :class:`HiC` h5dict analysis object, filters fends, and calculates the distance dependence function. This script can be called as follows::

  > python create_hic_set.py DATA_FILE HIC_FILE MIN_INTERACTIONS MAX_DIST MIN_SIZE NUM_BINS SMOOTHED

**Arguments:**

* **DATA_FILE** (*str.*) - File name of :class:`HiCData` h5dict to link with analysis.
* **OUT_FILE** (*str.*) - File name to write :class:`HiC` h5dict to.
* **MIN_INTERACTIONS** (*int.*) - Minimum number of interactions needed for valid fend.
* **MAX_DIST** (*int.*) - The largest interaction distance to be included for filtering fends.
* **MIN_SIZE** (*int.*) - The smallest interaction distance bin size for distance function.
* **NUM_BINS** (*int.*) - The number of bins to partion interaction distance range into for distance function.
* **SMOOTHED** (*int.*) - Number of adjacent bins to include for smoothing of distance function line.

**Returns:**  None

.. note:: This function is MPI compatible.

--------------

learn_fivec_normalization.py
----------------------------

This script learns fragment correction values for a :class:`FiveC` analysis object. This script can be called as follows::

  > python learn_fivec_normalization.py FIVEC_FILE RATE BURNIN ANNEALING MAX_DIST RECALC DISPLAY

**Arguments:**

* **FIVEC_FILE** (*str.*) - File name of :class:`FiveC` h5dict to analyze.
* **RATE** (*float*) - Percent of gradient to use for updating parameter values.
* **BURNIN** (*int.*) - Number of iterations to run burn-in phase for.
* **ANNEALING** (*int.*) - Number of iterations to run annealing phase for.
* **MAX_DIST** (*int.*) - Maximum interaction distance to include in learning.
* **RECALC** (*int.*) - Number of iterations to wait between recalculating distance function parameters.
* **DISPLAY** (*int.*) - Number of iterations to wait before explicitly calculating cost and updating display.

**Returns:**  None

--------------

learn_fivec_normalization_express.py
------------------------------------

This script learns fragment correction values for a :class:`FiveC` analysis object using the approximation approach. This script can be called as follows::

  > python learn_fivec_normalization_express.py FIVEC_FILE ITERATIONS REMOVE_DIST RECALC

**Arguments:**

* **FIVEC_FILE** (*str.*) - File name of :class:`FiveC` h5dict to analyze.
* **ITERATIONS** (*int.*) - Number of iterations to run learning for.
* **REMOVE_DIST** (*bool.*) - Specifies whether to remove distance-dependent portion of the signal prior to learning.
* **RECALC** (*int.*) - Number of iterations to wait between recalculating distance function parameters.

**Returns:**  None

--------------

learn_hic_normalization.py
----------------------------

This script learns fend correction values for a :class:`HiC` analysis object. This script can be called as follows::

  > python learn_hic_normalization.py HIC_FILE BURNIN ANNEALING MAX_DIST RECALC RATE DISPLAY

**Arguments:**

* **HIC_FILE** (*str.*) - File name of :class:`HiC` h5dict to analyze.
* **BURNIN** (*int.*) - Number of iterations to run burn-in phase for.
* **ANNEALING** (*int.*) - Number of iterations to run annealing phase for.
* **MAX_DIST** (*int.*) - Maximum interaction distance to include in learning.
* **RECALC** (*int.*) - Number of iterations to wait between recalculating distance function parameters.
* **RATE** (*float*) - Percent of gradient to use for updating parameter values.
* **DISPLAY** (*int.*) - Number of iterations to wait before explicitly calculating cost and updating display.

**Returns:**  None

.. note:: This function is MPI compatible.

--------------

learn_hic_normalization_express.py
------------------------------------

This script learns fend correction values for a :class:`HiC` analysis objectusing the approximation approach. This script can be called as follows::

  > python learn_hic_normalization_express.py HIC_FILE ITERATIONS MIN_INT MIN_DIST USE_READS REMOVE_DISTANCE RECALC

**Arguments:**

* **HIC_FILE** (*str.*) - File name of :class:`HiC` h5dict to analyze.
* **ITERATIONS** (*int.*) - Number of iterations to run learning for.
* **MIN_INT** (*int.*) - Minimum number of interactions for fend filtering, if refiltering is required.
* **MIN_DIST** (*int.*) - Minimum interaction distance to include for learning.
* **USE_READS** (*str.*) - Which set of reads, 'cis', 'trans', or 'both', to use for learning.
* **REMOVE_DISTANCE** (*bool.*) - Specifies whether to remove distance-dependent portion of the signal prior to learning.
* **RECALC** (*int.*) - Number of iterations to wait between recalculating distance function parameters.

**Returns:**  None

.. note:: This function is MPI compatible.

--------------

create_hic_heatmap_h5dict.py
----------------------------

This script creates an h5dict file containing binned heatmaps from a :class:`HiC` h5dict. This script can be called as follows::

  > python create_hic_heatmap_h5dict.py HIC_FILE OUT_FILE BINSIZE INCLUDE_TRANS REMOVE_DISTANCE CHROMS

* **HIC_FILE** (*str.*) - File name of a :class:`HiC` h5dict to pull data from.
* **OUT_FILE** (*str.*) - File name of heatmap h5dict to write data to.
* **BINSIZE** (*int.*) - Size of bins, in base pairs, to group data into.
* **INCLUDE_TRANS** (*bool.*) - Specifies whether to find inter-chromosome interactions.
* **REMOVE_DISTANCE** (*bool.*) - Specifies whether to remove distance-dependent portion of signal.
* **CHROMS** (*str.*) - Comma-separated list of chromosomes to find heatmaps for.

**Returns:**  None

--------------

find_hic_BI.py
--------------

This scripts takes a :class:`HiC` file and calculates a set of BI scores. This script can be called as follows::

  > python find_hic_BI.py HIC_FILE OUT_FILE WIDTH WINDOW HEIGHT MINCOUNT SMOOTHING [CHROM_1,...,CHROM_N]

**Arguments:**

* **HIC_FILE** (*str.*) - H5dict created by the :class:`HiC` class.
* **OUT_FILE** (*str.*) - File name for the new :class:`BI` h5dict created by this script.
* **WIDTH** (*int.*) - Integer specifying the width about each boundary point.
* **HEIGHT** (*int.*) - Integer specifying the height of bins extending across each window.
* **WINDOW** (*int.*) - Integer specifying the window around each boundary point.
* **SMOOTHING** (*int.*) - Integer specifying the width of smoothing weights.
* **CHROM_1,...CHROM_N** (*str.*) - A comma-separated list of chromosome names to include in the analysis. Optional.

**Returns:**  None

.. note:: This function is MPI compatible.

--------------

combine_BIs.py
--------------

This script takes two :class:`BI` files with different coordinates, such as would be created by two different restriction enzymes, annd comnines the data to create a composite set of scores. The script can be called as follows::

  > python combine_BIs.py BI_FILE_1 BI_FILE_2 OUT_FILE SMOOTHING CHROM_1[,...,CHROM_N]

**Arguments:**

* **BI_FILE_1** (*str.*) - The path of an h5dict file created by a :class:`BI` object.
* **BI_FILE_2** (*str.*) - The path of an h5dict file created by a :class:`BI` object.
* **OUT_FILE** (*str.*) - The path to write the new h5dict file to.
* **SMOOTHING** (*int.*) - The width, in base pairs, for smoothing BI scores.
* **CHROM_1[..,CHROM_N]** (*str.*) - A comma-separated list of chromosome names to include in the analysis.

**Returns:**  None

.. note:: This function is MPI compatible.

--------------

model_single_chr_BI.py
----------------------

This script bins data from a :class:`HiC` h5dict using peaks calls from a :class:`BI` object to partition signal, dynamically bins the data, and creates a 3D model using a PCA dimensionality reduction. The script can be called as follows::

  > python model_single_chr_BI.py HIC_FILE BI_FILE OUT_PREFIX CHROM CUTOFF MIN_OBS

**Arguments:**

* **HIC_FILE** (*str.*) - File name of :class:`HiC` h5dict to pull data from.
* **BI_FILE** (*str.*) - File name of :class:`BI` h5dict to find boundaries for partitioning from.
* **OUT_PREFIX** (*str.*) - File prefix for all output files of script.
* **CHROM** (*str.*) - Name of chromosome to model.
* **CUTOFF** (*float*) - Criteria for calling BI peaks.
* **MIN_OBS** (*int.*) - Minimum number of observations for valid dynamic bins.

**Returns:**  None

--------------

model_single_chr_binned.py
--------------------------

This script bins data from a :class:`HiC` h5dict, dynamically bins the data, and creates a 3D model using a PCA dimensionality reduction. The script can be called as follows::

  > python model_single_chr_BI.py HIC_FILE BI_FILE OUT_PREFIX CHROM CUTOFF MIN_OBS

**Arguments:**

* **HIC_FILE** (*str.*) - File name of :class:`HiC` h5dict to pull data from.
* **OUT_PREFIX** (*str.*) - File prefix for all output files of script.
* **BIN_SIZE** (*str.*) - Size of bins, in base pairs, to group data into.
* **MIN_OBS** (*int.*) - Minimum number of observations for valid dynamic bins.
* **CHROM** (*str.*) - Name of chromosome to model.

**Returns:**  None

--------------

model_whole_genome_BI.py
------------------------

This script bins data from a :class:`HiC` h5dict using peaks calls from a :class:`BI` object to partition signal, dynamically bins the data, and creates a 3D model using a PCA dimensionality reduction. The script can be called as follows::

  > python model_whole_genome_BI.py HIC_FILE BI_FILE OUT_PREFIX CUTOFF MIN_OBS CIS_SCALING CHROMS

**Arguments:**

* **HIC_FILE** (*str.*) - File name of :class:`HiC` h5dict to pull data from.
* **BI_FILE** (*str.*) - File name of :class:`BI` h5dict to find boundaries for partitioning from.
* **OUT_PREFIX** (*str.*) - File prefix for all output files of script.
* **CUTOFF** (*float*) - Criteria for calling BI peaks.
* **MIN_OBS** (*int.*) - Minimum number of observations for valid dynamic bins.
* **CIS_SCALING** (*float*) - Scaling factor to adjust cis interactions by prior to modeling.
* **CHROMS** (*str.*) - Comma-separated list of names of chromosomes to model.

**Returns:**  None

.. note:: This function is MPI compatible.

--------------

model_whole_genome_binned.py
----------------------------

This script bins data from a :class:`HiC` h5dict, dynamically bins the data, and creates a 3D model using a PCA dimensionality reduction. The script can be called as follows::

  > python model_whole_genome_binned.py HIC_FILE BI_FILE OUT_PREFIX BIN_SIZE MIN_OBS CIS_SCALING CHROMS

**Arguments:**

* **HIC_FILE** (*str.*) - File name of :class:`HiC` h5dict to pull data from.
* **BI_FILE** (*str.*) - File name of :class:`BI` h5dict to find boundaries for partitioning from.
* **OUT_PREFIX** (*str.*) - File prefix for all output files of script.
* **BIN_SIZE** (*str.*) - Size of bins, in base pairs, to group data into.
* **MIN_OBS** (*int.*) - Minimum number of observations for valid dynamic bins.
* **CIS_SCALING** (*float*) - Scaling factor to adjust cis interactions by prior to modeling.
* **CHROMS** (*str.*) - Comma-separated list of names of chromosomes to model.

**Returns:**  None

.. note:: This function is MPI compatible.
