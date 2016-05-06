.. _command_line:

*************
Using HiFive
*************

HiFive contains a variety of subcommands for handling different aspects of 5C and HiC analysis. To learn about a specific subcommand, you can use the following command::

  > hifive <subcommand> -h

The available subcommands are:

:fragments:               Create a fragment file from a BED file containing targeted RE fragment data.
:5c-data:                 Create a data file from mapped BAM or fragment-pair counts files.
:5c-project:              Create a project file, filter fragments, and estimate distance-dependence.
:5c-normalize:            Find correction parameter values using one of the available algorithms.
:5c-complete:             Perform all of the steps of the subcommands fragments, 5c-data, 5c-project, and 5c-normalization in one command.
:5c-heatmap:              Using an already created 5C project, generate an HDF5-formatted heatmap file and optional image.
:5c-interval:             Using an already created 5C project, generate a tabular genomic-interval file for a specified region and optional image.
:5c-combine-replicates:   Combine multiple 5C data files into a single file without needing to reload the data.
:fends:                   Create a fend file from either a BED or HiCPipe-style fend file containing RE fragment data or create an arbitrarily-binned interval file from chromosome length file.
:hic-data:                Create a data file from mapped BAM, MAT, or paired coordinate text (RAW) files or from binned matrix files.
:hic-project:             Create a project file, filter fends, and estimate distance-dependence.
:hic-normalize:           Find correction parameter values using one of the available algorithms.
:hic-complete:            Perform all of the steps of the subcommands fends, hic-data, hic-project, and hic-normalization in one command.
:hic-heatmap:             Using an already created HiC project, generate an HDF5-formatted heatmap file and optional image.
:hic-interval:            Using an already created HiC project, generate a tabular genomic-interval or matrix file for a specified region and optional image.
:hic-combine-replicates:  Combine multiple HiC data files into a single file without needing to reload the data.
:hic-mrheatmap:           Create a multi-resolution heatmap file from a HiFive HiC project file.

.. _5c_subcommands:

5C Subcommands
===================

Each subcommand uses its own set of options although there is some overlap of options between subcommands.

.. _fragments:

fragments
+++++++++

::

  > hifive fragments [-h] [-r RE] [-g GENOME] [-q] bed output

Arguments:

:bed:  A BED file with targeted restriction enzyme fragment coordinates and the associated name of the primer used to target each fragment. Additional fields can be included with RE fragment characterstics (such as GC content) in columns after the strand column. Including additional features requires a header with the name of each feature.
:output: A filename to write the HiFive Fragment file to.

Options:

-h/--help, -r/--re, -g/--genome, -q/--quiet

.. _5c_data:

5c-data
++++++++

::

  > hifive 5c-data [-h] (-B BAM BAM | -C COUNT) [-q] fragment output

Arguments:

:fragment: The HiFive Fragment file to associate with this dataset.
:output: A filename to write the HiFive 5C dataset to.

Options:

-h/--help, -B/--bam, -C/--count, -q/--quiet

.. _5c_project:

5c-project
++++++++++

::

  > hifive 5c-project [-h] [-f MININT] [-m MINDIST] [-x MAXDIST] [-q]
        data output

Arguments:

:data: The HiFive 5C dataset file to associate with this project.
:output: A filename to write the HiFive 5C project to.

Options:

-h/--help, -f/--min-interactions, -m/--min-distance, -x/--max-distance, -q/--quiet

.. _5c_normalize:

5c-normalize
++++++++++++

::

  > hifive 5c-normalize <SUBCOMMAND> [-h] [-m MINDIST] [-x MAXDIST]
        [-r REGIONS] [-o OUTPUT] [-q] [normalization options] project

HiFive's 5C normalization subcommand allows the user to select the normalization approach to use. Each approach has its own set of options.

Arguments:

:project: The HiFive 5C project to find fragment corrections for.

Options:

-h/--help, -m/--min-distance, -x/--max-distance, -r/--regions, -o/--output, -q/--quiet

Subcommands:

probability, express, binning, probability-binning, express-binning, binning-probability, binning-express

.. _5c_complete:

5c-complete
+++++++++++

::

  > hifive 5c-complete <SUBCOMMAND> [-h] [-r RE] [-g GENOME]
        (-B BAM BAM | -C COUNT) [-f MININT] [-m MINDIST] [-x MAXDIST]
        [-r REGIONS] (-o OUTPUT OUTPUT OUTPUT | -P PREFIX) [-q]
        [normalization options] bed

HiFive's complete 5C analysis subcommand allows the user to select the normalization approach to use. Each approach has its own set of options.

Arguments:

:bed:  A BED file with targeted restriction enzyme fragment coordinates and the associated name of the primer used to target each fragment.

Options:

-h/--help, -r/--re, -g/--genome, -B/--bam, -C/--count, -f/--min-interactions, -m/--min-distance, -x/--max-distance, -r/--regions, -o/--output, -P/--prefix -q/--quiet

Subcommands:

probability, express, binning, probability-binning, express-binning, binning-probability, binning-express

.. _5c_heatmap:

5c-heatmap
++++++++++

::

  > hifive 5c-heatmap [-h] [-b BINSIZE] [-t] [-r REGIONS]
        [-d {raw,fragment,distance,enrichment,expected}]
        [-a {compact,full}] [-y] [-x EXPBINSIZE] [-f MINOBS]
        [-g SEARCH] [-v] [-i IMAGE] [-p] [-l] [-n]
        [-k KEYWORDS] [-q] project output

Arguments:

:project: The HiFive 5C project to create a heatmap for.
:output: The filename to write the HiFive 5C heatmap to. 

Options:

-h/--help, -b/--binsize, -t/--trans, -r/--regions, -d/--datatype, -a/arraytype, -y/--dynamically-bin, -x/--expansion-binsize, -f/--minobservations, -g/--search-distance, -v/--remove-failed, -i/--image, -p/--pdf, -l/--legend, -n/--names, -k/--keyword, -q/--quiet

.. _5c_interval:

5c-interval
+++++++++++

::

  > hifive 5c-interval [-h] -c REGION [-s START] [-e STOP] [--region2 REGION2]
        [--start2 START2] [--stop2 STOP2] [-b BINSIZE]
        [-d {raw,fragment,distance,enrichment,expected}] [-y] [-x EXPBINSIZE]
        [-f MINOBS] [-g SEARCH] [-v] [-i IMAGE] [-p] [-r] [-t] [-l]
        [-k KEYWORDS] [-q] project output

Arguments:

:project: The HiFive 5C project to create a heatmap for
:output: The filename to write the HiFive 5C genomic interval file to. 

Options:

-h/--help, -c/--region, -s/--start, -e/--stop, -b/--binsize, -d/--datatype, -y/--dynamically-bin, -x/--expansion-binsize, -f/--minobservations, -g/--search-distance, -v/--remove-failed, -i/--image, -p/--pdf, -r/--rotate, -t/--ticks, -l/--legend, -k/--keyword, -q/--quiet

.. _5c_combine_replicates:

5c-combine-replicates
+++++++++++++++++++++

::

  > hifive 5c-combine-replicates [-h] [-q] output replicate
        [replicate ...]

Arguments:

:output: The filename to write the new HiFive 5C dataset file to. 
:replicate: A HiFive 5C dataset file.

Options:

-h/--help, -q/--quiet

.. _5c_options:

5C Options
===========

Universal Options:

-h, --help   Display the help message and command/subcommand options and arguments and exit.
-q, --quiet  Suppress all messages generated during HiFive processing.

5C Fend Options:

-r, --re str      The name of the restriction enzyme.
-g, --genome str  The name of the genome.

5C Data Options:

-B, --bam FILES   A pair of BAM filenames separated by spaces corresponding to the two independently-mapped ends of a set of reads. Multiple file pairs may be passed by calling this argument more than once. This option is mutually exclusive with -C/--count.
-C, --count FILE  A tabular text file containing pairs of fragment primer names and their associated read count (see `Loading 5C Data <loading_data.html>`_ for more information). This option is mutually exclusive with -B/--bam.

5C Project Options:

-f, --min-interactions int  The minimum number of interactions with valid fragments to keep a fragment in the analysis. [20]
-m, --min-distance int      The minimum distance between fragment midpoints to include in calculating numbers of interactions for fragment filtering and (if called by 5c-normalization or 5c-complete) the minimum interaction distance included in learning correction parameter values. [0]
-x, --max-distance int      The maximum distance between fragment midpoints to include in calculating numbers of interactions for fragment filtering and (if called by 5c-normalization or 5c-complete) the maximum interaction distance included in learning correction parameter values. A value of zero indicates no maximum distance cutoff. [0]

5C Normalization Options:

-r, --regions str   A comma-separated list of region numbers to include fragments from when calculating correction parameter values. [all regions]
-o, --output FILE   An optional filename to save the updated HiFive project to, leaving the original unchanged. [None]

5C Complete Options:

-o, --output FILES  A set of three filenames separated by spaces to save the newly-created HiFive fragment, dataset, and project files to. Mutually exclusive with -P/--prefix.
-P, --prefix str    A prefix for the output filenames. The file extensions .frags, .fcd, and .fcp will be used for the fragment, dataset, and project files, respectively. This option is mutually exclusive with -o/--output.

5C Normalization Algorithms
+++++++++++++++++++++++++++

5C Probability Options:

-b, --max-iterations int     The maximum number of iterations to run the learning process for. [1000]
-g, --min-change dec         The minimum allowable absolute gradient size to coninute learning process. [0.0005]
-p, --precalculate           Prior to beginning learning, set initial guesses for each correction value to be learned to the fragment's mean difference between its log-counts and predicted distance-dependence signal.
-l, --learning-step dec      The scaling factor for decreasing learning rate by if step doesn't meet Armijo criterion. [0.5]

5C Express Options:

-e, --express-iterations int  The number of iterations to run the learning process for. [1000]
-d, --remove-distance         Calculate and subtract out the predicted distance-dependence signal from each log-count prior to learning correction parameters.
-w, --express-reads str       Which set reads to use for learning correction parameter values, cis, trans, or all. [cis]
-k, --logged                  Use log-counts instead of counts for learning.
-z, --knight-ruiz             Use the Knight Ruiz matrix balancing algorithm instead of weighted matrix balancing. This option ignores 'iterations' and 'logged'.

5C Binning Options:

-i, --binning-iterations int  The maximum number of iterations to run the learning process for. [1000]
-t, --learning-threshold dec  The maximum change in log-likelihood necessary to stop the learning process early. [1.0]
-y, --binning-reads str       Which set of reads to use for learning correction parameter values, cis, trans, or all. [cis]
-v, --model str               A comma-separated list of fragment features to calculate corrections for. Acceptable features are len (length) and any features loaded in the BED file used to create the HiFive fragment file. [len]
-n, --model-bins str          A comma-separated list of numbers of bins to partition fragment features into for modeling. [10]
-u, --parameter-types str     A comma-separated list of model parameter types. Acceptable values are even, fixed, even-const, and fixed-const. Even means that fragment features are partitioned such that each bin has approximately even numbers of fragments. Fixed means that the range of the feature is divided into fixed-width bins. The -const suffix indicates that the correction values are held at their seed-values and not updated. [even]

5C Interaction Binning Options
++++++++++++++++++++++++++++++

5C Heatmap Options:

-b, --binsize int            The width of bins (in basepairs) to partition data into. A value of zero indicates that each bin is to correspond with a single fragment. [10000]
-t, --trans                  Calculate and include trans interactions in heatmaps.
-r, --regions str            A comma-separated list if region numbers to include in the heatmaps. [all regions]
-d, --datatype str           Type of data to produce for the heatmaps. Valid options are raw, fragment (only fragment corrections applied), distance (only distance-dependence signal removed), enrichment (both fragment correction and distance-dependence signal removed), and expected (only predicted signal). [fragment]
-a, --arraytype str          If data is unbinned, this option specifies whether the heatmaps should be full or compact. Full means that there is a row and column for every fragment, while compact means that rows are forward fragments only and columns are reverse fragments only. [full]
-y, --dynamically-bin        Dynamically bin heatmap.
-x, --expansion-binsize int  The size of bins, in base pairs, to group data into for expanding under-populated bins. [10000]
-f, --minobservations int    The minimum number of observed reads in a bin for it to be considered valid. [20]
-g, --search-distance int    The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [0]
-v, --remove-failed          If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.

5C Interval Options:

-c, --region int             The index of the region to pull data from.
-b, --binsize int            The width of bins (in basepairs) to partition data into. A value of zero indicates that each bin is to correspond with a single fragment.
-s, --start int              The first coordinate of the region to pull data from. None indicates the beginning of the region. [None]
-e, --stop int               The last coordinate + 1 of the region to pull data from. None indicates the end of the region. [None]
-y, --dynamically-bin        Dynamically bin heatmap.
-x, --expansion-binsize int  The size of bins, in base pairs, to group data into for expanding under-populated bins. [10000]
-f, --minobservations int    The minimum number of observed reads in a bin for it to be considered valid. [20]
-g, --search-distance int    The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [0]
-v, --remove-failed          If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.

5C Plotting Options:

-i, --image FILE    Generate an image from the region or regions for which heatmap data is being calculated. [None]
-p, --pdf           Format the image as a pdf. [None]
-r, --rotate        Rotate the image 45 degrees so the chromosome axis is horizontal and only plot the triangle above this axis. This option can only be used with a full arraytype.
-t, --ticks         Add coordinate ticks and labels to heatmap. This option can only be used if a pdf is requested.
-l, --legend        Add a color scale bar corresponding to interaction strength. This option can only be used if a pdf is requested.
-n, --names         Add region names to the plot. This option can only be used if a pdf is requested.
-k, --keyword str   Pass additional plotting options accepted by the :mod:`plotting <hifive.plotting>` module. Arguments should be of the format KEYWORD=VALUE. This option can be passed multiple times. [None]

.. _hic_subcommands:

HiC Subcommands
===================

Each subcommand uses its own set of options although there is some overlap of options between subcommands.

.. _fends:

fends
+++++++++

::

  > hifive fends [-h] (-F FEND | -B BED | -L LENGTH) [--binned] [-r RE] [-g GENOME] [-q] output

Arguments:

:output: A filename to write the HiFive Fend file to.

Options:

-h/--help, -F/--fend, -B/--bed, -L/--length, --binned, -r/--re, -g/--genome, -q/--quiet

.. _hic_data:

hic-data
++++++++

::

  > hifive hic-data [-h] (-S BAM BAM | -R RAW | -M MAT | -X MATRIX) [-i INSERT] [--skip-duplicate-filtering] [-q]
        fend output

Arguments:

:fragment: The HiFive Fend file to associate with this dataset.
:output: A filename to write the HiFive HiC dataset to.

Options:

-h/--help, -S/--bam, -R/--raw, -M/--mat, -X/--matrix, -i/--insert, --skip-duplicate-filtering, -q/--quiet

.. _hic_project:

hic-project
++++++++++++

This command is MPI-compatible.

::

  > [mpirun -np NP] hifive hic-project [-h] [-f MININT] [-m MINDIST]
                        [-x MAXDIST] [-j MINBIN] [-n NUMBINS] [-q] data
                        output

Arguments:

:data: The HiFive HiC dataset file to associate with this project.
:output: A filename to write the HiFive HiC project to.

Options:

-h/--help, -f/--min-interactions, -m/--min-distance, -x/--max-distance, -j/--min-binsize, -n/--num-bins, -q/--quiet

.. _hic_normalize:

hic-normalize
++++++++++++++

This command is MPI-compatible.

::

  > [mpirun -np NP] hifive hic-normalize <SUBCOMMAND> [-h] [-m MINDIST]
                        [-x MAXDIST] [-c CHROMS] [-o OUTPUT] [-q]
                        [normalization options] project

HiFive's HiC normalization subcommand allows the user to select the normalization approach to use. Each approach has its own set of options.

Arguments:

:project: The HiFive HiC project to find fragment corrections for.

Options:

-h/--help, -m/--min-distance, -x/--max-distance, -c/--chromosomes, -o/--output, -q/--quiet

Subcommands:

probability, express, binning, binning-probability, binning-express

.. _hic_complete:

hic-complete
+++++++++++++

This command is MPI-compatible.

::

   > [mpirun -np NP] hifive hic-complete <SUBCOMMAND> [-h]
                        (-F FEND | -B BED | -L LENGTH) [--binned]
                        [-r RE] [-g GENOME]
                        (-S BAM BAM | -R RAW | -M MAT | -X matrix)
                        [-i INSERT] [--skip-duplicate-filtering]
                        [-f MININT] [-m MINDIST] [-x MAXDIST]
                        [-j MINBIN] [-n NUMBINS] [-c CHROMS]
                        (-o OUTPUT OUTPUT OUTPUT | -P PREFIX) [-q]
                        [normalization options]

HiFive's complete HiC analysis subcommand allows the user to select the normalization approach to use. Each approach has its own set of options.

Options:

-h/--help, -F/--fend, -B/--bed, -L,--length, --binned, -r/--re, -g/--genome, -S/--bam, -R/--RAW, -M/--mat, -X/--matrix, -i/--insert, --skip-duplicate-filtering, -f/--min-interactions, -m/--min-distance, -x/--max-distance, -j/--min-binsize, -n/--num-bins, -c/--chromosomes, -o/--output, -P/--prefix -q/--quiet

Subcommands:

probability, express, binning, binning-probability, binning-express

.. _hic_heatmap:

hic-heatmap
++++++++++++

This command is MPI-compatible.

::

  > [mpirun -np NP] hifive hic-heatmap [-h] [-b BINSIZE] [-t]
                        [-c CHROMS]
                        [-d {raw,fend,distance,enrichment,expected}]
                        [-F {hdf5,txt,npz}] [-y] [-x EXPBINSIZE]
                        [-f MINOBS] [-a SEARCH] [-v]  [-i IMAGE]
                        [-p] [-l] [-n] [-k KEYWORDS]
                        [-q] project output

Arguments:

:project: The HiFive HiC project to create a heatmap for.
:output: The filename to write the HiFive HiC heatmap to. 

Options:

-h/--help, -b/--binsize, -t/--trans, -c/--chromosomes, -d/--datatype, -F/--format, -y/--dynamically-bin, -x/--expansion-binsize, -f/--minobservations, -a/--search-distance, -v/--remove-failed, -i/--image, -p/--pdf, -l/--legend, -n/--names, -k/--keyword, -q/--quiet

.. _hic_interval:

hic-interval
+++++++++++++

::

  > hifive hic-interval [-h] -c CHROM [-s START] [-e STOP] [-b BINSIZE]
        [-m MAXDIST] [-d {raw,fend,distance,enrichment,expected}] [-M]
        [-y] [-x EXPBINSIZE] [-f MINOBS] [-a SEARCH] [-v] [-i IMAGE] [-p]
        [-r] [-t] [-l] [-k KEYWORDS] [-q] project output

Arguments:

:project: The HiFive HiC project to create a heatmap for.
:output: The filename to write the HiFive HiC genomic interval file to. 

Options:

-h/--help, -c/--chromosome, -s/--start, -e/--stop, -b/--binsize, -m/--max-distance, -d/--datatype, -M/--matrix, -y/--dynamically-bin, -x/--expansion-binsize, -f/--minobservations, -a/--search-distance, -v/--remove-failed, -i/--image, -p/--pdf, -r/--rotate, -t/--ticks, -l/--legend, -k/--keyword, -q/--quiet

.. _hic_combine_replicates:

hic-combine-replicates
+++++++++++++++++++++++

::

  > hifive hic-combine-replicates [-h] [-q] replicate1 replicate2 output

Arguments:

:replicate1: The first HiFive HiC dataset file to be combined.
:replicate2: The second HiFive HiC dataset file to be combined.
:output: The filename to write the new HiFive HiC dataset file to. 

Options:

-h/--help, -q/--quiet

.. _hic_mrheatmap:

hic-mrheatmap
++++++++++++++

::

 > hifive hic-mrheatmap [-h] [-t] [-c CHROMS] [-f MINOBS] [-B MAXBIN]
       [-b MINBIN] [-R MAXTRANSBIN] [-r MINTRANSBIN] [-m MIDBIN]
       [-d {raw,fend,distance,enrichment}] [-q] project output

Arguments:

:project: The HiFive HiC project to create a multi-resolution heatmap for.
:output: The filename to write the multi-resolution heatmap to.

Options:

-h/--help, -q/--qiuet, -t/--trans, -c/--chromosomes, -f/--minobservations, -B/--maximum-binsize, -b/--minimum-binsize, -R/--maximum-trans-binsize, -r/--minimum-trans-binsize, -m/--mid-binsize, -d/--datatype, 

.. _hic_options:

HiC Options
===========

Universal Options:

-h, --help   Display the help message and command/subcommand options and arguments and exit.
-q, --quiet  Suppress all messages generated during HiFive processing.

HiC Fend Options:

-F, --fend FILE    A tabular file in a format compatible with HiCPipe containing fragment and fend indices, fragment length, start or end position, and any additional fragment features desired (see `Loading HiC Fends <fragment handling.html>`_ for more information).
-B, --bed FILE     A BED file containing either restriction enzyme fragment coordinates or retriction enzyme cutsite coordinates. Fragment features may be included in columns after the strand column. Features should be formatted with one feature per column and two values per feature separated by a comma. If the coordinates are of RE fragment boundaries, the feature values should correspond to the upstream end of the fragment followed by the downstream end of the fragment. If the coordinates are of RE cutsites, the values should correspond to the sequence just upstream of the cutsite followed by the sequence just downstream of the cutsite. If additional features are included, the bed file must have a header line identifying the features.
-L, --length FILE  A tab-separated text file containing chromosome names and lengths. Must be used in conjunction with a positive value of 'binned'.
--binned int       Indicates what size bins to break genome into. If None is passed, fend-level resolution is kept.
-r, --re str       The name of the restriction enzyme.
-g, --genome str   The name of the genome.

HiC Data Options:

-S, --bam FILES             A pair of BAM filenames separated by spaces corresponding to the two independently-mapped ends of a set of reads. Multiple file pairs may be passed by calling this argument more than once. This option is mutually exclusive with -R/--raw and -M/--mat.
-R, --raw FILE              A tabular file containing pairs of mapped read positions (see `Loading HiC Data <loading_data.html>`_ for more information).
-M, --mat FILE              A tabular file containing pairs of fend indices and their corresponding numbers of reads (see `Loading HiC Data <loading_data.html>`_ for more information).
-X, --matrix FILE           A tab-separated binned matrix containing summed fend interactions.
-i, --insert int            The maximum allowable insert size, as measured by the sum of both read end mapping positions to the nearest RE cutsite in the direction of alignment.
--skip-duplicate-filtering  Skip filtering of PCR duplicates (only applicable to raw and bam files).

HiC Project Options:

-f, --min-interactions int  The minimum number of interactions with valid fends to keep a fend in the analysis. [20]
-m, --min-distance int      The minimum distance between fend midpoints to include in calculating numbers of interactions for fend filtering and (if called by hic-normalization or hic-complete) the minimum interaction distance included in learning correction parameter values. [0]
-x, --max-distance int      The maximum distance between fend midpoints to include in calculating numbers of interactions for fend filtering and (if called by hic-normalization or hic-complete) the maximum interaction distance included in learning correction parameter values. A value of zero indicates no maximum distance cutoff. [0]
-j, --min-binsize int       The cutoff size limit for the smallest distance bin used for estimating the distance dependence (see `HiC Distance Dependence Estimation <distance_dependence.html>`_ for more information). [1000]
-n, --num-bins int          The number of bins to partition the interaction size ranges into for estimating the distance dependence function (see `HiC Distance Dependence Estimation <distance_dependence.html>`_ for more information). A value of zero indicates that finding the distance dependence function should be skipped.

HiC Normalization Options:

-c, --chromosomes str   A comma-separated list of chromosome names to include fends from when calculating correction parameter values. [all chromosomes]
-o, --output FILE   An optional filename to save the updated HiFive project to, leaving the original unchanged. [None]

HiC Complete Options:

-o, --output FILES  A set of three filenames separated by spaces to save the newly-created HiFive fend, dataset, and project files to. Mutually exclusive with -P/--prefix.
-P, --prefix str    A prefix for the output filenames. The file extensions .fends, .hcd, and .hcp will be used for the fragment, dataset, and project files, respectively. This option is mutually exclusive with -o/--output.

HiC Normalization Algorithms
+++++++++++++++++++++++++++++

HiC Probability Options:

-b, --max-iterations int        The maximum number of iterations to run the learning process for. [1000]
-g, --min-change dec            The minimum allowable absolute gradient size to coninute learning process. [0.0005]
-p, --precalculate              Prior to beginning learning, set initial guesses for each correction value to be learned to the fragment's mean difference between its log-counts and predicted distance-dependence signal.
-l, --learning-step dec         The scaling factor for decreasing learning rate by if step doesn't meet Armijo criterion. [0.5]
-a, --probability-model         Which probability model to use for normalization (binomial or poisson).

HiC Express Options:

-e, --express-iterations int  The number of iterations to run the learning process for. [1000]
-d, --remove-distance         Calculate and divide out the predicted distance-dependence signal from each count prior to learning correction parameters.
-w, --express-reads str       Which set reads to use for learning correction parameter values, cis, trans, or all. [cis]
-g, --min-change              The minimum mean change in fend correction parameter values needed to keep running past 'iterations' number of iterations. If using the Knight-Ruiz algorithm this is the residual cutoff.
-f, --min-interations int     The minimum number of interactions for fend filtering if refiltering is required due to distance cutoff parameters or selected reads to be used. [20]
-k, --binary bool             Use binary indicator instead of counts.
-z, --knight-ruiz bool        Use the Knight Ruiz matrix balancing algorithm instead of weighted matrix balancing. This option ignores 'iterations'.

HiC Binning Options:

-r, --binning-iterations int  The maximum number of iterations to run the learning process for. [1000]
-t, --learning-threshold dec  The maximum change in log-likelihood necessary to stop the learning process early. [1.0]
-y, --binning-reads str       Which set of reads to use for learning correction parameter values, cis, trans, or all. [cis]
-v, --model str               A comma-separated list of fend features to calculate corrections for. Acceptable features are len (length), distance, and any features loaded in the BED or FEND file used to create the HiFive fend file. [len,distance]
-s, --model-bins str          A comma-separated list of numbers of bins to partition fend features into for modeling. [20,20]
-u, --parameter-types str     A comma-separated list of model parameter types. Acceptable values are even, fixed, even-const, and fixed-const. Even means that fend features are partitioned such that each bin has approximately even numbers of fends. Fixed means that the range of the feature is divided into fixed-width bins. The -const suffix indicates that the correction values are held at their seed-values and not updated. [even,fixed-const]
--pseudocounts int            The number of pseudo-counts to add to each bin prior to seeding and learning normalization values. [None]

HiC Interaction Binning Options
++++++++++++++++++++++++++++++++

HiC Heatmap Options:

-b, --binsize int            The width of bins (in basepairs) to partition data into. A value of zero indicates that each bin is to correspond with a single fend. [10000]
-t, --trans                  Calculate and include trans interactions in heatmaps.
-c, --chromosomes str        A comma-separated list if chromosome names to include in the heatmaps. [all chromosomes]
-d, --datatype str           Type of data to produce for the heatmaps. Valid options are raw, fend (only fend corrections applied), distance (only distance-dependence signal removed), enrichment (both fend correction and distance-dependence signal removed), and expected (only predicted signal). [fend]
-F, --format str             The format of the output heatmap. Valid options are hdf5, txt, and npz. [hdf5]
-M, --matrix                 Store output as a tab-separated matrix of values.
-y, --dynamically-bin        Dynamically bin heatmap.
-x, --expansion-binsize int  The size of bins, in base pairs, to group data into for expanding under-populated bins. [10000]
-f, --minobservations int    The minimum number of observed reads in a bin for it to be considered valid. [20]
-a, --search-distance int    The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [0]
-v, --remove-failed          If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.

HiC Interval Options:

-c, --chromosome str    The chromosome to pull data from.
-b, --binsize int       The width of bins (in basepairs) to partition data into. A value of zero indicates that each bin is to correspond with a single fend.
-s, --start int         The first coordinate of the chromosome to pull data from. None indicates the beginning of the chromosome. [None]
-e, --stop int          The last coordinate + 1 of the chromosome to pull data from. None indicates the end of the chromosome. [None]
-m, --max-distance int  The largest interaction distance to include in the interval file. A value of zero indicates no upper limit. [0]
-d, --datatype str      Type of data to produce for the heatmaps. Valid options are raw, fend (only fend corrections applied), distance (only distance-dependence signal removed), enrichment (both fend correction and distance-dependence signal removed), and expected (only predicted signal). [fend]
-y, --dynamically-bin        Dynamically bin heatmap.
-x, --expansion-binsize int  The size of bins, in base pairs, to group data into for expanding under-populated bins. [10000]
-f, --minobservations int    The minimum number of observed reads in a bin for it to be considered valid. [20]
-a, --search-distance int    The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [0]
-v, --remove-failed          If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.

HiC Plotting Options:

-i, --image FILE    Generate an image from the region or regions for which heatmap data is being calculated. [None]
-p, --pdf           Format the image as a pdf. [None]
-r, --rotate        Rotate the image 45 degrees so the chromosome axis is horizontal and only plot the triangle above this axis.
-t, --ticks         Add coordinate ticks and labels to heatmap. This option can only be used if a pdf is requested.
-l, --legend        Add a color scale bar corresponding to interaction strength. This option can only be used if a pdf is requested.
-n, --names         Add chromosome names to the plot. This option can only be used if a pdf is requested.
-k, --keyword str   Pass additional plotting options accepted by the :mod:`plotting <hifive.plotting>` module. Arguments should be of the format KEYWORD=VALUE. This option can be passed multiple times. [None]

HiC Multi-Resolution Heatmap Options:

-t, --trans                      Calculate and include trans interactions in heatmaps.
-c, --chromosomes str            A comma-separated list if chromosome names to include in the heatmaps. [all chromosomes]
-f, --minobservations int        The minimum number of observed reads in a bin for it to be considered valid. [20]
-B, --maximum-binsize int        The largest sized bin to use (minimum resolution) in, base pairs. [1280000]
-b, --minimum-binsize int        The smallest sized bin to use (maximum resolution) in, base pairs. [10000]
-R, --maximum-trans-binsize int  The largest sized bin to use (minimum resolution) for inter-chromosomal interactions, in base pairs. If not specified, this defaults to the value of the -B option. [use -B value]
-r, --minimum-trans-binsize int  The smallest sized bin to use (maximum resolution) for inter-chromosomal interactions, in base pairs. If not specified, this defaults to the value of the -b option. [use -b value]
-m, --mid-binsize                The smalled sized bin to use for binning the entire chromosome, in base pairs. This is used to balance memory usage vs. speed and does not affect the output. [40000]
-d, --datatype str               Type of data to produce for the heatmaps. Valid options are raw, fend (only fend corrections applied), distance (only distance-dependence signal removed), enrichment (both fend correction and distance-dependence signal removed), and expected (only predicted signal). [fend]
