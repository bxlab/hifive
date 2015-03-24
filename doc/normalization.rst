.. _normalization:

*****************
Normalizing Data
*****************

HiFive offers several different approaches to normalizing 5C and HiC data. These include three different algorithmic approaches that can be used individually or chained together.

.. _probability algorithm:

Probability Algorithm
=============================

HiFive's probability-based algorithm models counts as arising from a probability distribution (Poisson and Log-Normal for HiC and 5C, respectively) with a mean equal to the product of the contributions from the distance-dependent signal (as determined by the distance-dependence estimation function) and the bias contributions from the fragment or fend at each end of the interaction. In order to learn the fragment or fend bias terms, HiFive pre-calculates the distance-dependence contribution for each cis interaction (intra-chromosomal or intra-regional for HiC and 5C, respectively). Each bias correction can be pre-estimated as the sum of counts (or log-counts in the case of 5C) for a fragment or fend divided by the sum of estimated distance-dependent signal for all interactions involving that fragment or fend.

Parameter optimization is done via a two-stage gradient descent. In the first phase (burn-in), the learning rate is kept constant. For HiC data, HiFive also allows a minimum parameter change cutoff to ensure that all parameter values have stabilized. This works by finding the parameter change from one iteration to the next. The first learning phase continues until this change is less than the specified threshold. The second phase (annealing) uses a linearly-decreasing learning rate ranging from the learning rate from the burn-in phase down to zero.

This algorithm can be limited to learning with a subset of data based on upper and lower distance cutoffs. This is useful for avoiding highly significant structures such as those occurring at short ranges. In the case of HiC data, this also allows a sub-sampling of data to fit within memory limitations for higher-resolution data. This is important as the HiC probability algorithm uses all interactions (both observed and unobserved) to learn correction parameters meaning that memory requirements expand very quickly.

.. _express algorithm:

Express Algorithm
==========================

HiFive's Express algorithm is a modified version of matrix balancing. Unlike the standard matrix balancing approach, takes into account the number of possible interactions for each fragment or fend such that differences between the number of considered interaction partners do not lead to bias in the correction values. This is particularly relevant as HiFive allows which values are used to be limited by upper and lower interaction distance cutoffs, resulting in potentially very different numbers of interaction partners for different fragments or fends depending on the restriction site density variation. Further, this allows use of the matrix-balancing approach with a sparse matrix, making it memory efficient and fast.

.. _binning algorithm:

Binning Algorithm
===========================

HiFive's binning algorithm is based directly on the multivariate binning probability modeling approach used by `HiCPipe <http://www.ncbi.nlm.nih.gov/pubmed/22001755>`_. A set of fragment characteristics (GC content, mappability - HiC only, and length) are broken into bins correction values for all bin combinations for each parameter are learned using the Broyden–Fletcher–Goldfarb–Shanno algorithm from the Scipy package. Distance can be used as an additional parameter, although this is broken into a one-dimensional set of bins. Distance bin parameters are used only during the algorithm learning stage to remove bias and are not kept for later normalization, which is instead taken care of by the distance-dependence estimation function.

This approach uses a binary readout of the data (observed / unobserved) rather than counts. Learning can be accomplished using cis, trans, or all interactions. In the case of cis interactions, lower and upper interaction distance cutoffs can be used to limit the set of interactions included in the learning process. If both cis and trans data are used and distance bins are part of the model, the last distance bin is reserved for trans interactions and cis interactions are assigned to the n - 1 remaining bins.

.. _chaining normalization:

Chaining Normalization Approaches
==================================
Because the corrections learned by the binning algorithm are fundamentally different from those found by the probability and express approaches, they can be combined to account for different aspects of systematic bias within the data. This is accomplished using a two step learning process. Initially the binning algorithm is performed, followed by either the probability or express algorithm using predicted values adjusted based on the binning corrections. This allows both approaches (combinatorial and multiplicative) to be used to best advantage.
