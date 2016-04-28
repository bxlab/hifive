.. _normalization:

*****************
Normalizing Data
*****************

HiFive offers several different approaches to normalizing 5C and HiC data. These include three different algorithmic approaches that can be used individually or chained together.

.. _probability algorithm:

Probability Algorithm
=============================

HiFive's probability-based algorithm models counts as arising from a probability distribution (Binomial or Poisson and Log-Normal for HiC and 5C, respectively) with a mean equal to the product of the contributions from the distance-dependent signal (as determined by the distance-dependence estimation function) and the bias contributions from the fragment or fend at each end of the interaction. In order to learn the fragment or fend bias terms, HiFive pre-calculates the distance-dependence contribution for each cis interaction (intra-chromosomal or intra-regional for HiC and 5C, respectively). Each bias correction can be pre-estimated as the sum of log-counts for a fragment divided by the sum of estimated distance-dependent signal for all interactions involving that fragment for 5C or the sum of interactions (binary, not counts) divided by the sum of the estimated binary distance-dependence signal for all of a fend's interactions.

Parameter optimization is done via a backtracking line gradient descent. This runs until either a maximum number of iterations is reached or all of the gradient absolute values fall below a cutoff threshold.

This algorithm can be limited to learning with a subset of data based on upper and lower distance cutoffs. This is useful for avoiding highly significant structures such as those occurring at short ranges. In the case of HiC data, this also allows a sub-sampling of data to fit within memory limitations for higher-resolution data. This is important as the HiC probability algorithm uses all interactions (both observed and unobserved) to learn correction parameters meaning that memory requirements expand very quickly.

.. _express algorithm:

Express Algorithm
==========================

HiFive's Express algorithm is actually two different but very similar algorithms. The first is a modified version of matrix balancing that weights correction updates by the number of valid interactions the associated fragment or fend is involved in (this can be different between fragments, especially when using distance range cutoffs). The other express algorithm is true matrix balancing and uses the `Knight and Ruiz algorithm <http://imajna.oxfordjournals.org/content/early/2012/10/26/imanum.drs019>`_ for learning corrections. Both versions allow the use of the matrix-balancing approach with a sparse matrix, making it memory efficient and fast. In addition, both allow correcting for distance dependence prior to learning for better results. In HiC data, values can be considered as true counts or as binary indicators of observed/unobserved. For 5C data, the data are far less sparse so the binary approach does not yield suitable results. Instead, counts may be learned either as is or as log-counts. In addition, for the Knight and Ruiz algorithm it is necessary to incorporate psuedo-counts along the diagonal in order achieve convergence.

.. _binning algorithm:

Binning Algorithm
===========================

HiFive's binning algorithm is based directly on the multivariate binning modeling approach used by `HiCPipe <http://www.ncbi.nlm.nih.gov/pubmed/22001755>`_. A set of fragment characteristics (GC content, length, mappability - HiC only, and distance - HiC only) are broken into bins correction values for all bin combinations for each parameter are learned using the Broyden–Fletcher–Goldfarb–Shanno algorithm from the Scipy package. In HiC, distance can be used as an additional parameter, although this is broken into a one-dimensional set of bins. Distance bin parameters are used only during the algorithm learning stage to remove bias and are not kept for later normalization, which is instead taken care of by the distance-dependence estimation function.

For HiC data, this approach uses a binary readout of the data (observed / unobserved) rather than counts and observation probability is assumed to be binomially distributed.  Learning can be accomplished using cis, trans, or all interactions. In the case of cis interactions, lower and upper interaction distance cutoffs can be used to limit the set of interactions included in the learning process. If both cis and trans data are used and distance bins are part of the model, the last distance bin is reserved for trans interactions and cis interactions are assigned to the n - 1 remaining bins. You also have the option to add pseudo-counts to the bins. This is done for each bin combination, meaning that for each model parameter, bin1 x bin2 will have twice as many counts added as bin1 x bin1 (to account for both possible combinations). This should be used with caution, as pseudo-counts tend to negatively impact the quality of normalization for cis interactions and HiC datasets typically have enough observations to alleviate the need for pseudo-counts.

For 5C data, only nonzero counts are used and are assumed to be log-normally distributed. For each read, the prior is taken as the predicted distance-dependence signal. Learning can be accomplished using cis, trans, or all interactions. In the case of cis interactions, lower and upper interaction distance cutoffs can be used to limit the set of interactions included in the learning process.

.. note:: Currently this algorithm is not supported by binned HiC data.

.. _chaining normalization:

Chaining Normalization Approaches
==================================
Because the corrections learned by the binning algorithm are fundamentally different from those found by the probability and express approaches, they can be combined to account for different aspects of systematic bias within the data. This is accomplished using a two step learning process. Initially the binning algorithm is performed, followed by either the probability or express algorithm using predicted values adjusted based on the binning corrections. This allows both approaches (combinatorial and multiplicative) to be used to best advantage.
