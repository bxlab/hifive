.. _binning:

*******************
Binning Strategies
*******************

HiFive uses several different binning formats, depending on the nature of the data being handled. This is both for ease of manipulation and for size issues.

++++++++++++++++++
The 'upper' format
++++++++++++++++++

The 'upper' format is used when the data can be represented in a symmetric matrix (bin i,j is equal to j,i). In these cases, HiFive flattens the matrix into a 1D array, using only the upper triangle of the matrix, not including the diagonal. This results in an N * (N - 1) / 2 sized matrix from an N x N square matrix. This is the format used for all cis genomic regions that are binned and have no distance bound on the length of included interactions.

.. note:: For 'binned' HiC data, because within-bin interactions are included, the array is of size N * (N + 1) / 2.

++++++++++++++++++
The 'full' format
++++++++++++++++++

The 'full' format is a 2D array that can represent either symmetric or non-symmetric data. This is used to represent trans interaction data and can be requested as the return format of cis interactions. In the case of unbinned 5C data, the 'full' format contains one row and column for every fragment. However, because only forward-reverse combinations are valid the resulting matrix is fairly sparse, akin to a checkerboard pattern.

++++++++++++++++++++++++
The HiC 'compact' format
++++++++++++++++++++++++

Because HiC data is often worked with over short ranges, HiFive has a 'compact' format that has an upper bound on the length of interactions included in the array. The array size is equal to the number of bins (or fends in unbinned data) by the maximum number of bins apart that interactions are included for minus one. This means that if we have ten 10 Kb bins and include interactions up to 30 Kb apart, we would need a 10 x 3 array. Position i,j would give us interactions from bin i interacting with bin i + j + 1. This format allows a more efficient way of holding data up to the point that the maximum interaction distance is half of the span of the region or greater, in which case it becomes more efficient to use the 'upper' format. This is only used internally by HiFive or through the API.

.. note:: For 'binned' HiC data, because within-bin interactions are included, the array is one larger along the second axis.

++++++++++++++++++++++++++
The 5C 'compact' format
++++++++++++++++++++++++++

Because of the orientation associated with 5C fragment probes, HiFive uses a 'compact' format for unbinned 5C data that has forward fragments along the first axis and reverse fragments along the second axis of the data array. For trans interactions, a pair of regions is represented by two 'compact' arrays, one for each region's set of forward primers by the other region's reverse primers.

++++++++++++++++++++++++++
Dynamic Binning
++++++++++++++++++++++++++

HiFive offers a unique feature to address the variety in read density across 5C and HiC experiments in the form of dynamic binning. Rather than selecting an arbitrary bin size and contending with bins that don't contain enough or any observations, HiFive has an option for expanding each bin to meet a user-defined threshold of observations. This expansion can draw from unbinned data or binned data and allows each bin in a heatmap to contain sufficient information to be reasonably reliable. For each bin, the bounds expand in all directions equally, incorporating new observations based on expansion bin or fragment midpoints. After each incorporation, the observation threshold is checked and if satisfied, the expansion is ceased. A maximum search distance is also used to prevent excessive running time. If binned data is used for the expansion, the bin will encounter new observations in all directions simultaneously. If unbinned data are used, each expansion will include the nearest fragment or fend and its associated observed and expected values.
