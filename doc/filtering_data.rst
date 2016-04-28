
.. _filtering_data:

************************
Filtering Data
************************

.. _fivec_filtering:

5C Filtering
============

HiFive filters reads from 5C experiments in two stages.

1. When a data object is created and data are loaded, only valid reads are read in. In 5C, this means only reads that represent an interaction between a forward-probed restriction fragment and a reverse-probed restriction fragment.

2. HiFive calls a function (:func:`filter_fragments <hifive.fivec.FiveC.filter_fragments>`) after a project is created and loaded with data that iteratively removes reads to produce a set of fragments that meet a minimum threshold of interactions with other valid fragments within the same region. This filtering is done using the number of interacting fragments rather than reads to account for the fact that fragment biases can produce wildly different numbers of interactions between fragment pairs. This filtering can also be limited to interactions falling within some range limitation of interaction distance. This can be used to ensure that fragments have sufficient reads for learning correction parameters when distance range constraints are used for normalization.

.. _hic_filtering:

HiC Filtering
=============

HiFive filters reads from HiC experiments in two stages.

1. Read-based filtering
------------------------

When a data object is created and data are loaded, only valid reads are read in. In HiC, this means that reads with identical mapping coordinates at both ends are assumed to be PCR duplicates and only one read is counted.

**Restriction fragment-based filtering**

 Reads originating from restriction enzyme-digested DNA are also filtered by the insert size as determined by the sum of the distances from the mapped positions of both ends to the nearest downstream target restriction site for each. Because of the possibility of incomplete restriction enzyme digestion and fragment circularization, reads with ends mapping to the same fragment and reads with ends mapping to adjacent fragments on opposite strands are also excluded from the data object.

 .. image:: _static/hic_filter.png

.. _bin filtering:

**Arbitrary or uniform bin-based filtering**

 Reads produced through non-specific breaking of the DNA (nuclease, random shearing, DNase, etc.) must be filtered in a different way than those produced by the standard HiC approach. In order to address the fact that some reads may originate from fragments that never underwent ligation or fragments that were circularized, HiFive allows a 'maxinsert' flag to be passed that subjects reads of lengths shorter than this flag to be filtered based on strand orienation. This filter removes reads whose ends map in opposite orienations, suggesting that possibility of a non-ligated or self-ligated fragment. Fragments with ends mapped to the same orienation are kept as a breakage and ligation even is gaurenteed if the read is not the result of a mapping error.

2. Fend or bin-based filtering
-------------------------------

 HiFive calls a function (:func:`filter_fragments <hifive.hic.HiC.filter_fends>`) after a project is created and loaded with data that iteratively removes reads to produce a set of fends that meet a minimum threshold of interactions with other valid fends on the same chromosome. This filtering is done using the number of interacting fends rather than reads to account for the fact that fend biases can produce wildly different numbers of interactions between fend pairs. This filtering can also be limited to interactions falling within some range limitation of interaction distance. This can be used to ensure that fends have sufficient reads for learning correction parameters when distance range constraints are used for normalization.
