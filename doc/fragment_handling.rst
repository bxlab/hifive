.. _restriction fragments:

*****************************
Restriction Fragment Handling
*****************************

HiFive handles all data at the level of the restriction fragment via the :class:`Fragment <hifive.fragment.Fragment>` and :class:`Fend <hifive.fend.Fend>` classes. In the case of 5C data, this means the set of fragments that are specifically targeted by primers, whereas for HiC data, each restriction fragment is treated as two separate units, one for each half of the fragment. The reason for this is that each fragment end (fend) is associated with its own characteristics such as GC content and mappability. Further, because there is an inverse relationship between fragment length and overall interaction signal strength, successful ligation and sequencing is biased towards interaction sites (and therefore cross-linking sites) that occur closer to the end being ligated. Effectively this creates two different populations of interactions, those associated with each half of a restriction fragment.

All normalization and binning are accomplished from the fragment or fend level for 5C and HiC data, respectively. This means that once normalization is complete, it can be used for any level of resolution.

==========================
Loading 5C Fragments
==========================

Fragments associated with 5C data are loaded from a BED file that contains the chromosome, start and ending coordinates, name, and strand for each primer-targeted restriction fragment from the 5C experiment. Additionally, a FASTA file containing the primer sequences can also be provided to associate the GC content of each primer with its genome position data (used in HiFive's :ref:`binning algorithm`). To reduce the storage space and processing time, only fragment data associated with primer targets are used in HiFive's 5C Fragment objects, meaning that a different Fragment object is needed for each experimental design, but experiments sharing the same targeted fragments may share a Fragment object file.


=========================
Loading HiC Fends
=========================

Fends associated with HiC data can be loaded from either a BED file or a HiCPipe-compatible tabular fend file. If using a BED file, the file should contain either the chromosome, start, and stop position of each restriction fragment or the chromosome, start, and stop position of the restriction enzyme recognition sites for the target genome. HiFive will infer which type of BED file is given based on the coordinate intervals. If a HiCPipe-compatible fend file is given, it should contain a header line, fend, fragment, chromosome, coordinate and fragment length information for each fend in the target genome.

::

  fend    frag    chr    coord     frag_len
  1       1       1      3002506   3372
  2       1       1      3005877   3372
  3       2       1      3005878   389
  4       2       1      3006266   389

.. note::
  HiCPipe-style fend files are 1-indexed, meaning that the first fragment and first fend both are labeled with a 1. This convention is used in HiFive only for these files to maintain compatibility with HiCPipe files.

The header line should contain the exact labels as seen above since HiFive uses them to determine which columns contain what information. In addition to the above characteristics, the tabular fend file may also contain the columns 'frag_gc' and 'map_score'. These fend characteristic values are used in HiFive's :ref:`binning algorithm`, although are not needed for either probability or express normalization.