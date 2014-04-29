#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

"""
This is a module class for handling FiveC paired end data.

Input data
----------

This class loads data from a set of text files containing pairs of
chromosomes, coordinates, and strand, pairs of fends and counts, or directly
from a set of bam files. Data are filtered by fragment insert size at time of
loading. Using the "Fend" class, reads are assigned to fragment-end (fend)
pairs.

Concepts
--------

Data are stored in h5dicts to allow easy access, fast retrieval, and reduce
memory requirements.

-----------------------------------------------------------------------------

API documentation
-----------------



"""

import os
import sys
from glob import glob

import numpy
import h5py
try:
    import pysam
except:
    pass

from ..fragment import Fragment


class FiveCData(object):
    """Base class for handling count data for FiveC experiments.

    This class stores mapped paired-end reads, indexing them by fragment
    number, in an h5dict."""

    def __init__(self, filename, mode='r'):
        """
        __init__ method

        Initialize counts dataset and create an h5dict.

        Parameters
        ----------
        filename : string
            A filename specifying where to store the data dictionary.
        mode : string, optional
            Specifies how to open the h5dict, depending on whether data is to
            be written or read.
        """
        self.data = h5py.File(filename, mode)
        return None

    def load_data_from_counts(self, fragfilename, filelist):
        """
        load_data_from_counts method

        Read counts from text files and place in h5dict.

        Parameters
        ----------
        fragfilename : string
            This specifies the filename of the fragment object.
        filelist : list
            A list containing all of the file names of counts text files to be
            included in the dataset.
        """
        # determine if fragment file exists and if so, load it
        if not os.path.exists(fragfilename):
            print >> sys.stderr, \
              ("The fragment file %s was not found. No data was loaded.\n") % (fragfilename.split('/')[-1]),
            return None
        self.frags = fragfilename
        frags = h5py.File(fragfilename, 'r')
        strands = frags['fragments']['strand'][...]
        chr2int = {}
        for i, j in enumerate(frags['chromosomes'][:]):
            chr2int[j] = i
        # create fragment name dictionary
        names = {}
        for i in range(frags['fragments'].shape[0]):
            names[frags['fragments']['name'][i]] = i
        # load data from all files, skipping if name not in the fragment file.
        if isinstance(filelist, str):
            filelist = [filelist]
        total_reads = 0
        data = {}
        for fname in filelist:
            reads = 0
            if not os.path.exists(fname):
                print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                continue
            print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
            input = open(fname, 'r')
            for line in input:
                temp = line.strip('\n').split('\t')
                if temp[0] not in names or temp[1] not in names or temp[0] == temp[1]:
                    continue
                frag1 = names[temp[0]]
                frag2 = names[temp[1]]
                # if both in same orientation, skip
                if strands[frag1] == strands[frag2]:
                    continue
                pair = (min(frag1, frag2), max(frag1, frag2))
                if pair not in data:
                    data[pair] = 0
                data[pair] = int(temp[2])
                reads += int(temp[2])
            input.close()
            print >> sys.stderr, ("%i validly-mapped reads loaded.\n") % (reads),
            total_reads += reads
        if len(data) == 0:
            print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
        print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(data)),
        # write fragment pairs to h5dict
        self._write_fragment_pairs(frags, data)
        return None

    def load_data_from_bam(self, fragfilename, filelist):
        """
        load_data_from_bam method

        Read counts from pairs of bam files and place in h5dict.

        Parameters
        ----------
        fragfilename : string
            This specifies the filename of the fragment object.
        filelist : list
            A list containing all of the bam file prefices to be included in
            the dataset. All files containing each prefix will be loaded.
        """
        if 'pysam' not in sys.modules.keys():
            print >> sys.stderr, ("The pysam module must be installed to use this function.")
            return None
        # determine if fragment file exists and if so, load it
        if not os.path.exists(fragfilename):
            print >> sys.stderr, ("The fragment file %s was not found. No data was loaded.\n") %\
                                 (fragfilename.split('/')[-1]),
            return None
        self.frags = fragfilename
        frags = h5py.File(fragfilename, 'r')
        strands = frags['fragments']['strand'][...]
        chr2int = {}
        for i, j in enumerate(frags['chromosomes'][:]):
            chr2int[j] = i
        # create fragment name dictionary
        names = {}
        for i in range(frags['fragments'].shape[0]):
            names[frags['fragments']['name'][i]] = i
        # load data from all files, skipping if either fragment not not in the fragment file.
        if isinstance(filelist, str):
            filelist = [filelist]
        total_reads = 0
        data = {}
        for prefix in filelist:
            reads = 0
            # determine which files have both mapped ends present
            fnames = glob("%s*_1.bam*" % prefix)
            if len(fnames) == 0:
                print >> sys.stderr, ("No files found for prefix %s") % (prefix.split('/')[-1]),
                continue
            for fname in fnames:
                if not os.path.exists(fname.replace('_1.bam', '_2.bam')):
                    print >> sys.stderr, ("%s does not have a paired end.") % (fname.split('/')[-1]),
                    del fnames[fnames.index(fname)]
            unpaired = {}
            # load first half of paired ends
            for fname in fnames:
                print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
                input = pysam.Samfile(fname, 'rb')
                for read in input.fetch(until_eof=True):
                    # Only consider reads with an alignment
                    if read.is_unmapped:
                        continue
                    # if mapping name not in fragment names, skip
                    seq_name = input.getrname(read.tid)
                    if seq_name not in names:
                        continue
                    # skip multiply-aligned reads
                    for tag in read.tags:
                        if tag[0] == 'XS':
                            break
                    else:
                        unpaired[read.qname] = names[seq_name]
                input.close()
                print >> sys.stderr, ("Done\n"),
            # load second half of paired ends
            for fname1 in fnames:
                fname = fname1.replace('_1.bam', '_2.bam')
                print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
                input = pysam.Samfile(fname, 'rb')
                for read in input.fetch(until_eof=True):
                    # Only consinder reads whose paired end was valid
                    if read.qname not in unpaired:
                        continue
                    # Only consider reads with an alignment
                    if read.is_unmapped:
                        continue
                    # if mapping name not in fragment names, skip
                    seq_name = input.getrname(read.tid)
                    if seq_name not in names:
                        continue
                    # skip multiply-aligned reads
                    for tag in read.tags:
                        if tag[0] == 'XS':
                            break
                    else:
                        # if both ends map to the same orientation, skip
                        if strands[unpaired[read.qname]] != strands[names[seq_name]]:
                            pair = (min(unpaired[read.qname], names[seq_name]),
                                    max(unpaired[read.qname], names[seq_name]))
                            if pair not in data:
                                data[pair] = 0
                            data[pair] += 1
                            reads += 1
                        del unpaired[read.qname]
                input.close()
                print >> sys.stderr, ("Done\n"),
            print >> sys.stderr, ("Read %i validly_mapped read paired.\n") % (reads),
            total_reads += reads           
        if len(data) == 0:
            print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
        print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(data)),
        self._write_fragment_pairs(frags, data)
        return None

    def _write_fragment_pairs(self, frags, frag_pairs):
        """Separate frag pairs into cis (within region) and trans (between region) interactions and write to
        h5dict with index arrays."""
        print >> sys.stderr, ("Writing fragment pair data to file..."),
        cis = {}
        trans = {}
        for i in range(frags['regions'].shape[0]):
            # Find pairs from same region
            for j in range(frags['regions']['start_frag'][i], frags['regions']['stop_frag'][i] - 1):
                for k in range(j + 1, frags['regions']['stop_frag'][i]):
                    if (j, k) in frag_pairs:
                        cis[(j, k)] = frag_pairs[(j, k)]
            # Find pairs from different regions
            for j in range(frags['regions']['start_frag'][i], frags['regions']['stop_frag'][i]):
                for k in range(frags['regions']['stop_frag'][i], frags['fragments'].shape[0]):
                    if (j, k) in frag_pairs:
                        trans[(j, k)] = frag_pairs[(j, k)]
        # convert data into arrays
        cis_data = numpy.empty((len(cis), 3), dtype=numpy.int32)
        trans_data = numpy.empty((len(trans), 3), dtype=numpy.int32)
        keys = cis.keys()
        keys.sort()
        for i in range(len(keys)):
            cis_data[i, 0] = keys[i][0]
            cis_data[i, 1] = keys[i][1]
            cis_data[i, 2] = cis[keys[i]]
        keys = trans.keys()
        keys.sort()
        for i in range(len(keys)):
            trans_data[i, 0] = keys[i][0]
            trans_data[i, 1] = keys[i][1]
            trans_data[i, 2] = trans[keys[i]]
        # find first instance of each fragment for cis and trans data
        if cis_data.shape[0] > 0:
            cis_indices = numpy.r_[0, numpy.bincount(cis_data[:, 0],
                                   minlength=frags['fragments'].shape[0])].astype(numpy.int32)
            for i in range(1, cis_indices.shape[0]):
                cis_indices[i] += cis_indices[i - 1]
        if trans_data.shape[0] > 0:
            trans_indices = numpy.r_[0, numpy.bincount(trans_data[:, 0],
                                     minlength=frags['fragments'].shape[0])].astype(numpy.int32)
            for i in range(1, trans_indices.shape[0]):
                trans_indices[i] += trans_indices[i - 1]
        # write data to h5dict
        self.data.attrs['fragfilename'] = self.frags
        if cis_data.shape[0] > 0:
            self.data.create_dataset('cis_data', data=cis_data)
            self.data.create_dataset('cis_indices', data=cis_indices)
        if trans_data.shape[0] > 0:
            self.data.create_dataset('trans_data', data=trans_data)
            self.data.create_dataset('trans_indices', data=trans_indices)
        print >> sys.stderr, ("Done\n"),
        return None
