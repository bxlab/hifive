#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

"""
This is a module class for handling Hi-C paired end data.

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

import numpy
import h5py
import pysam
from glob import glob

from ..fend import Fend


class HiCData(object):
    """Base class for handling count data for HiC experiments.

    This class stores mapped paired-end reads, indexing them by fend
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

    def load_data_from_raw(self, fendfilename, filelist, maxinsert):
        """
        load_data_from_raw method

        Read counts from raw text files and place in h5dict.

        Parameters
        ----------
        fendfilename : string
            This specifies the filename of the fend object.
        filelist : list
            A list containing all of the file names of raw text files to be
            included in the dataset.
        maxinsert : int
            A cutoff for filtering paired reads whose total distance to their
            respective restriction sites exceeds this value.
        """
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            print >> sys.stderr, \
              ("The fend file %s was not found. No data was loaded.\n") % (fendfilename.split('/')[-1]),
            return None
        self.fends = fendfilename
        self.maxinsert = maxinsert
        fends = h5py.File(fendfilename, 'r')
        chr2int = {}
        for i, j in enumerate(fends['chromosomes'][:]):
            chr2int[j] = i
        # load data from all files, skipping if chromosome not in the fend file.
        strand = {'+': 0, '-': 1}
        if isinstance(filelist, str):
            filelist = [filelist]
        fend_pairs = {}
        total_reads = 0
        for fname in filelist:
            data = {}
            if not os.path.exists(fname):
                print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                continue
            print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
            input = open(fname, 'r')
            for line in input:
                temp = line.strip('\n').split('\t')
                if temp[0] not in chr2int or temp[3] not in chr2int:
                    continue
                data[(chr2int[temp[0]], int(temp[1]), strand[temp[2]],
                      chr2int[temp[3]], int(temp[4]), strand[temp[5]])] = 0
            input.close()
            data = numpy.array(data.keys(), dtype=numpy.int32)
            print >> sys.stderr, ("%i validly-mapped reads loaded.\n") % (data.shape[0]),
            total_reads += data.shape[0]
            # map data to fends, filtering as needed
            if data.shape[0] > 0:
                self._find_fend_pairs(data, fends, fend_pairs)
        if len(fend_pairs) == 0:
            print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
        print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(fend_pairs)),
        # write fend pairs to h5dict
        self._write_fend_pairs(fends, fend_pairs)
        return None

    def load_data_from_bam(self, fendfilename, filelist, maxinsert):
        """
        load_data_from_bam method

        Read counts from pairs of bam files and place in h5dict.

        Parameters
        ----------
        fendfilename : string
            This specifies the filename of the fend object.
        filelist : list
            A list containing all of the bam file prefices to be included in
            the dataset. All files containing each prefix will be loaded.
        maxinsert : int
            A cutoff for filtering paired reads whose total distance to their
            respective restriction sites exceeds this value.
        """
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") %\
                                 (fendfilename.split('/')[-1]),
            return None
        self.fends = fendfilename
        self.maxinsert = maxinsert
        fends = h5py.File(fendfilename, 'r')
        chr2int = {}
        for i, j in enumerate(fends['chromosomes'][:]):
            chr2int[j] = i
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist, str):
            filelist = [filelist]
        total_reads = 0
        fend_pairs = {}
        for prefix in filelist:
            data = {}
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
                    # if chromosome not in chr2int, skip
                    chrom = input.getrname(read.tid).strip('chr')
                    if chrom not in chr2int:
                        continue
                    # skip multiply-aligned reads
                    for tag in read.tags:
                        if tag[0] == 'XS':
                            break
                    else:
                        if read.is_reverse:
                            strand = 1
                            end = str(read.pos + len(read.seq))
                        else:
                            strand = 0
                            end = str(read.pos)
                        unpaired[read.qname] = [chr2int[chrom], end, strand]
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
                    # if chromosome not in chr2int, skip
                    chrom = input.getrname(read.tid).strip('chr')
                    if chrom not in chr2int:
                        continue
                    # skip multiply-aligned reads
                    for tag in read.tags:
                        if tag[0] == 'XS':
                            break
                    else:
                        if read.is_reverse:
                            strand = 1
                            end = str(read.pos + len(read.seq))
                        else:
                            strand = 0
                            end = str(read.pos)
                        data[tuple([chr2int[chrom], end, strand] + unpaired[read.qname])] = 0
                        del unpaired[read.qname]
                input.close()
                print >> sys.stderr, ("Done\n"),
            data = numpy.array(data.keys(), dtype=numpy.int32)
            print >> sys.stderr, ("Read %i validly_mapped read paired.\n") % (data.shape[0]),
            total_reads += data.shape[0]
            # map data to fends, filtering as needed
            if data.shape[0] > 0:
                self._find_fend_pairs(data, fends, fend_pairs)
        if len(fend_pairs) == 0:
            print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
        print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(fend_pairs)),
        self._write_fend_pairs(fends, fend_pairs)
        return None

    def load_data_from_mat(self, fendfilename, filename, maxinsert=0):
        """
        load_data_from_mat method

        Read counts from mat file and place in h5dict.

        Parameters
        ----------
        fendfilename : string
            This specifies the filename of the fend object.
        filename : string
            Filename of mat file containing fend pair and count data.
        maxinsert : int, optional
            Cutoff value used in filtering total read lengths in mat file.
        """
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") %\
                                 (fendfilename.split('/')[-1]),
            return None
        self.fends = fendfilename
        self.maxinsert = maxinsert
        fends = h5py.File(fendfilename, 'r')
        # load data from mat file. This assumes that the mat data was mapped
        # using the same fend numbering as in the fend file.
        fend_pairs = {}
        if not os.path.exists(filename):
            print >> sys.stderr, ("%s not found... no data loaded.\n") % (filename.split('/')[-1]),
            return None
        print >> sys.stderr, ("Loading data from mat file..."),
        input = open(filename, 'r')
        for line in input:
            temp = line.strip('\n').split('\t')
            try:
                fend1 = int(temp[0]) - 1
            except ValueError:
                continue
            fend_pairs[(fend1, int(temp[1]) - 1)] = int(temp[2])
        input.close()
        if len(fend_pairs) == 0:
            print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
        # remove fend pairs from same fragment or opposite strand adjacents
        for i in range(fends['chromosomes'].shape[0]):
            for j in range(fends['chr_indices'][i], fends['chr_indices'][i + 1] - 2, 2):
                # same fend
                name = (j, j)
                if name in fend_pairs:
                    del fend_pairs[name]
                name = (j + 1, j + 1)
                if name in fend_pairs:
                    del fend_pairs[name]
                # same fragment
                name = (j, j + 1)
                if name in fend_pairs:
                    del fend_pairs[name]
                # adjacent fragments, opposite strands
                name = (j, j + 3)
                if name in fend_pairs:
                    del fend_pairs[name]
                name = (j + 1, j + 2)
                if name in fend_pairs:
                    del fend_pairs[name]
            j = fends['chr_indices'][i + 1] - 2
            # same fend
            name = (j, j)
            if name in fend_pairs:
                del fend_pairs[name]
            name = (j + 1, j + 1)
            if name in fend_pairs:
                del fend_pairs[name]
            # same fragment
            name = (j, j + 1)
            if name in fend_pairs:
                del fend_pairs[name]
        print >> sys.stderr, ("%i valid fend pairs loaded.\n") % (fend_pairs.shape[0]),
        # write fend pairs to h5dict
        self._write_fend_pairs(fends, fend_pairs)
        return None

    def _find_fend_pairs(self, data, fends, fend_pairs):
        """Return array with lower fend, upper fend, and count for pair."""
        mapped_fends = numpy.empty((data.shape[0], 2), dtype=numpy.int32)
        mapped_fends.fill(-1)
        distances = numpy.zeros((data.shape[0], 3), dtype=numpy.int32)
        # assign fends on a per-chromosome basis
        for i, chrom in enumerate(fends['chromosomes']):
            print >> sys.stderr, ("\rMapping first  fend for %s") % (chrom.ljust(10)),
            cuts = numpy.r_[fends['fends']['start'][fends['chr_indices'][i]:(fends['chr_indices'][i + 1])][::2],
                            fends['fends']['stop'][fends['chr_indices'][i + 1] - 1]]
            # determine valid first fends for chromosome
            where = numpy.where((data[:, 0] == i) * (data[:, 1] >= cuts[0]) * (data[:, 1] < cuts[-1]))[0]
            if where.shape[0] > 0:
                indices = numpy.searchsorted(cuts, data[where, 1], side='right')
                mapped_fends[where, 0] = (indices * 2 - 1 - data[where, 2]) + fends['chr_indices'][i]
                # find first distance from cutsite to mapped coordinate
                distances[where, 0] = numpy.abs(data[where, 1] - cuts[indices - data[where, 2]])
            print >> sys.stderr, ("\rMapping second fend for %s") % (chrom.ljust(10)),
            # determine valid second fends for chromosome
            where = numpy.where((data[:, 3] == i) * (data[:, 4] >= cuts[0]) * (data[:, 4] < cuts[-1]))[0]
            if where.shape[0] > 0:
                indices = numpy.searchsorted(cuts, data[where, 4], side='right')
                mapped_fends[where, 1] = (indices * 2 - 1 - data[where, 5]) + fends['chr_indices'][i]
                # find second distance from cutsite to mapped coordinate
                distances[where, 1] = numpy.abs(data[where, 4] - cuts[indices - data[where, 5]])
        print >> sys.stderr, ("\r%s\rCounting fend pairs...") % (' ' * 50),
        # arrange so first fend is always smaller number
        mapped_fends.sort(axis=1)
        # find validly mapped pairs
        distances[:, 2] = numpy.sum(distances[:, :2], axis=1)
        # remove reads not mapped to fends
        where = numpy.where((mapped_fends[:, 0] != -1) * (mapped_fends[:, 1] != -1) *
                            # remove reads with total insert size too large
                            (distances[:, 2] <= self.maxinsert))[0]
        # count pair occurences
        for i in where:
            name = (mapped_fends[i, 0], mapped_fends[i, 1])
            if name not in fend_pairs:
                fend_pairs[name] = 1
            else:
                fend_pairs[name] += 1
        # remove fend pairs from same fragment or opposite strand adjacents
        for i in range(fends['chromosomes'].shape[0]):
            for j in range(fends['chr_indices'][i], fends['chr_indices'][i + 1] - 2, 2):
                # same fend
                name = (j, j)
                if name in fend_pairs:
                    del fend_pairs[name]
                name = (j + 1, j + 1)
                if name in fend_pairs:
                    del fend_pairs[name]
                # same fragment
                name = (j, j + 1)
                if name in fend_pairs:
                    del fend_pairs[name]
                # adjacent fragments, opposite strands
                name = (j, j + 3)
                if name in fend_pairs:
                    del fend_pairs[name]
                name = (j + 1, j + 2)
                if name in fend_pairs:
                    del fend_pairs[name]
            j = fends['chr_indices'][i + 1] - 2
            # same fend
            name = (j, j)
            if name in fend_pairs:
                del fend_pairs[name]
            name = (j + 1, j + 1)
            if name in fend_pairs:
                del fend_pairs[name]
            # same fragment
            name = (j, j + 1)
            if name in fend_pairs:
                del fend_pairs[name]
        print >> sys.stderr, ("Done\n"),
        return None

    def _write_fend_pairs(self, fends, fend_pairs):
        """Separate fend pairs into cis and trans interactions and write to
        h5dict with index arrays."""
        print >> sys.stderr, ("Writing fend pair data to file..."),
        # convert counts into array
        fend_array = numpy.empty((len(fend_pairs), 3), dtype=numpy.int32)
        i = 0
        for name, count in fend_pairs.iteritems():
            fend_array[i, 0] = name[0]
            fend_array[i, 1] = name[1]
            fend_array[i, 2] = count
            i += 1
        del fend_pairs
        # reorder by first fend, then second fend
        fend_array = fend_array[numpy.lexsort((fend_array[:, 1], fend_array[:, 0])), :]
        # find which pairs are from the same chromosome
        cis = numpy.zeros(fend_array.shape[0], dtype=numpy.bool)
        for i in range(fends['chromosomes'].shape[0]):
            where = numpy.where((fend_array[:, 0] >= fends['chr_indices'][i]) *
                                (fend_array[:, 0] < fends['chr_indices'][i + 1]) *
                                (fend_array[:, 1] >= fends['chr_indices'][i]) *
                                (fend_array[:, 1] < fends['chr_indices'][i + 1]))[0]
            cis[where] = True
        # separate cis from trans data
        cis_data = fend_array[cis, :]
        trans_data = fend_array[numpy.invert(cis), :]
        if cis_data.shape[0] > 0:
            cis_indices = numpy.r_[0, numpy.bincount(cis_data[:, 0],
                                   minlength=fends['fends'].shape[0])].astype(numpy.int32)
            for i in range(1, cis_indices.shape[0]):
                cis_indices[i] += cis_indices[i - 1]
        if trans_data.shape[0] > 0:
            trans_indices = numpy.r_[0, numpy.bincount(trans_data[:, 0],
                                     minlength=fends['fends'].shape[0])].astype(numpy.int32)
            for i in range(1, trans_indices.shape[0]):
                trans_indices[i] += trans_indices[i - 1]
        self.data.attrs['maxinsert'] = self.maxinsert
        self.data.attrs['fendfilename'] = self.fends
        if cis_data.shape[0] > 0:
            self.data.create_dataset('cis_data', data=cis_data)
            self.data.create_dataset('cis_indices', data=cis_indices)
        if trans_data.shape[0] > 0:
            self.data.create_dataset('trans_data', data=trans_data)
            self.data.create_dataset('trans_indices', data=trans_indices)
        print >> sys.stderr, ("Done\n"),
        return None

    def export_to_mat(self, outfilename):
        """
        export_to_mat method

        Write reads loaded in data object to text file in mat format.

        Parameters
        ----------
        outfilename : string
            This specifies the file to save data to.
        """
        print >> sys.stderr, ("Writing data to mat file..."),
        output = open(outfilename, 'w')
        print >> output, "fend1\tfend2\tcount"
        cis_indices = self.data['cis_indices'][...]
        trans_indices = self.data['trans_indices'][...]
        cis_data = self.data['cis_data'][...]
        trans_data = self.data['trans_data'][...]
        for i in range(cis_indices.shape[0] - 1):
            for j in range(cis_indices[i], cis_indices[i + 1]):
                # One is added to indices so numbering starts from one.
                print >> output, "%i\t%i\t%i" % (i + 1, cis_data[j, 1] + 1, cis_data[j, 2])
            for j in range(trans_indices[i], trans_indices[i + 1]):
                print >> output, "%i\t%i\t%i" % (i + 1, trans_data[j, 1] + 1, trans_data[j, 2])
        output.close()
        print >> sys.stderr, ("Done\n"),
        return None
