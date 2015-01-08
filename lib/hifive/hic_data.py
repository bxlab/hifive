#!/usr/bin/env python

import os
import sys
from glob import glob

import numpy
import h5py
try:
    import pysam
except:
    pass

from fend import Fend


class HiCData(object):
    """This class handles interaction count data for HiC experiments.

    This class stores mapped paired-end reads, indexing them by fragment-end (fend) number, in an h5dict.

    .. note::
      This class is also available as hifive.HiCData

    When initialized, this class creates an h5dict in which to store all data associated with this object.
    
    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`HiCData` class object.
    """

    def __init__(self, filename, mode='r', silent=False):
        """
        Create a :class:`HiCData` object.
        """
        self.filename = os.path.abspath(filename)
        self.data = h5py.File(filename, mode)
        self.silent = silent
        return None

    def load_data_from_raw(self, fendfilename, filelist, maxinsert):
        """
        Read interaction counts from a text file(s) and place in h5dict.

        Files should contain both mapped ends of a read, one read per line, separated by tabs. Each line should be in the following format::

          chromosome1    coordinate1  strand1   chromosome2    coordinate2  strand2

        where strands are given by the characters '+' and '-'.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filelist: A list containing all of the file names of mapped read text files to be included in the dataset. If only one file is needed, this may be passed as a string.
        :type filelist: list
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :returns: None
        """
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            return None
        self.fends = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                os.path.dirname(self.filename)), os.path.basename(fendfilename))
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
                if not self.silent:
                    print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                continue
            if not self.silent:
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
            if not self.silent:
                print >> sys.stderr, ("%i validly-mapped reads loaded.\n") % (data.shape[0]),
            total_reads += data.shape[0]
            # map data to fends, filtering as needed
            if data.shape[0] > 0:
                self._find_fend_pairs(data, fends, fend_pairs)
        if len(fend_pairs) == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(fend_pairs)),
        # write fend pairs to h5dict
        self._write_fend_pairs(fends, fend_pairs)
        return None

    def load_data_from_bam(self, fendfilename, filelist, maxinsert):
        """
        Read interaction counts from pairs of BAM-formatted alignment file(s) and place in h5dict.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filelist: A list containing all of the bam file prefices to be included in the dataset. All files containing each prefix will be loaded. If only one pair of files is needed, the prefix may be passed as a string.
        :type filelist: list
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :returns: None
        """
        if 'pysam' not in sys.modules.keys():
            if not self.silent:
                print >> sys.stderr, ("The pysam module must be installed to use this function.")
            return None
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            return None
        self.fends = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                os.path.dirname(self.filename)), os.path.basename(fendfilename))
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
                if not self.silent:
                    print >> sys.stderr, ("No files found for prefix %s") % (prefix.split('/')[-1]),
                continue
            for fname in fnames:
                if not os.path.exists(fname.replace('_1.bam', '_2.bam')):
                    if not self.silent:
                        print >> sys.stderr, ("%s does not have a paired end.") % (fname.split('/')[-1]),
                    del fnames[fnames.index(fname)]
            unpaired = {}
            # load first half of paired ends
            for fname in fnames:
                if not self.silent:
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
                if not self.silent:
                    print >> sys.stderr, ("Done\n"),
            # load second half of paired ends
            for fname1 in fnames:
                fname = fname1.replace('_1.bam', '_2.bam')
                if not self.silent:
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
                if not self.silent:
                    print >> sys.stderr, ("Done\n"),
            data = numpy.array(data.keys(), dtype=numpy.int32)
            if not self.silent:
                print >> sys.stderr, ("Read %i validly_mapped read paired.\n") % (data.shape[0]),
            total_reads += data.shape[0]
            # map data to fends, filtering as needed
            if data.shape[0] > 0:
                self._find_fend_pairs(data, fends, fend_pairs)
        if len(fend_pairs) == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                                 (total_reads,len(fend_pairs)),
        self._write_fend_pairs(fends, fend_pairs)
        return None

    def load_data_from_mat(self, fendfilename, filename, maxinsert=0):
        """
        Read interaction counts from a :mod:`HiCPipe`-compatible 'mat' text file and place in h5dict.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filename: File name of a 'mat' file containing fend pair and interaction count data.
        :type filename: str.
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :returns: None
        """
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            return None
        self.fends = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                os.path.dirname(self.filename)), os.path.basename(fendfilename))
        self.maxinsert = maxinsert
        fends = h5py.File(fendfilename, 'r')
        # load data from mat file. This assumes that the mat data was mapped
        # using the same fend numbering as in the fend file.
        fend_pairs = {}
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("%s not found... no data loaded.\n") % (filename.split('/')[-1]),
            return None
        if not self.silent:
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
            if not self.silent:
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
        if not self.silent:
            print >> sys.stderr, ("%i valid fend pairs loaded.\n") % (len(fend_pairs)),
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
            if not self.silent:
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
            if not self.silent:
                print >> sys.stderr, ("\rMapping second fend for %s") % (chrom.ljust(10)),
            # determine valid second fends for chromosome
            where = numpy.where((data[:, 3] == i) * (data[:, 4] >= cuts[0]) * (data[:, 4] < cuts[-1]))[0]
            if where.shape[0] > 0:
                indices = numpy.searchsorted(cuts, data[where, 4], side='right')
                mapped_fends[where, 1] = (indices * 2 - 1 - data[where, 5]) + fends['chr_indices'][i]
                # find second distance from cutsite to mapped coordinate
                distances[where, 1] = numpy.abs(data[where, 4] - cuts[indices - data[where, 5]])
        if not self.silent:
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
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None

    def _write_fend_pairs(self, fends, fend_pairs):
        """Separate fend pairs into cis and trans interactions and write to
        h5dict with index arrays."""
        if not self.silent:
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
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None

    def export_to_mat(self, outfilename):
        """
        Write reads loaded in data object to text file in :mod:`HiCPipe`-compatible 'mat' format.

        :param outfilename: Specifies the file to save data in.
        :type outfilename: str.
        :returns: None
        """
        if not self.silent:
            print >> sys.stderr, ("Writing data to mat file..."),
        output = open(outfilename, 'w')
        if not self.silent:
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
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None

    def convert_to_binned(self, binnedfilename, outfilename):
        """
        Create new dataset file by reassigning reads to bins specified by a binned :class:`Fend` file.

        :param binnedfilename: Speficies the :class:`Fend` file to use to determine bin boundaries.
        :type binnedfilename: str.
        :param outfilename: Filename to write the new :class:`HiCData` object to.
        :type outfilename: str.
        :returns: None
        """
        # determine if binned fend file exists and if so, load it
        if not os.path.exists(binnedfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not found. No data was loaded.\n") % (binnedfilename),
            return None
        filename = os.path.abspath(self.filename)
        fendfilename = self.data.attrs['fendfilename']
        if fendfilename[:2] == './':
            fendfilename = fendfilename[2:]
        parent_count = fendfilename.count('../')
        fendfilename = '/'.join(os.path.abspath(filename).split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
        fends = h5py.File(fendfilename, 'r')
        bins = h5py.File(binnedfilename, 'r')
        outfile = h5py.File(outfilename, 'w')
        outfile.attrs['fendfilename'] = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(binnedfilename)),
                                                   os.path.dirname(self.filename)), os.path.basename(binnedfilename))
        outfile.attrs['original_fendfilename'] = self.data.attrs['fendfilename']
        outfile.attrs['maxinsert'] = self.data.attrs['maxinsert']
        outfile.attrs['binned'] = True
        chr2int = {}
        fend_chr2int = {}
        fend_chroms = fends['chromosomes'][...]
        bin_chroms = bins['chromosomes'][...]
        for i, j in enumerate(bin_chroms):
            if j in fend_chroms:
                chr2int[j] = i
                fend_chr2int[j] = numpy.where(fend_chroms == j)[0][0]
        mapping = numpy.zeros(fends['fends'].shape[0], dtype=numpy.int32) - 1
        bin_indices = bins['chr_indices'][...]
        fend_indices = fends['chr_indices'][...]
        for j in chr2int.keys():
            fstart = fend_indices[fend_chr2int[j]]
            fstop = fend_indices[fend_chr2int[j] + 1]
            bstart = bin_indices[chr2int[j]]
            bstop = bin_indices[chr2int[j] + 1]
            mapping[fstart:fstop] = numpy.searchsorted(bins['fends']['stop'][bstart:bstop],
                                                       fends['fends']['mid'][fstart:fstop],
                                                       side='right') + bstart
        cis_data = self.data['cis_data'][...]
        valid = numpy.where((mapping[cis_data[:, 0]] != -1) * (mapping[cis_data[:, 1]] != -1))[0]
        cis_data = cis_data[valid, :]
        cis_data[:, 0] = mapping[cis_data[:, 0]]
        cis_data[:, 1] = mapping[cis_data[:, 1]]
        valid = numpy.where(cis_data[:, 0] != cis_data[:, 1])[0]
        cis_data = cis_data[valid, :]
        indices = cis_data[:, 0].astype(numpy.int64) * bin_indices[-1] + cis_data[:, 1]
        unique_indices = numpy.unique(indices)
        index_mapping = numpy.searchsorted(unique_indices, indices)
        counts = numpy.bincount(index_mapping, weights=cis_data[:, 2], minlength=unique_indices.shape[0])
        del index_mapping
        cis_data = numpy.zeros((counts.shape[0], 3), dtype=numpy.int32)
        cis_data[:, 2] = counts
        cis_data[:, 0] = unique_indices / bin_indices[-1]
        cis_data[:, 1] = unique_indices % bin_indices[-1]
        del indices
        del unique_indices
        cis_indices = numpy.r_[0, numpy.bincount(cis_data[:, 0],
                               minlength=bin_indices[-1])].astype(numpy.int32)
        for i in range(cis_indices.shape[0] - 1):
            cis_indices[i + 1] += cis_indices[i]
        outfile.create_dataset(name='cis_data', data=cis_data)
        outfile.create_dataset(name='cis_indices', data=cis_indices)
        del cis_data
        del cis_indices
        trans_data = self.data['trans_data'][...]
        valid = numpy.where((mapping[trans_data[:, 0]] != -1) * (mapping[trans_data[:, 1]] != -1))[0]
        trans_data = trans_data[valid, :]
        trans_data[:, 0] = mapping[trans_data[:, 0]]
        trans_data[:, 1] = mapping[trans_data[:, 1]]
        indices = trans_data[:, 0].astype(numpy.int64) * bin_indices[-1] + trans_data[:, 1]
        unique_indices = numpy.unique(indices)
        index_mapping = numpy.searchsorted(unique_indices, indices)
        counts = numpy.bincount(index_mapping, weights=trans_data[:, 2])
        del index_mapping
        trans_data = numpy.zeros((counts.shape[0], 3), dtype=numpy.int32)
        trans_data[:, 2] = counts
        trans_data[:, 0] = unique_indices / bin_indices[-1]
        trans_data[:, 1] = unique_indices % bin_indices[-1]
        del indices
        del unique_indices
        trans_indices = numpy.bincount(trans_data[:, 0], minlength=(bins['fends'].shape[0] + 1))
        for i in range(trans_indices.shape[0] - 1):
            trans_indices[i + 1] += trans_indices[i]
        outfile.create_dataset(name='trans_data', data=trans_data)
        outfile.create_dataset(name='trans_indices', data=trans_indices)
        del trans_data
        del trans_indices
        outfile.close()
        fends.close()
        bins.close()
        return None
