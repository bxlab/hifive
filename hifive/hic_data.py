#!/usr/bin/env python

"""A class for handling HiC read data."""

import os
import sys

import numpy
import h5py
try:
    import pysam
except:
    pass


class HiCData(object):

    """This class handles interaction count data for HiC experiments.

    This class stores mapped paired-end reads, indexing them by fend-end (fend) number, in an h5dict.

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

    :Attributes: * **file** (*str.*) A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcomes.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a :class:`HiCData` object."""
        self.file = os.path.abspath(filename)
        self.silent = silent
        self.history = ''
        if mode != "w":
            self.load()
        return None

    def __getitem__(self, key):
        """Dictionary-like lookup."""
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        """Dictionary-like value setting."""
        self.__dict__[key] = value
        return None

    def save(self):
        """
        Save analysis parameters to h5dict.

        :returns: None
        """
        self.history.replace("'None'", "None")
        datafile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['file', 'chr2int', 'fends', 'silent']:
                continue
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                datafile.attrs[key] = self[key]
        datafile.close()
        return None

    def load(self):
        """
        Load data from h5dict specified at object creation.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :returns: None
        """
        datafile = h5py.File(self.file, 'r')
        for key in datafile.keys():
            self[key] = numpy.copy(datafile[key])
        for key in datafile['/'].attrs.keys():
            self[key] = datafile['/'].attrs[key]
        # ensure fend h5dict exists
        if 'fendfilename' in self.__dict__:
            fendfilename = self.fendfilename
            if fendfilename[:2] == './':
                fendfilename = fendfilename[2:]
            parent_count = fendfilename.count('../')
            fendfilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(fendfilename):
                if not self.silent:
                    print >> sys.stderr, ("Could not find %s. No fends loaded.\n") % (fendfilename),
            else:
                self.fends = h5py.File(fendfilename, 'r')
                # create dictionary for converting chromosome names to indices
                self.chr2int = {}
                for i, chrom in enumerate(self.fends['chromosomes']):
                    self.chr2int[chrom] = i
        datafile.close()
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

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_raw(fendfilename='%s', filelist=%s, maxinsert=%i) - " % (fendfilename, str(filelist), maxinsert)
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.maxinsert = maxinsert
        self.fends = h5py.File(fendfilename, 'r')
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        for i, j in enumerate(self.fends['chromosomes'][:]):
            self.chr2int[j] = i
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
                self.history += "'%s' not found, " % fname
                continue
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
            input = open(fname, 'r', 1)
            for line in input:
                temp = line.strip('\n').split('\t')
                if temp[0] not in self.chr2int or temp[3] not in self.chr2int:
                    continue
                data[(self.chr2int[temp[0]], int(temp[1]), strand[temp[2]],
                      self.chr2int[temp[3]], int(temp[4]), strand[temp[5]])] = 0
            input.close()
            data = numpy.array(data.keys(), dtype=numpy.int32)
            if not self.silent:
                print >> sys.stderr, ("%i validly-mapped reads pairs loaded.\n") % (data.shape[0]),
            new_pairs = data.shape[0]
            total_reads += new_pairs
            # map data to fends, filtering as needed
            if new_pairs > 0:
                self._find_fend_pairs(data, fend_pairs)
        self._clean_fend_pairs(fend_pairs)
        if len(fend_pairs) == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(fend_pairs)),
        # write fend pairs to h5dict
        self._parse_fend_pairs(fend_pairs)
        self.history += 'Success\n'
        return None

    def load_data_from_bam(self, fendfilename, filelist, maxinsert):
        """
        Read interaction counts from pairs of BAM-formatted alignment file(s) and place in h5dict.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filelist: A list containing lists of paired end bam files. If only one pair of files is needed, the list may contain both file path strings.
        :type filelist: list of mapped sequencing runs. Each run should be a list of the first and second read end bam files ([[run1_1, run1_2], [run2_1, run2_2]...])
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_bam(fendfilename='%s', filelist=%s, maxinsert=%i) - " % (fendfilename, str(filelist), maxinsert)
        if 'pysam' not in sys.modules.keys():
            if not self.silent:
                print >> sys.stderr, ("The pysam module must be installed to use this function.")
            self.history += 'Error: pysam module missing\n'
            return None
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.maxinsert = maxinsert
        self.fends = h5py.File(fendfilename, 'r')
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        for i, j in enumerate(self.fends['chromosomes'][:]):
            self.chr2int[j] = i
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist[0], str):
            filelist = [[filelist[0], filelist[1]]]
        total_reads = 0
        fend_pairs = {}
        for filepair in filelist:
            data = []
            # determine which files have both mapped ends present
            present = True
            if not os.path.exists(filepair[0]):
                if not self.silent:
                    print >> sys.stderr, ("%s could not be located.") % (filepair[0]),
                self.history += "'%s' not found, " % filepair[0]
                present = False
            if not os.path.exists(filepair[1]):
                if not self.silent:
                    print >> sys.stderr, ("%s could not be located.") % (filepair[1]),
                self.history += "'%s' not found, " % filepair[1]
                present = False
            if not present:
                if not self.silent:
                    print >> sys.stderr, ("No data for one or both ends could be located. Skipping this run.\n")
                continue
            unpaired = {}
            # load first half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[0].split('/')[-1]),
            input = pysam.Samfile(filepair[0], 'rb')
            idx2int = {}
            for i in range(len(input.header['SQ'])):
                chrom = input.header['SQ'][i]['SN'].strip('chr')
                if chrom in self.chr2int:
                    idx2int[i] = self.chr2int[chrom]
            for read in input.fetch(until_eof=True):
                # Only consider reads with an alignment
                if read.is_unmapped:
                    continue
                # if chromosome not in chr2int, skip
                if read.tid not in idx2int:
                    continue
                if read.is_reverse:
                    strand = 1
                    end = read.pos + len(read.seq)
                else:
                    strand = 0
                    end = read.pos
                unpaired[read.qname] = [idx2int[read.tid], end, strand]
            input.close()
            if not self.silent:
                print >> sys.stderr, ("Done\n"),
            # load second half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[1].split('/')[-1]),
            input = pysam.Samfile(filepair[1], 'rb')
            idx2int = {}
            for i in range(len(input.header['SQ'])):
                chrom = input.header['SQ'][i]['SN'].strip('chr')
                if chrom in self.chr2int:
                    idx2int[i] = self.chr2int[chrom]
            for read in input.fetch(until_eof=True):
                # Only consider reads with an alignment
                if read.is_unmapped:
                    continue
                # if chromosome not in chr2int, skip
                if read.tid not in idx2int:
                    continue
                if read.is_reverse:
                    strand = 1
                    end = read.pos + len(read.seq)
                else:
                    strand = 0
                    end = read.pos
                if read.qname in unpaired:
                    data.append(tuple([idx2int[read.tid], end, strand] + unpaired[read.qname]))
            input.close()
            del unpaired
            data = numpy.array(list(set(data)), dtype=numpy.int32)
            new_pairs = data.shape[0]
            if not self.silent:
                print >> sys.stderr, ("Done\nRead %i validly-mapped read pairs.\n") % (new_pairs),
            total_reads += new_pairs
            # map data to fends, filtering as needed
            if new_pairs > 0:
                self._find_fend_pairs(data, fend_pairs)
            del data
        self._clean_fend_pairs(fend_pairs)
        if len(fend_pairs) == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                                 (total_reads, len(fend_pairs)),
        self._parse_fend_pairs(fend_pairs)
        self.history += "Success\n"
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

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_mat(fendfilename='%s', filename='%s', maxinsert=%i) - " % (fendfilename, filename, maxinsert)
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.maxinsert = maxinsert
        self.fends = h5py.File(fendfilename, 'r')
        self.history = self.fends['/'].attrs['history'] + self.history
        # load data from mat file. This assumes that the mat data was mapped
        # using the same fend numbering as in the fend file.
        fend_pairs = {}
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("%s not found... no data loaded.\n") % (filename.split('/')[-1]),
            self.history += "Error: '%s' not found\n" % fendfilename
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
        self._clean_fend_pairs(fend_pairs)
        if len(fend_pairs) == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i valid fend pairs loaded.\n") % (len(fend_pairs)),
        # write fend pairs to h5dict
        self._parse_fend_pairs(fend_pairs)
        self.history += "Success\n"
        return None

    def _find_fend_pairs(self, data, fend_pairs):
        """Return array with lower fend, upper fend, and count for pair."""
        # split data amongst nodes
        mapped_fends = numpy.empty((data.shape[0], 2), dtype=numpy.int32)
        mapped_fends.fill(-1)
        distances = numpy.zeros((data.shape[0], 2), dtype=numpy.int32)
        # assign fends on a per-chromosome basis
        for i, chrom in enumerate(self.fends['chromosomes']):
            if not self.silent:
                print >> sys.stderr, ("\rMapping first  fend for %s") % (chrom.ljust(10)),
            cuts = numpy.r_[
                self.fends['fends']['start'][self.fends['chr_indices'][i]:(self.fends['chr_indices'][i + 1])][::2],
                self.fends['fends']['stop'][self.fends['chr_indices'][i + 1] - 1]]
            # determine valid first fends for chromosome
            temp = numpy.where(data[:, 0] == i)[0]
            where = temp[numpy.where((data[temp, 1] >= cuts[0]) * (data[temp, 1] < cuts[-1]))]
            if where.shape[0] > 0:
                indices = numpy.searchsorted(cuts, data[where, 1], side='right')
                mapped_fends[where, 0] = (indices * 2 - 1 - data[where, 2]) + self.fends['chr_indices'][i]
                # find first distance from cutsite to mapped coordinate
                distances[where, 0] = numpy.abs(data[where, 1] - cuts[indices - data[where, 2]])
            if not self.silent:
                print >> sys.stderr, ("\rMapping second fend for %s") % (chrom.ljust(10)),
            # determine valid second fends for chromosome
            temp = numpy.where(data[:, 3] == i)[0]
            where = temp[numpy.where((data[temp, 4] >= cuts[0]) * (data[temp, 4] < cuts[-1]))]
            if where.shape[0] > 0:
                indices = numpy.searchsorted(cuts, data[where, 4], side='right')
                mapped_fends[where, 1] = (indices * 2 - 1 - data[where, 5]) + self.fends['chr_indices'][i]
                # find second distance from cutsite to mapped coordinate
                distances[where, 1] = numpy.abs(data[where, 4] - cuts[indices - data[where, 5]])
        # remove fends with too great an insert distance
        distances[:, 0] += distances[:, 1]
        where = numpy.where(distances[:, 0] > self.maxinsert)[0]
        mapped_fends[where, :] = -1
        if not self.silent:
            print >> sys.stderr, ("\r%s\rCounting fend pairs...") % (' ' * 50),
        # arrange so first fend is always smaller number
        mapped_fends.sort(axis=1)
        # remove reads not mapped to fends
        where = numpy.where((mapped_fends[:, 0] != -1) * (mapped_fends[:, 1] != -1))[0]
        # count pair occurences
        for i in where:
            name = (mapped_fends[i, 0], mapped_fends[i, 1])
            if name not in fend_pairs:
                fend_pairs[name] = 1
            else:
                fend_pairs[name] += 1
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None

    def _clean_fend_pairs(self, fend_pairs):
        # remove fend pairs from same fend or opposite strand adjacents
        for i in range(self.fends['chromosomes'].shape[0]):
            for j in range(self.fends['chr_indices'][i], self.fends['chr_indices'][i + 1] - 2, 2):
                # same fend
                name = (j, j)
                if name in fend_pairs:
                    del fend_pairs[name]
                name = (j + 1, j + 1)
                if name in fend_pairs:
                    del fend_pairs[name]
                # same fend
                name = (j, j + 1)
                if name in fend_pairs:
                    del fend_pairs[name]
                # adjacent fends, opposite strands
                name = (j, j + 3)
                if name in fend_pairs:
                    del fend_pairs[name]
                name = (j + 1, j + 2)
                if name in fend_pairs:
                    del fend_pairs[name]
            j = self.fends['chr_indices'][i + 1] - 2
            # same fend
            name = (j, j)
            if name in fend_pairs:
                del fend_pairs[name]
            name = (j + 1, j + 1)
            if name in fend_pairs:
                del fend_pairs[name]
            # same fend
            name = (j, j + 1)
            if name in fend_pairs:
                del fend_pairs[name]
        return None

    def _parse_fend_pairs(self, fend_pairs):
        """Separate fend pairs into cis and trans interactions index."""
        if not self.silent:
            print >> sys.stderr, ("Parsing fend pairs..."),
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
        for i in range(self.fends['chromosomes'].shape[0]):
            temp = numpy.where((fend_array[:, 0] >= self.fends['chr_indices'][i]) *
                               (fend_array[:, 0] < self.fends['chr_indices'][i + 1]))[0]
            where = temp[numpy.where((fend_array[temp, 1] >= self.fends['chr_indices'][i]) *
                                     (fend_array[temp, 1] < self.fends['chr_indices'][i + 1]))]
            cis[where] = True
        # separate cis from trans data
        self.cis_data = fend_array[cis, :]
        self.trans_data = fend_array[numpy.invert(cis), :]
        if self.cis_data.shape[0] > 0:
            cis_indices = numpy.r_[0, numpy.bincount(self.cis_data[:, 0],
                                   minlength=self.fends['fends'].shape[0])].astype(numpy.int64)
            for i in range(1, cis_indices.shape[0]):
                cis_indices[i] += cis_indices[i - 1]
            self.cis_indices = cis_indices
        if self.trans_data.shape[0] > 0:
            trans_indices = numpy.r_[0, numpy.bincount(self.trans_data[:, 0],
                                     minlength=self.fends['fends'].shape[0])].astype(numpy.int64)
            for i in range(1, trans_indices.shape[0]):
                trans_indices[i] += trans_indices[i - 1]
            self.trans_indices = trans_indices
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
        for i in range(self.cis_indices.shape[0] - 1):
            for j in range(self.cis_indices[i], self.cis_indices[i + 1]):
                # One is added to indices so numbering starts from one.
                print >> output, "%i\t%i\t%i" % (i + 1, self.cis_data[j, 1] + 1, self.cis_data[j, 2])
            for j in range(self.trans_indices[i], self.trans_indices[i + 1]):
                print >> output, "%i\t%i\t%i" % (i + 1, self.trans_data[j, 1] + 1, self.trans_data[j, 2])
        output.close()
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None
