#!/usr/bin/env python

"""A class for handling 5C read data."""

import os
import sys

import numpy
import h5py
try:
    import pysam
except:
    pass



class FiveCData(object):

    """
    This class handles interaction count data for 5C experiments.

    The FiveCData class contains all of the interaction information for a 5C experiment, including pairs of fragment indices and their associated counts.
   
    .. note::
      This class is also available as hifive.FiveCData

    When initialized, this class creates an h5dict in which to store all data associated with this object.
    
    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`FiveCData` class object.

    :Attributes: * **file** (*str.*) A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcomes.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a :class:`FiveCData` object."""
        self.file = os.path.abspath(filename)
        self.silent = silent
        self.history = ''
        self.filetype = 'fivec_data'
        if mode != 'w':
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
            if key in ['file', 'chr2int', 'frags', 'silent']:
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
        if 'fragfilename' in self.__dict__:
            fragfilename = self.fragfilename
            if fragfilename[:2] == './':
                fragfilename = fragfilename[2:]
            parent_count = fragfilename.count('../')
            fragfilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                fragfilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(fragfilename):
                if not self.silent:
                    print >> sys.stderr, ("Could not find %s. No fragments loaded.\n") % (fragfilename),
            else:
                self.frags = h5py.File(fragfilename, 'r')
                # create dictionary for converting chromosome names to indices
                self.chr2int = {}
                for i, chrom in enumerate(self.frags['chromosomes']):
                    self.chr2int[chrom] = i
        datafile.close()
        return None

    def load_data_from_counts(self, fragfilename, filelist):
        """
        Read interaction counts from a text file(s) and place in h5dict.

        :param fragfilename: This specifies the file name of the :class:`Fragment` object to associate with the dataset.
        :type fragfilename: str.
        :param filelist: A list containing all of the file names of counts text files to be included in the dataset. If only one file is needed, this may be passed as a string.
        :type filelist: list
        :returns: None

        :Attributes: * **fragfilename** (*str.*) - A string containing the relative path of the fragment file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-regional fragment pairings observed in the data. The first column contains the fragment index (from the 'fragments' array in the Fragment object) of the upstream fragment, the second column contains the idnex of the downstream fragment, and the third column contains the number of reads observed for that fragment pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fragments + 1. Each position contains the first entry for the correspondingly-indexed fragment in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fragment at index 5 in the Fragment object 'fragments' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-regional fragment pairings observed in the data. The first column contains the fragment index (from the 'fragments' array in the Fragment object) of the upstream fragment (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fragment, and the third column contains the number of reads observed for that fragment pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fragments + 1. Each position contains the first entry for the correspondingly-indexed fragment in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fragment at index 5 in the Fragment object 'fragments' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*ndarray*) - A filestream to the hdf5 Fragment file such that all saved Fragment attributes can be accessed through this class attribute.

        When data is loaded the 'history' attribute is updated to include the history of the Fragment file that becomes associated with it.
        """
        self.history += "FiveCData.load_data_from_counts(fragfilename='%s', filelist=%s) - " % (fragfilename, str(filelist))
        # determine if fragment file exists and if so, load it
        if not os.path.exists(fragfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fragment file %s was not found. No data was loaded.\n") % (fragfilename),
            self.history += "Error: '%s' not found\n" % fragfilename
            return None
        self.fragfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fragfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fragfilename))
        self.frags = h5py.File(fragfilename, 'r')
        self.history = self.frags['/'].attrs['history'] + self.history
        strands = self.frags['fragments']['strand'][...]
        chr2int = {}
        for i, j in enumerate(self.frags['chromosomes'][:]):
            chr2int[j] = i
        # create fragment name dictionary
        names = {}
        for i in range(self.frags['fragments'].shape[0]):
            names[self.frags['fragments']['name'][i]] = i
        # load data from all files, skipping if name not in the fragment file.
        if isinstance(filelist, str):
            filelist = [filelist]
        total_reads = 0
        data = {}
        for fname in filelist:
            reads = 0
            if not os.path.exists(fname):
                if not self.silent:
                    print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                self.history += "Error: '%s' not found, " % fname
                continue
            if not self.silent:
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
            if not self.silent:
                print >> sys.stderr, ("%i validly-mapped reads loaded.\n") % (reads),
            total_reads += reads
        if len(data) == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            return None
            self.history += "Error: no valid data loaded\n"
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(data)),
        # write fragment pairs to h5dict
        self._parse_fragment_pairs(data)
        self.history += "Success\n"
        return None

    def load_data_from_bam(self, fragfilename, filelist):
        """
        Read interaction counts from pairs of BAM files and place in h5dict.

        :param fragfilename: This specifies the file name of the :class:`Fragment` object to associate with the dataset.
        :type fragfilename: str.
        :param filelist: A list containing lists of paired read end files.
        :type filelist: list
        :returns: None

        :Attributes: * **fragfilename** (*str.*) - A string containing the relative path of the fragment file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-regional fragment pairings observed in the data. The first column contains the fragment index (from the 'fragments' array in the Fragment object) of the upstream fragment, the second column contains the idnex of the downstream fragment, and the third column contains the number of reads observed for that fragment pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fragments + 1. Each position contains the first entry for the correspondingly-indexed fragment in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fragment at index 5 in the Fragment object 'fragments' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-regional fragment pairings observed in the data. The first column contains the fragment index (from the 'fragments' array in the Fragment object) of the upstream fragment (upstream also refers to the lower indexed chromosome in this context), the second column contains the idnex of the downstream fragment, and the third column contains the number of reads observed for that fragment pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fragments + 1. Each position contains the first entry for the correspondingly-indexed fragment in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fragment at index 5 in the Fragment object 'fragments' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **frags** (*filestream*) - A filestream to the hdf5 Fragment file such that all saved Fragment attributes can be accessed through this class attribute.

        When data is loaded the 'history' attribute is updated to include the history of the Fragment file that becomes associated with it.
        """
        self.history += "FiveCData.load_data_from_counts(fragfilename='%s', filelist=%s) - " % (fragfilename, str(filelist))
        if 'pysam' not in sys.modules.keys():
            if not self.silent:
                print >> sys.stderr, ("The pysam module must be installed to use this function.")
            self.history += 'Error: pysam module missing\n'
            return None
        # determine if fragment file exists and if so, load it
        if not os.path.exists(fragfilename):
            if not self.silent:
                print >> sys.stderr, ("The fragment file %s was not found. No data was loaded.\n") % (fragfilename),
            self.history += "Error: '%s' not found\n" % fragfilename
            return None
        self.fragfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fragfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fragfilename))
        self.frags = h5py.File(fragfilename, 'r')
        strands = self.frags['fragments']['strand'][...]
        chr2int = {}
        for i, j in enumerate(self.frags['chromosomes'][:]):
            chr2int[j] = i
        # create fragment name dictionary
        names = {}
        for i in range(self.frags['fragments'].shape[0]):
            names[self.frags['fragments']['name'][i]] = i
        # load data from all files, skipping if either fragment not not in the fragment file.
        if isinstance(filelist[0], str):
            filelist = [[filelist[0], filelist[1]]]
        total_reads = 0
        data = {}
        for filepair in filelist:
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
            reads = 0
            unpaired = {}
            # load first half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[0].split('/')[-1]),
            input = pysam.Samfile(filepair[0], 'rb')
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
            if not self.silent:
                print >> sys.stderr, ("Done\n"),
            # load second half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[1].split('/')[-1]),
            input = pysam.Samfile(filepair[1], 'rb')
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
            if not self.silent:
                print >> sys.stderr, ("Done\n"),
            if not self.silent:
                print >> sys.stderr, ("Read %i validly_mapped read paired.\n") % (reads),
            total_reads += reads           
        if len(data) == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i unique pairs\n") %\
                             (total_reads,len(data)),
        self._parse_fragment_pairs(data)
        self.history += 'Success\n'
        return None

    def _parse_fragment_pairs(self, frag_pairs):
        """Separate frag pairs into cis (within region) and trans (between region) interactions and write to h5dict with index arrays."""
        if not self.silent:
            print >> sys.stderr, ("Writing fragment pair data to file..."),
        cis = {}
        trans = {}
        for i in range(self.frags['regions'].shape[0]):
            # Find pairs from same region
            for j in range(self.frags['regions']['start_frag'][i], self.frags['regions']['stop_frag'][i] - 1):
                for k in range(j + 1, self.frags['regions']['stop_frag'][i]):
                    if (j, k) in frag_pairs:
                        cis[(j, k)] = frag_pairs[(j, k)]
            # Find pairs from different regions
            for j in range(self.frags['regions']['start_frag'][i], self.frags['regions']['stop_frag'][i]):
                for k in range(self.frags['regions']['stop_frag'][i], self.frags['fragments'].shape[0]):
                    if (j, k) in frag_pairs:
                        trans[(j, k)] = frag_pairs[(j, k)]
        # convert data into arrays
        self.cis_data = numpy.empty((len(cis), 3), dtype=numpy.int32)
        self.trans_data = numpy.empty((len(trans), 3), dtype=numpy.int32)
        keys = cis.keys()
        keys.sort()
        for i in range(len(keys)):
            self.cis_data[i, 0] = keys[i][0]
            self.cis_data[i, 1] = keys[i][1]
            self.cis_data[i, 2] = cis[keys[i]]
        keys = trans.keys()
        keys.sort()
        for i in range(len(keys)):
            self.trans_data[i, 0] = keys[i][0]
            self.trans_data[i, 1] = keys[i][1]
            self.trans_data[i, 2] = trans[keys[i]]
        # find first instance of each fragment for cis and trans data
        if self.cis_data.shape[0] > 0:
            self.cis_indices = numpy.r_[0, numpy.bincount(self.cis_data[:, 0],
                                        minlength=self.frags['fragments'].shape[0])].astype(numpy.int64)
            for i in range(1, self.cis_indices.shape[0]):
                self.cis_indices[i] += self.cis_indices[i - 1]
        if self.trans_data.shape[0] > 0:
            self.trans_indices = numpy.r_[0, numpy.bincount(self.trans_data[:, 0],
                                          minlength=self.frags['fragments'].shape[0])].astype(numpy.int64)
            for i in range(1, self.trans_indices.shape[0]):
                self.trans_indices[i] += self.trans_indices[i - 1]
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None
