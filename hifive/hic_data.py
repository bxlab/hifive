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
        self.filetype = 'hic_data'
        self.stats = {
            'total_reads': 0,
            'valid_cis_reads': 0,
            'valid_trans_reads': 0,
            'valid_cis_pairs': 0,
            'valid_trans_pairs': 0,
            'pcr_duplicates': 0,
            'out_of_bounds': 0,
            'insert_size': 0,
            'same_fragment': 0,
            'failed_cut': 0,
            'chr_not_in_fends': 0,
        }
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
        chroms = self.fends['chromosomes'][...]
        for key in self.__dict__.keys():
            if key in ['file', 'chr2int', 'fends', 'silent', 'cuts']:
                continue
            elif key == 'stats':
                stats = []
                for name, count in self.stats.iteritems():
                    stats.append((name, count))
                stats = numpy.array(stats, dtype=numpy.dtype([('name', 'S20'), ('count', numpy.int32)]))
                datafile.create_dataset(name='stats', data=stats)
            elif self[key] is None:
                continue
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif isinstance(self[key], list):
                if isinstance(self[key][0], numpy.ndarray):
                    for i in range(len(self[key])):
                        datafile.create_dataset("%s.%s" % (key, chroms[i]), data=self[key][i])
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
            if key == 'stats':
                stats = datafile[key][...]
                for i in range(stats.shape[0]):
                    self.stats[stats['name'][i]] = stats['count'][i]
            else:
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

    def _find_cut_sites(self):
        self.cuts = []
        chroms = self.fends['chromosomes'][...]
        chr_indices = self.fends['chr_indices'][...]
        for i in range(len(chroms)):
            self.cuts.append(numpy.r_[
                self.fends['fends']['start'][chr_indices[i]:(chr_indices[i + 1])][::2],
                self.fends['fends']['stop'][chr_indices[i + 1] - 1]])
        return None

    def load_data_from_raw(self, fendfilename, filelist, maxinsert, skip_duplicate_filtering=False):
        """
        Read interaction counts from a text file(s) and place in h5dict.

        Files should contain both mapped ends of a read, one read per line, separated by tabs. Each line should be in the following format::

          chromosome1    coordinate1  strand1   chromosome2    coordinate2  strand2

        where strands are given by the characters '+' and '-'.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filelist: A list containing all of the file names of mapped read text files to be included in the dataset. If only one file is needed, this may be passed as a string.
        :type filelist: list
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value. If data was produced without a restriction enzyme (fend object has no fend data, only bin data), this integer specifies the maximum intra-chromosomal insert size that strandedness is considered for filtering. Fragments below the maxinsert size are only kept if they occur on the same orientation strand. This filtering is skipped is maxinsert is None.
        :type maxinsert: int.
        :param skip_duplicate_filtering: Do not remove PCR duplicates. This allows much lower memoer requirements since files can be processed in chunks.
        :type skip_duplicate_filtering: bool.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **fends** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_raw(fendfilename='%s', filelist=%s, maxinsert=%i, skip_duplicate_filtering=%s) - " % (fendfilename, str(filelist), maxinsert, str(skip_duplicate_filtering))
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
        if 'binned' in self.fends['/'].attrs and self.fends['/'].attrs['binned'] is not None:
            self.binned = True
        else:
            self.binned = False
        if 'fends' in self.fends and self.fends['fends'] is not None:
            self.re = True
        else:
            self.re = False
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        chroms = self.fends['chromosomes'][...]
        for i, j in enumerate(chroms):
            self.chr2int[j] = i
        self.insert_distribution = numpy.zeros((182, 2), dtype=numpy.int32)
        self.insert_distribution[1:, 1] = numpy.round(numpy.exp(numpy.linspace(3.8, 12.8, 181))).astype(numpy.int32)
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist, str):
            filelist = [filelist]
        raw_filelist = []
        for filename in filelist:
            raw_filelist.append("%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                           os.path.dirname(self.file)), os.path.basename(filename)))
        self.raw_filelist = ",".join(raw_filelist)
        fend_pairs = []
        for i in range(len(chroms)):
            fend_pairs.append([])
            for j in range(i + 1):
                fend_pairs[i].append({})
        total_reads = 0
        for fname in filelist:
            data = []
            for i in range(len(chroms)):
                data.append([])
                for j in range(i + 1):
                    data[i].append({})
            if not os.path.exists(fname):
                if not self.silent:
                    print >> sys.stderr, ("The file %s was not found...skipped.\n") % (fname.split('/')[-1]),
                self.history += "'%s' not found, " % fname
                continue
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (fname.split('/')[-1]),
            a = 0
            input = open(fname, 'r', 1)
            new_reads = 0
            for line in input:
                temp = line.strip('\n').split('\t')
                if temp[0] not in self.chr2int or temp[3] not in self.chr2int:
                    self.stats['chr_not_in_fends'] += 1
                    continue
                chrint1 = self.chr2int[temp[0]]
                chrint2 = self.chr2int[temp[3]]
                start1 = int(temp[1])
                start2 = int(temp[4])
                if temp[2] == '-':
                    start1 = -start1
                if temp[5] == '-':
                    start2 = -start2
                if chrint1 == chrint2:
                    if abs(start1) > abs(start2):
                        start1, start2 = start2, start1
                elif chrint2 < chrint1:
                    chrint1, chrint2, start1, start2 = chrint2, chrint1, start2, start1
                key = (start1, start2)
                if key not in data[chrint2][chrint1]:
                    if skip_duplicate_filtering:
                        data[chrint2][chrint1][key] = 1
                    else:
                        data[chrint2][chrint1][key] = None
                else:
                    if skip_duplicate_filtering:
                        data[chrint2][chrint1][key] += 1
                    else:
                        self.stats['pcr_duplicates'] += 1
                a += 1
                if skip_duplicate_filtering and a >= 10000000:
                    for i in range(len(data)):
                        for j in range(len(data[i])):
                            temp = numpy.zeros((len(data[i][j]), 3), dtype=numpy.int32)
                            k = 0
                            for key, count in data[i][j].iteritems():
                                temp[k, :] = (key[0], key[1], count)
                                self.stats['pcr_duplicates'] += count - 1
                                k += 1
                            data[i][j] = temp
                            new_reads += numpy.sum(data[i][j][:, 2])
                    if self.re:
                        self._find_fend_pairs(data, fend_pairs, skip_duplicate_filtering)
                    else:
                        self._find_bin_pairs(data, fend_pairs, skip_duplicate_filtering)
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rLoading data from %s...") % (' '*50, fname.split('/')[-1]),
                    for i in range(len(data)):
                        for j in range(len(data[i])):
                            data[i][j] = {}
                    a = 0
            input.close()
            if not skip_duplicate_filtering or a > 0:
                for i in range(len(data)):
                    for j in range(len(data[i])):
                        if skip_duplicate_filtering:
                            temp = numpy.zeros((len(data[i][j]), 3), dtype=numpy.int32)
                            k = 0
                            for key, count in data[i][j].iteritems():
                                temp[k, :] = (key[0], key[1], count)
                                self.stats['pcr_duplicates'] += count - 1
                                k += 1
                            data[i][j] = temp
                            new_reads += numpy.sum(data[i][j][:, 2])
                        else:
                            data[i][j] = numpy.array(data[i][j].keys(), dtype=numpy.int32)
                            new_reads += data[i][j].shape[0]
                # map data to fends, filtering as needed
                if new_reads > 0 and (a > 0 or not skip_duplicate_filtering):
                    if self.re:
                        self._find_fend_pairs(data, fend_pairs, skip_duplicate_filtering)
                    else:
                        self._find_bin_pairs(data, fend_pairs, skip_duplicate_filtering)
            total_reads += new_reads
            if not self.silent:
                print >> sys.stderr, ("\r%s\r%i validly-mapped reads pairs loaded.\n") % (' ' * 50, new_reads),
        if skip_duplicate_filtering:
            self.stats['total_reads'] = total_reads + self.stats['chr_not_in_fends']
        else:
            self.stats['total_reads'] = total_reads + self.stats['chr_not_in_fends'] + self.stats['pcr_duplicates']
        if self.re:
            self._clean_fend_pairs(fend_pairs)
        total_fend_pairs = 0
        for i in range(len(fend_pairs)):
            for j in range(len(fend_pairs[i])):
                total_fend_pairs += len(fend_pairs[i][j])
        if total_fend_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i valid fend pairs\n") %\
                             (total_reads, total_fend_pairs),
        # write fend pairs to h5dict
        if self.re and self.binned:
            self._parse_binned_fend_pairs(fend_pairs)
        else:
            self._parse_fend_pairs(fend_pairs)
        self.history += 'Success\n'
        return None

    def load_data_from_bam(self, fendfilename, filelist, maxinsert, skip_duplicate_filtering=False):
        """
        Read interaction counts from pairs of BAM-formatted alignment file(s) and place in h5dict.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filelist: A list containing lists of paired end bam files. If only one pair of files is needed, the list may contain both file path strings.
        :type filelist: list of mapped sequencing runs. Each run should be a list of the first and second read end bam files ([[run1_1, run1_2], [run2_1, run2_2]...])
        :param maxinsert: A cutoff for filtering paired end reads whose total distance to their respective restriction sites exceeds this value.
        :type maxinsert: int.
        :param skip_duplicate_filtering: Do not remove PCR duplicates. This allows much lower memoer requirements since files can be processed in chunks.
        :type skip_duplicate_filtering: bool.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **fends** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.
                     * **maxinsert** (*int.*) - An interger denoting the maximum included distance sum between both read ends and their downstream RE site.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_bam(fendfilename='%s', filelist=%s, maxinsert=%i, skip_duplicate_filtering=%s) - " % (fendfilename, str(filelist), maxinsert, str(skip_duplicate_filtering))
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
        if 'binned' in self.fends['/'].attrs and self.fends['/'].attrs['binned'] is not None:
            self.binned = True
        else:
            self.binned = False
        if 'fends' in self.fends and self.fends['fends'] is not None:
            self.re = True
        else:
            self.re = False
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        chroms = self.fends['chromosomes'][...]
        for i, j in enumerate(chroms):
            self.chr2int[j] = i
        self.insert_distribution = numpy.zeros((92, 2), dtype=numpy.int32)
        self.insert_distribution[1:, 1] = numpy.round(numpy.exp(numpy.linspace(3.8, 8.3, 91))).astype(numpy.int32)
        # load data from all files, skipping if chromosome not in the fend file.
        if isinstance(filelist[0], str):
            filelist = [[filelist[0], filelist[1]]]
        bam_filelist = []
        for filenames in filelist:
            bam_filelist.append("%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filenames[0])),
                                           os.path.dirname(self.file)), os.path.basename(filenames[0])))
            bam_filelist.append("%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filenames[1])),
                                           os.path.dirname(self.file)), os.path.basename(filenames[1])))
        self.bam_filelist = ",".join(bam_filelist)
        total_reads = 0
        fend_pairs = []
        for i in range(len(chroms)):
            fend_pairs.append([])
            for j in range(i + 1):
                fend_pairs[i].append({})
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
                continue
            unpaired = {}
            # load first half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[0].split('/')[-1]),
            input = pysam.Samfile(filepair[0], 'rb')
            idx2int = {}
            for i in range(len(input.header['SQ'])):
                chrom = input.header['SQ'][i]['SN']
                if chrom in self.chr2int:
                    idx2int[i] = self.chr2int[chrom]
            for read in input.fetch(until_eof=True):
                # Only consider reads with an alignment
                if read.is_unmapped:
                    continue
                # if chromosome not in chr2int, skip
                if read.tid not in idx2int:
                    self.stats['chr_not_in_fends'] += 1
                    continue
                if read.is_reverse:
                    end = -(read.pos + len(read.seq))
                else:
                    end = read.pos
                unpaired[read.qname] = (idx2int[read.tid], end)
            input.close()
            if not self.silent:
                print >> sys.stderr, ("Done\n"),
            # load second half of paired ends
            if not self.silent:
                print >> sys.stderr, ("Loading data from %s...") % (filepair[1].split('/')[-1]),
            data = []
            for i in range(len(chroms)):
                data.append([])
                for j in range(i + 1):
                    data[i].append({})
            input = pysam.Samfile(filepair[1], 'rb')
            idx2int = {}
            for i in range(len(input.header['SQ'])):
                chrom = input.header['SQ'][i]['SN']
                if chrom in self.chr2int:
                    idx2int[i] = self.chr2int[chrom]
            a = 0
            new_reads = 0
            for read in input.fetch(until_eof=True):
                # Only consider reads with an alignment
                if read.is_unmapped:
                    continue
                if read.qname not in unpaired:
                    continue
                # if chromosome not in chr2int, skip
                if read.tid not in idx2int:
                    self.stats['chr_not_in_fends'] += 1
                    continue
                if read.is_reverse:
                    start2 = -(read.pos + len(read.seq))
                else:
                    start2 = read.pos
                chr1, start1 = unpaired[read.qname]
                chr2 = idx2int[read.tid]
                if chr1 == chr2:
                    if abs(start1) > abs(start2):
                        start1, start2 = start2, start1
                elif chr2 < chr1:
                    chr1, chr2, start1, start2 = chr2, chr1, start2, start1
                key = (start1, start2)
                if key not in data[chr2][chr1]:
                    if skip_duplicate_filtering:
                        data[chr2][chr1][key] = 1
                    else:
                        data[chr2][chr1][key] = None
                else:
                    if skip_duplicate_filtering:
                        data[chr2][chr1][key] += 1
                    else:
                        self.stats['pcr_duplicates'] += 1
                a += 1
                if skip_duplicate_filtering and a >= 10000000:
                    for i in range(len(data)):
                        for j in range(len(data[i])):
                            temp = numpy.zeros((len(data[i][j]),3), dtype=numpy.int32)
                            k = 0
                            for key, count in data[i][j].iteritems():
                                temp[k, :] = (key[0], key[1], count)
                                self.stats['pcr_duplicates'] += count - 1
                                k += 1
                            data[i][j] = temp
                            new_reads += numpy.sum(data[i][j][:, 2])
                    if self.re:
                        self._find_fend_pairs(data, fend_pairs, skip_duplicate_filtering)
                    else:
                        self._find_bin_pairs(data, fend_pairs, skip_duplicate_filtering)
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rLoading data from %s...") % (' '*50, filepair[1].split('/')[-1]),
                    for i in range(len(data)):
                        for j in range(len(data[i])):
                            data[i][j] = {}
                    a = 0
            input.close()
            del unpaired
            if a > 0 or not skip_duplicate_filtering:
                for i in range(len(data)):
                    for j in range(len(data[i])):
                        if skip_duplicate_filtering:
                                temp = numpy.zeros((len(data[i][j]),3), dtype=numpy.int32)
                                k = 0
                                for key, count in data[i][j].iteritems():
                                    temp[k, :] = (key[0], key[1], count)
                                    self.stats['pcr_duplicates'] += count - 1
                                    k += 1
                                data[i][j] = temp
                                new_reads += numpy.sum(data[i][j][:, 2])
                        else:
                            data[i][j] = numpy.array(data[i][j].keys(), dtype=numpy.int32)
                            new_reads += data[i][j].shape[0]
                if new_reads > 0 and (a > 0 or not skip_duplicate_filtering):
                    if self.re:
                        self._find_fend_pairs(data, fend_pairs, skip_duplicate_filtering)
                    else:
                        self._find_bin_pairs(data, fend_pairs, skip_duplicate_filtering)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rRead %i validly-mapped read pairs.\n") % (' ' * 50, new_reads),
            total_reads += new_reads
            del data
        if skip_duplicate_filtering:
            self.stats['total_reads'] = total_reads + self.stats['chr_not_in_fends']
        else:
            self.stats['total_reads'] = total_reads + self.stats['chr_not_in_fends'] + self.stats['pcr_duplicates']
        if self.re:
            self._clean_fend_pairs(fend_pairs)
        total_fend_pairs = 0
        for i in range(len(fend_pairs)):
            for j in range(len(fend_pairs[i])):
                total_fend_pairs += len(fend_pairs[i][j])
        if total_fend_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i total validly-mapped read pairs loaded. %i valid fend pairs\n") %\
                                 (total_reads, total_fend_pairs),
        if self.re and self.binned:
            self._parse_binned_fend_pairs(fend_pairs)
        else:
            self._parse_fend_pairs(fend_pairs)
        self.history += "Success\n"
        return None

    def load_data_from_mat(self, fendfilename, filename):
        """
        Read interaction counts from a :mod:`HiCPipe`-compatible 'mat' text file and place in h5dict.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filename: File name of a 'mat' file containing fend pair and interaction count data.
        :type filename: str.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **fends** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_data_from_mat(fendfilename='%s', filename='%s') - " % (fendfilename, filename)
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.fends = h5py.File(fendfilename, 'r')
        self.maxinsert = None
        if 'binned' in self.fends['/'].attrs and self.fends['/'].attrs['binned'] is not None:
            self.binned = True
        else:
            self.binned = False
        if 'fends' in self.fends and self.fends['fends'] is not None:
            self.re = True
        else:
            self.re = False
            # if fend file indicates no re fragment data, assume this is a binned mat file
            self._load_binned_data_from_mat()
            return None
        self.history = self.fends['/'].attrs['history'] + self.history
        # load data from mat file. This assumes that the mat data was mapped
        # using the same fend numbering as in the fend file.
        chr_indices = self.fends['chr_indices'][...]
        fend_pairs = []
        for i in range(chr_indices.shape[0] - 1):
            fend_pairs.append([])
            for j in range(i + 1): 
                fend_pairs[i].append({})
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("%s not found... no data loaded.\n") % (filename.split('/')[-1]),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        if not self.silent:
            print >> sys.stderr, ("Loading data from mat file..."),
        self.matfile = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                  os.path.dirname(self.file)), os.path.basename(filename))
        input = open(filename, 'r')
        total_reads = 0
        for line in input:
            temp = line.strip('\n').split('\t')
            try:
                fend1 = int(temp[0]) - 1
            except ValueError:
                continue
            fend2 = int(temp[1]) - 1
            count = int(temp[2])
            total_reads += count
            chr1 = numpy.searchsorted(chr_indices, fend1, side='right') - 1
            chr2 = numpy.searchsorted(chr_indices, fend2, side='right') - 1
            if chr2 < chr1:
                fend_pairs[chr1][chr2][(fend2 - chr_indices[chr2], fend1 - chr_indices[chr1])] = count
            elif chr1 < chr2:
                fend_pairs[chr2][chr1][(fend1 - chr_indices[chr1], fend2 - chr_indices[chr2])] = count
            elif fend1 < fend2:
                fend_pairs[chr2][chr1][(fend1 - chr_indices[chr1], fend2 - chr_indices[chr2])] = count
            else:
                fend_pairs[chr1][chr2][(fend2 - chr_indices[chr2], fend1 - chr_indices[chr1])] = count
        input.close()
        self.stats['total_reads'] = total_reads
        self._clean_fend_pairs(fend_pairs)
        total_fend_pairs = 0
        for i in range(len(fend_pairs)):
            for j in range(len(fend_pairs[i])):
                total_fend_pairs += len(fend_pairs[i][j])
        if total_fend_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i valid fend pairs loaded.\n") % (total_fend_pairs),
        # write fend pairs to h5dict
        if self.re and self.binned:
            self._parse_binned_fend_pairs(fend_pairs)
        else:
            self._parse_fend_pairs(fend_pairs)
        self.history += "Success\n"
        return None

    def _load_binned_data_from_mat(self, filename):
        """
        Read binned interaction counts from a 'mat' text file and place in h5dict.

        :param filename: File name of a 'mat' file containing bin pair and interaction count data.
        :type filename: str.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend, the second column contains the idnex of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'cis_data'. For example, all of the downstream cis interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chroosomal fend pairings observed in the data. The first column contains the fend index (from the 'fends' array in the fend object) of the upstream fend (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream fend, and the third column contains the number of reads observed for that fend pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of fends + 1. Each position contains the first entry for the correspondingly-indexed fend in the first column of 'trans_data'. For example, all of the downstream trans interactions for the fend at index 5 in the fend object 'fends' array are in cis_data[cis_indices[5]:cis_indices[6], :].
                     * **fends** (*ndarray*) - A filestream to the hdf5 fend file such that all saved fend attributes can be accessed through this class attribute.

        Data in a binned 'mat' format should include three tab-separaterd columns: bin1, bin2, and count. The bin values should correspond to the bin indices in the fend object.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        # load data from mat file. This assumes that the mat data was mapped
        # using the same bin numbering as in the fend file.
        chr_indices = self.fends['bin_indices'][...]
        bin_pairs = []
        for i in range(chr_indices.shape[0] - 1):
            bin_pairs.append([])
            for j in range(i + 1): 
                bin_pairs[i].append({})
        if not self.silent:
            print >> sys.stderr, ("Loading data from mat file..."),
        self.matfile = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                  os.path.dirname(self.file)), os.path.basename(filename))
        input = open(filename, 'r')
        total_reads = 0
        for line in input:
            temp = line.strip('\n').split('\t')
            try:
                bin1 = int(temp[0])
            except ValueError:
                continue
            bin2 = int(temp[1])
            count = int(temp[2])
            total_reads += count
            chr1 = numpy.searchsorted(chr_indices, bin1, side='right') - 1
            chr2 = numpy.searchsorted(chr_indices, bin2, side='right') - 1
            if chr2 < chr1:
                bin_pairs[chr1][chr2][(bin2 - chr_indices[chr2], bin1 - chr_indices[chr1])] = count
            elif chr1 < chr2:
                bin_pairs[chr2][chr1][(bin1 - chr_indices[chr1], bin2 - chr_indices[chr2])] = count
            elif bin1 < bin2:
                bin_pairs[chr2][chr1][(bin1 - chr_indices[chr1], bin2 - chr_indices[chr2])] = count
            else:
                bin_pairs[chr1][chr2][(bin2 - chr_indices[chr2], bin1 - chr_indices[chr1])] = count
        input.close()
        self.stats['total_reads'] = total_reads
        total_bin_pairs = 0
        for i in range(len(bin_pairs)):
            for j in range(len(bin_pairs[i])):
                total_bin_pairs += len(bin_pairs[i][j])
        if total_bin_pairs == 0:
            if not self.silent:
                print >> sys.stderr, ("No valid data was loaded.\n"),
            self.history += "Error: no valid data loaded\n"
            return None
        if not self.silent:
            print >> sys.stderr, ("%i valid bin pairs loaded.\n") % (total_bin_pairs),
        # write fend pairs to h5dict
        self._parse_fend_pairs(bin_pairs)
        self.history += "Success\n"
        return None

    def load_binned_data_from_matrices(self, fendfilename, filename, format=None):
        """
        Read interaction counts from a tab-separated set of matrix files, one per chromosome, and place in h5dict.

        Each file is assumed to contain a complete matrix of integers divided into equal-width bins. If row and column names are present, the bin ranges will be taken from the labels in the format "XXX|XXX|chrX:XXX-XXX", where only the block of text following the last "|" is looked at. If no labels are present, bins are assumed to begin at coordinate zero.

        :param fendfilename: This specifies the file name of the :class:`Fend` object to associate with the dataset.
        :type fendfilename: str.
        :param filename: The file containing the data matrices. If data are in individual text files, this should be a filename template with an '*' in place of the chromosome(s) names. Each chromosome and chromosome pair in the fend file will be checked and loaded if present. If format is not passed, the format of the matrices will be inferred from the filename (a '*' will default to 'txt', otherwise the filename extension will be used). If data are in hdf5 or npz format, the individual matrices should either be named with the chromosome name or 'N.counts' where N is the chromosome name. For inter-chromosomal interactions, names should be 'N_by_M' for chromosomes N and M.
        :type filename: str.
        :param format: The format of the file(s) to load data from.
        :type format: str.
        :returns: None

        :Attributes: * **fendfilename** (*str.*) - A string containing the relative path of the fend file.
                     * **cis_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero intra-chromosomal bin pairings observed in the data. The first column contains the bin index (from the 'bins' array in the fend object) of the upstream bin, the second column contains the index of the downstream bin, and the third column contains the number of reads observed for that bin pair.
                     * **cis_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of bins + 1. Each position contains the first entry for the correspondingly-indexed bin in the first column of 'cis_data'. For example, all of the downstream cis interactions for the bin at index 5 in the fend object 'bins' array are in cis_data[cis_indices[5]:cis_indices[6], :]. 
                     * **trans_data** (*ndarray*) - A numpy array of type int32 and shape N x 3 where N is the number of valid non-zero inter-chromosomal bin pairings observed in the data. The first column contains the bin index (from the 'bins' array in the fend object) of the upstream bin (upstream also refers to the lower indexed chromosome in this context), the second column contains the index of the downstream bin, and the third column contains the number of reads observed for that bin pair.
                     * **trans_indices** (*ndarray*) - A numpy array of type int64 and a length of the number of bins + 1. Each position contains the first entry for the correspondingly-indexed bin in the first column of 'trans_data'. For example, all of the downstream trans interactions for the bin at index 5 in the fend object 'bins' array are in trans_data[trans_indices[5]:trans_indices[6], :].
                     * **fends** (*ndarray*) - A filestream to the hdf5 fend file such that all saved bin attributes can be accessed through this class attribute.

        When data is loaded the 'history' attribute is updated to include the history of the fend file that becomes associated with it.
        """
        self.history += "HiCData.load_binned_data_from_matrices(fendfilename='%s', filename='%s', format='%s') - " % (fendfilename, filename, str(format))
        # determine if fend file exists and if so, load it
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not found. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        if format is None:
            if filename.count('*'):
                format = 'txt'
            elif filename.split('.')[-1] == 'hdf5':
                format = 'hdf5'
            elif filename.split('.')[-1] == 'npz':
                format = 'npz'
        if format not in ['txt', 'hdf5', 'npz']:
            if not self.silent:
                print >> sys.stderr, ("Could not determine file format. No data loaded.\n"),
            self.history += "Error: File format not recognized.\n"
            return None
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        self.fends = h5py.File(fendfilename, 'r')
        if 'binned' not in self.fends['/'].attrs or self.fends['/'].attrs['binned'] is None:
            if not self.silent:
                print >> sys.stderr, \
                ("The fend file %s was not created for binned data. No data was loaded.\n") % (fendfilename),
            self.history += "Error: '%s' not binned\n" % fendfilename
            return None
        bins = self.fends['bins'][...]
        bin_indices = self.fends['bin_indices'][...]
        self.history = self.fends['/'].attrs['history'] + self.history
        self.chr2int = {}
        chroms = self.fends['chromosomes'][...]
        for i, j in enumerate(chroms):
            self.chr2int[j] = i
        # load matrix files
        data = []
        for i in range(chroms.shape[0]):
            data.append([])
            for j in range(i + 1):
                data[i].append([])
        cis_counts = 0
        trans_counts = 0
        if format == 'txt':
            cis_counts, trans_counts = self._load_txt_matrices(data, filename, chroms, bins, bin_indices)
        elif format == 'hdf5':
            cis_counts, trans_counts = self._load_hdf5_npz_matrices(data, filename, chroms, bins, bin_indices, 'hdf5')
        else:
            cis_counts, trans_counts = self._load_hdf5_npz_matrices(data, filename, chroms, bins, bin_indices, 'npz')
        self.stats['valid_cis_pairs'] = cis_counts
        self.stats['valid_trans_pairs'] = trans_counts
        if cis_counts > 0:
            self.cis_data = numpy.zeros((cis_counts, 3), dtype=numpy.int32)
            pos = 0
            for i in range(chroms.shape[0]):
                if len(data[i][i]) > 0:
                    self.cis_data[pos:(pos + data[i][i].shape[0]), :] = data[i][i]
                    pos += data[i][i].shape[0]
                    data[i][i] = None
            self.stats['valid_cis_reads'] = numpy.sum(self.cis_data[:, 2])
        else:
            self.cis_data = None
        if trans_counts > 0:
            self.trans_data = numpy.zeros((trans_counts, 3), dtype=numpy.int32)
            pos = 0
            for i in range(chroms.shape[0]):
                for j in range(i + 1, chroms.shape[0]):
                    if len(data[j][i]) > 0:
                        self.trans_data[pos:(pos + data[j][i].shape[0]), :] = data[j][i]
                        pos += data[j][i].shape[0]
                        data[j][i] = None
            self.stats['valid_trans_reads'] = numpy.sum(self.trans_data[:, 2])
        else:
            self.trans_data = None
        # create data indices
        if self.cis_data is not None:
            cis_indices = numpy.r_[0, numpy.bincount(self.cis_data[:, 0],
                                   minlength=self.fends['bins'].shape[0])].astype(numpy.int64)
            for i in range(1, cis_indices.shape[0]):
                cis_indices[i] += cis_indices[i - 1]
            self.cis_indices = cis_indices
        else:
            self.cis_indices = None
        if self.trans_data is not None:
            trans_indices = numpy.r_[0, numpy.bincount(self.trans_data[:, 0],
                                     minlength=self.fends['bins'].shape[0])].astype(numpy.int64)
            for i in range(1, trans_indices.shape[0]):
                trans_indices[i] += trans_indices[i - 1]
            self.trans_indices = trans_indices
        else:
            self.trans_indices = None
        # create interaction partner profiles for quality reporting
        cis_reads = 0
        trans_reads = 0
        if self.cis_data is not None:
            fend_profiles = numpy.bincount(self.cis_data[:, 0], minlength=bin_indices[-1])
            fend_profiles += numpy.bincount(self.cis_data[:, 1], minlength=bin_indices[-1])
            self.cis_interaction_distribution = numpy.bincount(fend_profiles)
            cis_reads = numpy.sum(self.cis_data[:, 2])
        if self.trans_data is not None:
            fend_profiles = numpy.bincount(self.trans_data[:, 0], minlength=bin_indices[-1])
            fend_profiles += numpy.bincount(self.trans_data[:, 1], minlength=bin_indices[-1])
            self.trans_interaction_distribution = numpy.bincount(fend_profiles)
            trans_reads = numpy.sum(self.trans_data[:, 2])
        if not self.silent:
            print >> sys.stderr, ("Done  %i cis reads, %i trans reads\n") % (cis_reads, trans_reads),
        return None

    def _load_txt_matrices(self, data, filename, chroms, bins, bin_indices):
        """Load count data from tab-separated matrix files."""
        cis_counts = 0
        trans_counts = 0
        for chrom in chroms:
            for chrom2 in chroms:
                if chrom == chrom2:
                    fname = filename.replace('*', chrom)
                    if not os.path.exists(fname):
                        fname = filename.replace('*', 'chr%s' % chrom)
                        if not os.path.exists(fname):
                            continue
                else:
                    fname = filename.replace('*', '%s_by_%s' % (chrom, chrom2))
                    if not os.path.exists(fname):
                        fname = filename.replace('*', 'chr%s_by_chr%s' % (chrom, chrom2))
                        if not os.path.exists(fname):
                            continue
                col_labels = None
                row_labels = []
                tempdata = []
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rLoading %s...") % (' '*80, fname),
                for line in open(fname):
                    temp = line.rstrip('\n').rstrip('\t').split('\t')
                    try:
                        temp[0] = int(temp[0])
                        temp[0] = str(temp[0])
                    except:
                        if col_labels is None:
                            col_labels = temp[1:]
                            continue
                        else:
                            row_labels.append(temp[0])
                            temp = temp[1:]
                    tempdata.append(numpy.fromstring(' '.join(temp), sep=' ', dtype=numpy.int32))
                if not self.silent:
                    print >> sys.stderr, ("\r%s\r") % (' '*80),
                tempdata = numpy.array(tempdata, dtype=numpy.int32)
                chrint1 = self.chr2int[chrom]
                mapping1 = numpy.zeros(tempdata.shape[0], dtype=numpy.int32) - 1
                if chrom != chrom2:
                    chrint2 = self.chr2int[chrom2]
                    mapping2 = numpy.zeros(tempdata.shape[1], dtype=numpy.int32) - 1
                else:
                    chrint2 = chrint1
                    mapping2 = mapping1
                if col_labels is None:
                    n = bin_indices[chrint1 + 1] - bin_indices[chrint1]
                    mapping1[:n] = numpy.arange(n)
                    valid1 = numpy.where(mapping1 >= 0)[0]
                    if chrom != chrom2:
                        m = bin_indices[chrint2 + 1] - bin_indices[chrint2]
                        mapping2[:m] = numpy.arange(m)
                        valid2 = numpy.where(mapping2 >= 0)[0]
                        tempdata = tempdata[valid1, :][:, valid2]
                    else:
                        tempdata = tempdata[valid1, :][:, valid1]
                else:
                    for i in range(len(row_labels)):
                        temp = row_labels[i].split('|')[-1].split(':')[1].split('-')
                        row_labels[i] = (int(temp[0]) + int(temp[1])) / 2
                    row_labels = numpy.array(row_labels)
                    starts = numpy.searchsorted(bins['start'][bin_indices[chrint1]:bin_indices[chrint1 + 1]],
                                                row_labels, side='right') - 1
                    stops = numpy.searchsorted(bins['stop'][bin_indices[chrint1]:bin_indices[chrint1 + 1]],
                                               row_labels, side='right')
                    valid = numpy.where(starts == stops)[0]
                    mapping1[valid] = starts[valid]
                    invalid1 = numpy.where(mapping1 < 0)[0]
                    if chrom != chrom2:
                        for i in range(len(col_labels)):
                            temp = col_labels[i].split('|')[-1].split(':')[1].split('-')
                            col_labels[i] = (int(temp[0]) + int(temp[1])) / 2
                        col_labels = numpy.array(col_labels)
                        starts = numpy.searchsorted(bins['start'][bin_indices[chrint2]:bin_indices[chrint2 + 1]],
                                                    col_labels, side='right') - 1
                        stops = numpy.searchsorted(bins['stop'][bin_indices[chrint2]:bin_indices[chrint2 + 1]],
                                                   col_labels, side='right')
                        valid = numpy.where(starts == stops)[0]
                        mapping2[valid] = starts[valid]
                        invalid2 = numpy.where(mapping2 < 0)[0]
                        tempdata[invalid1, :] = 0
                        tempdata[:, invalid2] = 0
                    else:
                        tempdata[invalid1, :] = 0
                        tempdata[:, invalid1] = 0
                if chrom != chrom2:
                        where = numpy.where(tempdata > 0)
                        counts = numpy.zeros((where[0].shape[0], 3), dtype=numpy.int32)
                        if chrint1 < chrint2:
                            counts[:, 0] = mapping1[where[0]] + bin_indices[chrint1]
                            counts[:, 1] = mapping2[where[1]] + bin_indices[chrint2]
                        else:
                            counts[:, 1] = mapping1[where[0]] + bin_indices[chrint1]
                            counts[:, 0] = mapping2[where[1]] + bin_indices[chrint2]
                        counts[:, 2] = tempdata[where]
                        trans_counts += counts.shape[0]
                else:
                    indices = numpy.triu_indices(tempdata.shape[0], 0)
                    tempdata = tempdata[indices]
                    where = numpy.where(tempdata > 0)[0]
                    counts = numpy.zeros((where.shape[0], 3), dtype=numpy.int32)
                    counts[:, 0] = mapping1[indices[0][where]] + bin_indices[chrint1]
                    counts[:, 1] = mapping1[indices[1][where]] + bin_indices[chrint1]
                    counts[:, 2] = tempdata[where]
                    cis_counts += counts.shape[0]
                data[chrint2][chrint1] = counts
        return cis_counts, trans_counts

    def _load_hdf5_npz_matrices(self, data, filename, chroms, bins, bin_indices, format):
        """Load count data from hdf5 file."""
        cis_counts = 0
        trans_counts = 0
        if format == 'hdf5':
            infile = h5py.File(filename, 'r')
        else:
            infile = numpy.load(filename)
        for chrom in chroms:
            chrint1 = self.chr2int[chrom]
            for chrom2 in chroms:
                chrint2 = self.chr2int[chrom2]
                if chrom == chrom2:
                    for template in ['*', 'chr*', '*.counts', 'chr*.counts', '*.observed', 'chr*.observed']:
                        name = template.replace('*', chrom)
                        if name in infile:
                            break
                else:
                    for template in ['*_by_+', 'chr*_by_chr+', '*_by_+.counts', 'chr*_by_chr+.counts', '*_by_+.observed', 'chr*_by_chr+.observed']:
                        name = template.replace('*', chrom).replace('+', chrom2)
                        if name in infile:
                            break
                if name not in infile:
                    continue
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rLoading %s...") % (' '*80, name),
                if '%s.positions' % chrom in infile or 'chr%s.positions' % chrom in infile:
                    if '%s.positions' % chrom in infile:
                        pos1 = infile['%s.positions' % chrom][:]
                    else:
                        pos1 = infile['chr%s.positions' % chrom][:]
                    pos1 = (pos1[:, 0] + pos1[:, 1]) / 2
                    if chrom != chrom2:
                        if '%s.positions' % chrom2 in infile:
                            pos2 = infile['%s.positions' % chrom2][:]
                        else:
                            pos2 = infile['chr%s.positions' % chrom2][:]
                        pos2 = (pos2[:, 0] + pos2[:, 1]) / 2
                    else:
                        pos2 = pos1
                else:
                    pos1 = None
                tempdata = infile[name][:]
                if chrom == chrom2 and len(tempdata.shape) == 1:
                    if pos1 is not None:
                        N = pos1.shape[0]
                        if tempdata.shape[0] == N * (N - 1) / 2:
                            diag = False
                        else:
                            diag = True
                    elif format == 'npz' or ('diagonal' in infile['/'].attrs and infile['/'].attrs['diagonal']):
                        diag = True
                        N = (-0.5 + (0.25 + 2 * tempdata.shape[0]))
                    else:
                        diag = False
                        N = (0.5 + (0.25 + 2 * tempdata.shape[0]))
                    temp = numpy.zeros((N, N), dtype=numpy.int32)
                    temp[numpy.triu_indices(N, 1 - int(diag))] = tempdata
                    tempdata = temp
                else:
                    tempdata = infile[name][...]
                    N, M = tempdata.shape
                    if chrom == chrom2:
                        tempdata[numpy.tril_indices(N, 1)] = 0
                if pos1 is None:
                    n = bin_indices[chrint1 + 1] - bin_indices[chrint1]
                    mapping1 = numpy.zeros(N, dtype=numpy.int32) - 1
                    mapping1[:min(n, N)] = numpy.arange(min(n, N))
                    if chrom == chrom2:
                        mapping2 = mapping1
                    else:
                        m = bin_indices[chrint2 + 1] - bin_indices[chrint2]
                        mapping2 = numpy.zeros(M, dtype=numpy.int32) - 1
                        mapping2[:min(m, M)] = numpy.arange(min(m, M))
                else:
                    starts = numpy.searchsorted(bins['start'][bin_indices[chrint1]:bin_indices[chrint1 + 1]],
                                                pos1, side='right') - 1
                    stops = numpy.searchsorted(bins['stop'][bin_indices[chrint1]:bin_indices[chrint1 + 1]],
                                               pos1, side='right')
                    valid1 = numpy.where(starts == stops)[0]
                    mapping1 = numpy.zeros(N, dtype=numpy.int32) - 1
                    mapping1[valid1] = starts[valid1]
                    invalid1 = numpy.where(mapping1 < 0)[0]
                    if chrom == chrom2:
                        mapping2 = mapping1
                        invalid2 = invalid1
                    else:
                        starts = numpy.searchsorted(bins['start'][bin_indices[chrint2]:bin_indices[chrint2 + 1]],
                                                    pos2, side='right') - 1
                        stops = numpy.searchsorted(bins['stop'][bin_indices[chrint2]:bin_indices[chrint2 + 1]],
                                                   pos2, side='right')
                        valid = numpy.where(starts == stops)[0]
                        mapping2 = numpy.zeros(M, dtype=numpy.int32) - 1
                        mapping2[valid] = starts[valid]
                        invalid2 = numpy.where(mapping2 < 0)[0]
                tempdata[invalid1, :] = 0
                tempdata[:, invalid2] = 0
                where = numpy.where(tempdata > 0)
                counts = numpy.zeros((where[0].shape[0], 3), dtype=numpy.int32)
                if chrint1 <= chrint2:
                    counts[:, 0] = mapping1[where[0]] + bin_indices[chrint1]
                    counts[:, 1] = mapping2[where[1]] + bin_indices[chrint2]
                else:
                    counts[:, 1] = mapping1[where[0]] + bin_indices[chrint1]
                    counts[:, 0] = mapping2[where[1]] + bin_indices[chrint2]
                counts[:, 2] = tempdata[where]
                if chrom == chrom2:
                    cis_counts += counts.shape[0]
                else:
                    trans_counts += counts.shape[0]
                data[chrint2][chrint1] = counts
        infile.close()
        if not self.silent:
            print >> sys.stderr, ("\r%s\r") % (' '*80),
        return cis_counts, trans_counts

    def _find_fend_pairs(self, data, fend_pairs, skip_duplicate_filtering=False):
        """Return array with lower fend, upper fend, and count for pair."""
        chroms = self.fends['chromosomes'][...]
        if 'cuts' not in self.__dict__.keys():
            self._find_cut_sites()
        # for looking at fragment end distribution across genome in fragments with >insert size
        """
        if 'invalid_distribution' not in self.__dict__.keys():
            self.invalid_starts = numpy.zeros(chroms.shape[0], dtype=numpy.int32)
            self.invalid_distribution = []
            self.invalid_binsize = 100
            for i, chrom in enumerate(chroms):
                start = (self.cuts[i][0] / self.invalid_binsize) * self.invalid_binsize
                stop = ((self.cuts[i][-1] - 1) / self.invalid_binsize + 1) * self.invalid_binsize
                self.invalid_starts[i] = start
                self.invalid_distribution.append(numpy.zeros(((stop - start) / self.invalid_binsize, 2),
                                                              dtype=numpy.int32))
        """
        # assign fends on a per-chromosome pair basis
        if 'non_cis_invalid_insert' not in self.__dict__.keys():
            self.non_cis_invalid_insert = 0
            self.different_fragment_invalid_insert = 0
        for i in range(len(data)):
            for j in range(len(data[i])):
                if data[i][j].shape[0] == 0:
                    continue
                mapped_fends = numpy.empty((data[i][j].shape[0], 2), dtype=numpy.int32)
                distances = numpy.zeros(data[i][j].shape[0], dtype=numpy.int32)
                signs = numpy.minimum(0, numpy.sign(data[i][j]))
                data[i][j] = numpy.abs(data[i][j])
                if not self.silent:
                    print >> sys.stderr, ("\rMapping fends for %s by %s") % (chroms[i].ljust(10), chroms[j].ljust(10)),
                mapped_fends[:, 0] = numpy.searchsorted(self.cuts[j], data[i][j][:, 0])
                mapped_fends[:, 1] = numpy.searchsorted(self.cuts[i], data[i][j][:, 1])
                # make sure coordinates are within first and last cutsites
                valid = numpy.where((mapped_fends[:, 0] > 0) * (mapped_fends[:, 0] < self.cuts[j].shape[0]) *
                                    (mapped_fends[:, 1] > 0) * (mapped_fends[:, 1] < self.cuts[i].shape[0]))[0]
                if skip_duplicate_filtering:
                    self.stats['out_of_bounds'] += numpy.sum(data[i][j][:, 2]) - numpy.sum(data[i][j][valid, 2])
                else:
                    self.stats['out_of_bounds'] += mapped_fends.shape[0] - valid.shape[0]
                # find distance from cutsite to mapped coordinate
                distances[valid] = (numpy.abs(data[i][j][valid, 0] - self.cuts[j][mapped_fends[valid, 0] +
                                    signs[valid, 0]]) + numpy.abs(data[i][j][valid, 1] -
                                    self.cuts[i][mapped_fends[valid, 1] + signs[valid, 1]]))
                # add distances to distance distribution counts
                self.insert_distribution[:, 0] += numpy.bincount(numpy.searchsorted(self.insert_distribution[1:, 1],
                    distances[valid], side='right'), minlength=self.insert_distribution.shape[0])
                # remove fends with too great an insert distance
                invalid = valid[numpy.where(distances[valid] > self.maxinsert)[0]]
                if i != j:
                    self.non_cis_invalid_insert += invalid.shape[0]
                # for looking at fragment end distribution across genome in fragments with >insert size
                """
                strand = numpy.where(signs[invalid, 1] == 0)[0]
                self.invalid_distribution[i][:, 0] += numpy.bincount((data[i][j][invalid[strand], 1] -
                    self.invalid_starts[i]) / self.invalid_binsize, minlength=self.invalid_distribution[i].shape[0])
                strand = numpy.where(signs[invalid, 0] == 0)[0]
                self.invalid_distribution[j][:, 0] += numpy.bincount((data[i][j][invalid[strand], 0] -
                    self.invalid_starts[j]) / self.invalid_binsize, minlength=self.invalid_distribution[j].shape[0])
                strand = numpy.where(signs[invalid, 1] == -1)[0]
                self.invalid_distribution[i][:, 1] += numpy.bincount((data[i][j][invalid[strand], 1] -
                    self.invalid_starts[i]) / self.invalid_binsize, minlength=self.invalid_distribution[i].shape[0])
                strand = numpy.where(signs[invalid, 0] == -1)[0]
                self.invalid_distribution[j][:, 1] += numpy.bincount((data[i][j][invalid[strand], 0] -
                    self.invalid_starts[j]) / self.invalid_binsize, minlength=self.invalid_distribution[j].shape[0])
                del invalid
                del strand
                """
                valid1 = numpy.where(distances[valid] <= self.maxinsert)[0]
                if skip_duplicate_filtering:
                    self.stats['insert_size'] += (numpy.sum(data[i][j][valid, 2]) -
                                                  numpy.sum(data[i][j][valid[valid1], 2]))
                else:
                    self.stats['insert_size'] += valid.shape[0] - valid1.shape[0]
                valid = valid[valid1]
                # convert from fragments to fends
                mapped_fends[valid, :2] = mapped_fends[valid, :2] * 2 - 1 + signs[valid, :2]
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rCounting fend pairs...") % (' ' * 50),
                for k in valid:
                    key = (mapped_fends[k, 0], mapped_fends[k, 1])
                    if key not in fend_pairs[i][j]:
                        if skip_duplicate_filtering:
                            fend_pairs[i][j][key] = data[i][j][k, 2]
                        else:
                            fend_pairs[i][j][key] = 1
                    else:
                        if skip_duplicate_filtering:
                            fend_pairs[i][j][key] += data[i][j][k, 2]
                        else:
                            fend_pairs[i][j][key] += 1
        if not self.silent and not skip_duplicate_filtering:
            print >> sys.stderr, ("Done\n"),
        return None

    def _find_bin_pairs(self, data, bin_pairs, skip_duplicate_filtering=False):
        """Return array with lower bin, upper bin, and count for pair."""
        chroms = self.fends['chromosomes'][...]
        bins = self.fends['bins'][...]
        bin_indices = self.fends['bin_indices'][...]
        if 'maxinsert' in self.__dict__ and self.maxinsert is not None:
            maxinsert = self.maxinsert
        else:
            maxinsert = None
        # assign fends on a per-chromosome pair basis
        for i in range(len(data)):
            for j in range(len(data[i])):
                if data[i][j].shape[0] == 0:
                    continue
                mapped_bins = numpy.zeros((data[i][j].shape[0], 2), dtype=numpy.int32) - 1
                signs = numpy.minimum(0, numpy.sign(data[i][j]))
                data[i][j] = numpy.abs(data[i][j])
                if not self.silent:
                    print >> sys.stderr, ("\rMapping bins for %s by %s") % (chroms[i].ljust(10), chroms[j].ljust(10)),
                starts = numpy.searchsorted(bins['start'][bin_indices[j]:bin_indices[j + 1]],
                                            data[i][j][:, 0], side='right') - 1
                stops = numpy.searchsorted(bins['stop'][bin_indices[j]:bin_indices[j + 1]], data[i][j][:, 0])
                valid = numpy.where(starts == stops)[0]
                mapped_bins[valid, 0] = starts[valid]
                starts = numpy.searchsorted(bins['start'][bin_indices[i]:bin_indices[i + 1]],
                                            data[i][j][:, 1], side='right') - 1
                stops = numpy.searchsorted(bins['stop'][bin_indices[i]:bin_indices[i + 1]], data[i][j][:, 1])
                valid = numpy.where(starts == stops)[0]
                mapped_bins[valid, 1] = starts[valid]
                valid = numpy.where((mapped_bins[:, 0] >= 0) * (mapped_bins[:, 1] >= 0))[0]
                if skip_duplicate_filtering:
                    self.stats['out_of_bounds'] += numpy.sum(data[i][j][:, 2]) - numpy.sum(data[i][j][valid, 2])
                else:
                    self.stats['out_of_bounds'] += mapped_bins.shape[0] - valid.shape[0]
                # if filtering by insert size, determine distances
                if i == j and maxinsert is not None:
                    distances = numpy.zeros(data[i][j].shape[0], dtype=numpy.int32)
                    distances[valid] = data[i][j][valid, 1] - data[i][j][valid, 0]
                    valid1 = numpy.where((distances[valid] > maxinsert) + (signs[valid, 0] == signs[valid, 1]))[0]
                    if skip_duplicate_filtering:
                        self.stats['insert_size'] += (numpy.sum(data[i][j][valid, 2]) -
                                                      numpy.sum(data[i][j][valid[valid1], 2]))
                    else:
                        self.stats['insert_size'] += valid.shape[0] - valid1.shape[0]
                    valid = valid[valid1]
                # convert to bin paired
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rCounting bin pairs...") % (' ' * 50),
                for k in valid:
                    key = (mapped_bins[k, 0], mapped_bins[k, 1])
                    if key not in bin_pairs[i][j]:
                        if skip_duplicate_filtering:
                            bin_pairs[i][j][key] = data[i][j][k, 2]
                        else:
                            bin_pairs[i][j][key] = 1
                    else:
                        if skip_duplicate_filtering:
                            bin_pairs[i][j][key] += data[i][j][k, 2]
                        else:
                            bin_pairs[i][j][key] += 1
        if not self.silent and not skip_duplicate_filtering:
            print >> sys.stderr, ("Done\n"),
        return None

    def _clean_fend_pairs(self, fend_pairs):
        chr_indices = self.fends['chr_indices'][...]
        # remove fend pairs from same fend or opposite strand adjacents
        for i in range(len(fend_pairs)):
            for j in range(0, chr_indices[i + 1] - chr_indices[i] - 2, 2):
                # same fend
                name = (j, j)
                if name in fend_pairs[i][i]:
                    self.stats['same_fragment'] += fend_pairs[i][i][name]
                    del fend_pairs[i][i][name]
                name = (j + 1, j + 1)
                if name in fend_pairs[i][i]:
                    self.stats['same_fragment'] += fend_pairs[i][i][name]
                    del fend_pairs[i][i][name]
                # same fend
                name = (j, j + 1)
                if name in fend_pairs[i][i]:
                    self.stats['same_fragment'] += fend_pairs[i][i][name]
                    del fend_pairs[i][i][name]
                # same fend
                name = (j + 1, j)
                if name in fend_pairs[i][i]:
                    self.stats['same_fragment'] += fend_pairs[i][i][name]
                    del fend_pairs[i][i][name]
                # adjacent fends, opposite strands
                name = (j, j + 3)
                if name in fend_pairs[i][i]:
                    self.stats['failed_cut'] += fend_pairs[i][i][name]
                    del fend_pairs[i][i][name]
                name = (j + 1, j + 2)
                if name in fend_pairs[i][i]:
                    self.stats['failed_cut'] += fend_pairs[i][i][name]
                    del fend_pairs[i][i][name]
            j = chr_indices[i + 1] - chr_indices[i] - 2
            # same fend
            name = (j, j)
            if name in fend_pairs[i][i]:
                self.stats['same_fragment'] += fend_pairs[i][i][name]
                del fend_pairs[i][i][name]
            name = (j + 1, j + 1)
            if name in fend_pairs[i][i]:
                self.stats['same_fragment'] += fend_pairs[i][i][name]
                del fend_pairs[i][i][name]
            name = (j, j + 1)
            if name in fend_pairs[i][i]:
                self.stats['same_fragment'] += fend_pairs[i][i][name]
                del fend_pairs[i][i][name]
            name = (j + 1, j)
            if name in fend_pairs[i][i]:
                self.stats['same_fragment'] += fend_pairs[i][i][name]
                del fend_pairs[i][i][name]
        return None

    def _parse_fend_pairs(self, fend_pairs):
        """Separate fend pairs into cis and trans interactions index."""
        if not self.silent:
            print >> sys.stderr, ("Parsing fend pairs..."),
        if self.binned and not self.re:
            chr_indices = self.fends['bin_indices'][...]
        else:
            chr_indices = self.fends['chr_indices'][...]
        # determine number of cis pairs
        cis_count = 0
        for i in range(len(fend_pairs)):
            cis_count += len(fend_pairs[i][i])
        self.stats['valid_cis_pairs'] = cis_count
        # create cis array
        self.cis_data = numpy.empty((cis_count, 3), dtype=numpy.int32)
        pos = 0
        # fill in each chromosome's cis interactions
        for i in range(len(fend_pairs)):
            if len(fend_pairs[i][i]) == 0:
                continue
            chr_start = pos
            for key, count in fend_pairs[i][i].iteritems():
                self.cis_data[pos, :2] = key
                self.cis_data[pos, 2] = count
                pos += 1
            fend_pairs[i][i] = None
            # sort interactions
            order = numpy.lexsort((self.cis_data[chr_start:pos, 1], self.cis_data[chr_start:pos, 0]))
            self.cis_data[chr_start:pos, :] = self.cis_data[order + chr_start, :]
            self.cis_data[chr_start:pos, :2] += chr_indices[i]
            del order
        self.stats['valid_cis_reads'] += numpy.sum(self.cis_data[:, 2])
        # determine number of trans pairs
        trans_count = 0
        for i in range(len(fend_pairs)):
            for j in range(i + 1, len(fend_pairs)):
                trans_count += len(fend_pairs[j][i])
        self.stats['valid_trans_pairs'] = trans_count
        # create trans array
        self.trans_data = numpy.empty((trans_count, 3), dtype=numpy.int32)
        pos = 0
        # fill in each chromosome's trans interactions
        for i in range(len(fend_pairs) - 1):
            chr1_start = pos
            for j in range(i + 1, len(fend_pairs)):
                if len(fend_pairs[j][i]) == 0:
                    continue
                chr2_start = pos
                for key, count in fend_pairs[j][i].iteritems():
                    self.trans_data[pos, :2] = key
                    self.trans_data[pos, 2] = count
                    pos += 1
                fend_pairs[j][i] = None
                self.trans_data[chr2_start:pos, 0] += chr_indices[i]
                self.trans_data[chr2_start:pos, 1] += chr_indices[j]
            # sort interactions
            order = numpy.lexsort((self.trans_data[chr1_start:pos, 1], self.trans_data[chr1_start:pos, 0]))
            self.trans_data[chr1_start:pos, :] = self.trans_data[order + chr1_start, :]
            del order
        self.stats['valid_trans_reads'] += numpy.sum(self.trans_data[:, 2])
        # create data indices
        if self.cis_data.shape[0] > 0:
            cis_indices = numpy.r_[0, numpy.bincount(self.cis_data[:, 0],
                                   minlength=chr_indices[-1])].astype(numpy.int64)
            for i in range(1, cis_indices.shape[0]):
                cis_indices[i] += cis_indices[i - 1]
            self.cis_indices = cis_indices
        else:
            self.cis_data = None
        if self.trans_data.shape[0] > 0:
            trans_indices = numpy.r_[0, numpy.bincount(self.trans_data[:, 0],
                                     minlength=chr_indices[-1])].astype(numpy.int64)
            for i in range(1, trans_indices.shape[0]):
                trans_indices[i] += trans_indices[i - 1]
            self.trans_indices = trans_indices
        else:
            self.trans_data = None
        # create interaction partner profiles for quality reporting
        cis_reads = 0
        trans_reads = 0
        if self.cis_data is not None:
            fend_profiles = numpy.bincount(self.cis_data[:, 0], minlength=chr_indices[-1])
            fend_profiles += numpy.bincount(self.cis_data[:, 1], minlength=chr_indices[-1])
            self.cis_interaction_distribution = numpy.bincount(fend_profiles)
            cis_reads = numpy.sum(self.cis_data[:, 2])
        if self.trans_data is not None:
            fend_profiles = numpy.bincount(self.trans_data[:, 0], minlength=chr_indices[-1])
            fend_profiles += numpy.bincount(self.trans_data[:, 1], minlength=chr_indices[-1])
            self.trans_interaction_distribution = numpy.bincount(fend_profiles)
            trans_reads = numpy.sum(self.trans_data[:, 2])
        if not self.silent:
            print >> sys.stderr, ("Done  %i cis reads, %i trans reads\n") % (cis_reads, trans_reads),
        return None

    def _parse_binned_fend_pairs(self, fend_pairs):
        """Separate fend pairs into cis and trans interactions index."""
        if not self.silent:
            print >> sys.stderr, ("Parsing fend pairs..."),
        chr_indices = self.fends['chr_indices'][...]
        # determine number of cis pairs
        cis_count = 0
        for i in range(len(fend_pairs)):
            cis_count += len(fend_pairs[i][i])
        self.stats['valid_cis_pairs'] = cis_count
        # determine number of bin pairs present
        cis_count = 0
        bin_indices = self.fends['bin_indices'][...]
        bins = self.fends['bins'][...]
        binsize = self.fends['/'].attrs['binned']
        indices = []
        for i in range(len(fend_pairs)):
            mapping = (self.fends['fends']['mid'][chr_indices[i]:chr_indices[i + 1]] -
                       bins['start'][bin_indices[i]]) / binsize
            n = bin_indices[i + 1] - bin_indices[i]
            data = numpy.zeros(n * (n + 1) / 2, dtype=numpy.int32)
            for key, count in fend_pairs[i][i].iteritems():
                bin1 = mapping[key[0]]
                bin2 = mapping[key[1]]
                index = bin1 * (n - 1) - (bin1 * (bin1 - 1) / 2) + bin2
                data[index] += count
            indices.append(numpy.where(data > 0)[0])
            cis_count += indices[-1].shape[0]
            fend_pairs[i][i] = data[indices[-1]]
        # create cis array
        self.cis_data = numpy.empty((cis_count, 3), dtype=numpy.int32)
        pos = 0
        # fill in each chromosome's cis interactions
        for i in range(len(fend_pairs)):
            if indices[i].shape[0] == 0:
                continue
            n = bin_indices[i + 1] - bin_indices[i]
            triu_indices = numpy.triu_indices(n, 0)
            self.cis_data[pos:(pos + indices[i].shape[0]), 0] = triu_indices[0][indices[i]] + bin_indices[i]
            self.cis_data[pos:(pos + indices[i].shape[0]), 1] = triu_indices[1][indices[i]] + bin_indices[i]
            self.cis_data[pos:(pos + indices[i].shape[0]), 2] = fend_pairs[i][i]
            pos += indices[i].shape[0]
            indices[i] = None
            fend_pairs[i][i] = None
        self.stats['valid_cis_reads'] += numpy.sum(self.cis_data[:, 2])
        # determine number of trans pairs
        trans_count = 0
        for i in range(len(fend_pairs)):
            for j in range(i + 1, len(fend_pairs)):
                trans_count += len(fend_pairs[j][i])
        self.stats['valid_trans_pairs'] = trans_count        
        # determine number of bin pairs present
        trans_count = 0
        indices = []
        for i in range(len(fend_pairs)):
            indices.append([])
        for i in range(len(fend_pairs)):
            mapping1 = (self.fends['fends']['mid'][chr_indices[i]:chr_indices[i + 1]] -
                        bins['start'][bin_indices[i]]) / binsize
            n = bin_indices[i + 1] - bin_indices[i]
            for j in range(i + 1, len(fend_pairs)):
                mapping2 = (self.fends['fends']['mid'][chr_indices[j]:chr_indices[j + 1]] -
                            bins['start'][bin_indices[j]]) / binsize
                m = bin_indices[j + 1] - bin_indices[j]
                data = numpy.zeros((n, m), dtype=numpy.int32)
                for key, count in fend_pairs[j][i].iteritems():
                    bin1 = mapping1[key[0]]
                    bin2 = mapping2[key[1]]
                    data[bin1, bin2] += count
                indices[j].append(numpy.where(data > 0))
                trans_count += indices[j][-1][0].shape[0]
                fend_pairs[j][i] = data[indices[j][i]]
        # create trans array
        self.trans_data = numpy.empty((trans_count, 3), dtype=numpy.int32)
        pos = 0
        # fill in each chromosome's cis interactions
        for i in range(len(fend_pairs)):
            for j in range(i + 1, len(fend_pairs)):
                if indices[j][i][0].shape[0] == 0:
                    continue
                self.trans_data[pos:(pos + indices[j][i][0].shape[0]), 0] = indices[j][i][0] + bin_indices[i]
                self.trans_data[pos:(pos + indices[j][i][0].shape[0]), 1] = indices[j][i][1] + bin_indices[j]
                self.trans_data[pos:(pos + indices[j][i][0].shape[0]), 2] = fend_pairs[j][i]
                pos += indices[j][i][0].shape[0]
                indices[j][i] = None
                fend_pairs[j][i] = None
        self.stats['valid_trans_reads'] += numpy.sum(self.trans_data[:, 2])
        # create data indices
        if self.cis_data.shape[0] > 0:
            cis_indices = numpy.r_[0, numpy.bincount(self.cis_data[:, 0],
                                   minlength=self.fends['fends'].shape[0])].astype(numpy.int64)
            for i in range(1, cis_indices.shape[0]):
                cis_indices[i] += cis_indices[i - 1]
            self.cis_indices = cis_indices
        else:
            self.cis_data = None
        if self.trans_data.shape[0] > 0:
            trans_indices = numpy.r_[0, numpy.bincount(self.trans_data[:, 0],
                                     minlength=self.fends['fends'].shape[0])].astype(numpy.int64)
            for i in range(1, trans_indices.shape[0]):
                trans_indices[i] += trans_indices[i - 1]
            self.trans_indices = trans_indices
        else:
            self.trans_data = None
        # create interaction partner profiles for quality reporting
        cis_reads = 0
        trans_reads = 0
        if self.cis_data is not None:
            fend_profiles = numpy.bincount(self.cis_data[:, 0], minlength=chr_indices[-1])
            fend_profiles += numpy.bincount(self.cis_data[:, 1], minlength=chr_indices[-1])
            self.cis_interaction_distribution = numpy.bincount(fend_profiles)
            cis_reads = numpy.sum(self.cis_data[:, 2])
        if self.trans_data is not None:
            fend_profiles = numpy.bincount(self.trans_data[:, 0], minlength=chr_indices[-1])
            fend_profiles += numpy.bincount(self.trans_data[:, 1], minlength=chr_indices[-1])
            self.trans_interaction_distribution = numpy.bincount(fend_profiles)
            trans_reads = numpy.sum(self.trans_data[:, 2])
        if not self.silent:
            print >> sys.stderr, ("Done  %i cis reads, %i trans reads\n") % (cis_reads, trans_reads),
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
