#!/usr/bin/env python

"""A class for handling HiC fend information."""

import os
import sys

import numpy
import h5py
#from reportlab.graphics.shapes import Drawing
#from reportlab.graphics.charts.lineplots import LinePlot
#from reportlab.lib.colors import Color


class Fend(object):

    """
    This class handles restriction enzyme digest-generated fragment data for HiC experiments.

    This class stores a list of chromosomes, a dictionary for converting from chromosome label to integer and back, fragment starts, stops, and chromosome number in an h5dict.
        
    .. note::
      This class is also available as hifive.Fend

    When initialized, this class creates an h5dict in which to store all data associated with this object.

    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`Fend <hifive.fend.Fend>` class object

    :attributes: * **file** (*str.*) - A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcome.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a Fend object."""
        self.file = filename
        self.silent = silent
        self.history = ''
        self.filetype = 'fend'
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
        Save fend data to h5dict.

        :returns: None
        """
        self.history.replace("'None'", "None")
        fendfile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['file', 'silent']:
                continue
            elif isinstance(self[key], numpy.ndarray):
                fendfile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                fendfile.attrs[key] = self[key]
        fendfile.close()
        return None

    def load(self):
        """
        Load fend data from h5dict specified at object creation.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :returns: None
        """
        fendfile = h5py.File(self.file, 'r')
        for key in fendfile.keys():
            self[key] = numpy.copy(fendfile[key])
        for key in fendfile['/'].attrs.keys():
            self[key] = fendfile['/'].attrs[key]
        fendfile.close()
        return None

    def load_fends(self, filename, genome_name=None, re_name=None, format=None):
        """
        Parse and store fend data in h5dict.

        :param filename: A file name to read restriction fragment data from. The file may be a 'mat' file compatible with HiCPipe, or a BED file containing RE fragment boundaries or cutsites.
        :type filename: str.
        :param genome_name: The name of the species and build. Optional.
        :type genome_name: str.
        :param re_name: The name of the restriction enzyme used to produce the fragment set. Optional.
        :type re_name: str.
        :param format: Format of the input file. If not specified, it will be inferred from the file extension. Optional.
        :type format: str.
        :returns: None

        :Attributes: * **chromosomes** (*ndarray*) - A numpy array containing chromosome names as strings. The position of the chromosome name in this array is referred to as the chromosome index.
             * **fends** (*ndarray*) - A numpy array of length N where N is the number of fends and containing the fields 'chr', 'start', 'stop', and 'mid'. All of these are of type int32. The 'chr' field contains the index of the chromosome. If the bed file or fend file used to create the Fend object contains additional columns, these features are also included as fields with names corresponding to the header names. These additional fields are of type float32. If produced from a bed file, fends are sorted by chromosome (the order in the 'chromosomes' array) and then by coordinates.
             * **chr_indices** (*ndarray*) - A numpy array with a length of the number of chromosomes in 'chromosomes' + 1. This array contains the first position in 'fragments' for the chromosome in the corresponding position in the 'chromosomes' array. The last position in the array contains the total number of fragments.
             * **genome_name** (*str.*) - A string (or None if not passed as argument) of the genome from which the fends originated.
             * **re_name** (*str.*) - A string (or None if not passed as argument) of the restriction enzyme used to produce the fends.
        """
        self.history += "Fend.load_fends(filename='%s', genome_name='%s', re_name='%s', format='%s') - " % (filename, genome_name, re_name, format)
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename),
            self.history += "Error: '%s' no located\n" % filename
            return None
        # if no genome name given, determine from filename
        if genome_name is None:
            genome_name = filename.split('/')[-1].split('_')[0]
        # if no re name given, determine from filename
        if re_name is None:
            re_name = filename.split('_')[-1].split('.')[0]
        if format == 'bed' or (format is None and filename.split('.')[-1] == 'bed'):
            fends, chromosomes = self._load_from_bed(filename)
        elif format == "fend" or (format is None and filename.split('.')[-1] == 'fend'):
            fends, chromosomes = self._load_from_fend(filename)
        else:
            if not self.silent:
                print >> sys.stderr, ("Unrecognized format.")
            self.history += "Error: Unrecognized fend format\n"
            return None
        # make note of chromosome positions in fend array
        chr_indices = numpy.zeros(chromosomes.shape[0] + 1, dtype=numpy.int32)
        chr_indices[1:] = numpy.bincount(fends['chr'])
        for i in range(1, chr_indices.shape[0]):
            chr_indices[i] += chr_indices[i - 1]
        # write all data to h5dict
        self.fends = fends
        self.re_name = re_name
        self.genome_name = genome_name
        self.chr_indices = chr_indices
        self.chromosomes = chromosomes
        # calculate statistics
        #self.fend_sizes = numpy.zeros((40, self.chromosomes.shape[0] + 1), dtype=numpy.float32)
        #sizes = numpy.log(self.fends['stop'] - self.fends['start'])
        #splits = numpy.linspace(numpy.amin(sizes), numpy.amax(sizes) + 1, 41)
        #for i in range(self.chr_indices.shape[0] - 1):
        #    self.fend_sizes[:, i] = numpy.histogram(sizes[chr_indices[i]:chr_indices[i + 1]], bins=splits)[0]
        #self.fend_sizes[:, -1] = (splits[1:] + splits[:-1]) / 2.0
        self.history += "Success\n"
        return None

    def _load_from_fend(self, fname):
        chromosomes = []
        chr2int = {}
        data = {}
        input = open(fname, 'r')
        chromosome_index = 2
        coordinate_index = 3
        length_index = 5
        feature_names = []
        for line in input:
            temp = line.strip('\n').split('\t')
            try:
                temp[0] == int(temp[0])
            except:
                fend_index = temp.index('fend')
                chromosome_index = temp.index('chr')
                coordinate_index = temp.index('coord')
                for i in range(6, len(temp)):
                    feature_names.append(temp[i])
                continue
            chrom = temp[chromosome_index].strip('chr')
            fend = int(temp[fend_index]) - 1
            if chrom not in chr2int:
                chr2int[chrom] = len(chromosomes)
                chromosomes.append(chrom)
            features = []
            for i in range(6, 6 + len(feature_names)):
                features.append(float(temp[i]))
            length = int(temp[length_index])
            if fend % 2 == 0:
                start = int(temp[coordinate_index])
                stop = start + length / 2
            else:
                stop = int(temp[coordinate_index]) + 1
                start = stop - length + length / 2
            data[fend] = [chr2int[chrom], start, stop] + features
        input.close()
        dtypes = [('chr', numpy.int32), ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32)]
        for i in range(len(feature_names)):
            dtypes.append((feature_names[i], numpy.float32))
        fends = numpy.empty(len(data), dtype=numpy.dtype(dtypes))
        for i in data.keys():
            fends['chr'][i] = data[i][0]
            fends['start'][i] = data[i][1]
            fends['stop'][i] = data[i][2]
            fends['mid'][i] = (data[i][1] + data[i][2]) / 2
            for j in range(3, 3 + len(feature_names)):
                fends[feature_names[j - 3]][i] = data[i][j]
        chromosomes = numpy.array(chromosomes)
        return [fends, chromosomes]

    def _load_from_bed(self, fname):
        chromosomes = []
        data = {}
        feature_names = []
        sizes = []
        input = open(fname, 'r')
        for line in input:
            temp = line.strip('\n').split('\t')
            try:
                start = int(temp[1])
            except:
                for i in range(6, len(temp)):
                    feature_names.append(temp[i])
                continue
            chrom = temp[0].strip('chr')
            if chrom not in data:
                data[chrom] = []
            stop = int(temp[2])
            features = []
            for i in range(6, 6 + len(feature_names)):
                temp1, temp2 = temp[i].split(',')
                features += [float(temp1), float(temp2)]
            data[chrom].append(tuple([start, stop] + features))
            sizes.append(stop - start)
        input.close()
        chromosomes = data.keys()
        chromosomes.sort()
        chromosomes = numpy.array(chromosomes)
        sizes = numpy.array(sizes)
        dtypes = [('chr', numpy.int32), ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32)]
        dtypes2 = [('start', numpy.int32), ('stop', numpy.int32)]
        for i in range(len(feature_names)):
            dtypes.append((feature_names[i], numpy.float32))
            dtypes2.append(("%s_1" % feature_names[i], numpy.float32))
            dtypes2.append(("%s_2" % feature_names[i], numpy.float32))
        for i, chrom in enumerate(chromosomes):
            data[chrom].sort()
            data[chrom] = numpy.array(data[chrom], dtype=numpy.dtype(dtypes2))
        if numpy.mean(sizes) < 10.0:
            fends = numpy.empty((sizes.shape[0] - chromosomes.shape[0]) * 2, dtype=numpy.dtype(dtypes))
            pos = 0
            for i, chrom in enumerate(chromosomes):
                data_len = (data[chrom].shape[0] - 1) * 2
                cuts = (data[chrom]['start'][:] + data[chrom]['stop'][:]) / 2
                mids = (cuts[:-1] + cuts[1:]) / 2
                fends['chr'][pos:(pos + data_len)] = i
                fends['start'][pos:(pos + data_len):2] = cuts[:-1]
                fends['stop'][pos:(pos + data_len):2] = mids
                fends['start'][(pos + 1):(pos + data_len):2] = mids
                fends['stop'][(pos + 1):(pos + data_len):2] = cuts[1:]
                for j in range(len(feature_names)):
                    fends[feature_names[j]][pos:(pos + data_len):2] = data[chrom]["%s_2" % feature_names[j]][:-1]
                    fends[feature_names[j]][(pos + 1):(pos + data_len):2] = data[chrom]["%s_1" % feature_names[j]][1:]
                pos += data_len
        else:
            fends = numpy.empty(sizes.shape[0] * 2, dtype=numpy.dtype(dtypes))
            pos = 0
            for i, chrom in enumerate(chromosomes):
                data_len = data[chrom].shape[0] * 2
                mids = (data[chrom]['start'][:] + data[chrom]['stop'][:]) / 2
                fends['chr'][pos:(pos + data_len)] = i
                fends['start'][pos:(pos + data_len):2] = data[chrom]['start'][:]
                fends['start'][(pos + 1):(pos + data_len):2] = mids
                fends['stop'][pos:(pos + data_len):2] = mids
                fends['stop'][(pos + 1):(pos + data_len):2] = data[chrom]['stop'][:]
                for j in range(len(feature_names)):
                    fends[feature_names[j]][pos:(pos + data_len):2] = data[chrom]["%s_1" % feature_names[j]]
                    fends[feature_names[j]][(pos + 1):(pos + data_len):2] = data[chrom]["%s_2" % feature_names[j]]
                pos += data_len
        fends['mid'][:] = (fends['start'][:] + fends['stop'][:]) / 2
        return [fends, chromosomes]

    def plot_statistics(self):
        width, height = 800, 400
        c = Drawing(width + 5, height + 5)
        data = []
        where = numpy.where(numpy.sum(self.fend_sizes[:, :-1], axis=0))[0]
        for i in where:
            data.append(tuple(zip(self.fend_sizes[:, -1],
                        self.fend_sizes[:, i] / float(numpy.sum(self.fend_sizes[:, i])))))
        data.append(tuple(zip(self.fend_sizes[:, -1], numpy.sum(self.fend_sizes[:, :-1], axis=1) /
                    float(numpy.sum(self.fend_sizes[:, :-1])))))
        lp = LinePlot()
        lp.x = 40
        lp.y = 20
        lp.height = height - 20
        lp.width = width - 40
        lp.data = data
        lp.joinedLines = 1
        def findcolor(i, n):
            H = i * 6.0 / float(n)
            X = 1.0 - abs(H % 2 - 1.0)
            if H < 1.0:
                return (1.0, X, 0., 1.0)
            elif H < 2.0:
                return (X, 1.0, 0., 1.0)
            elif H < 3.0:
                return (0., 1.0, X, 1.0)
            elif H < 4.0:
                return (0., X, 1.0, 1.0)
            elif H < 5.0:
                return (X, 0., 1.0, 1.0)
            else:
                return (1.0, 0., X, 1.0)

        def axisstring(i):
            return "%0.1E" % numpy.exp(i)

        lp.xValueAxis.labelTextFormat = axisstring
        for i in range(where.shape[0]):
            lp.lines[i].strokeColor = Color(*findcolor(i, where.shape[0]))
        c.add(lp)
        #c.drawString(220, 0, "Length of fend (bp)")
        #c.drawString(0, 110, "Percent of fends")
        return c


if __name__ == '__main__':
    filename = sys.argv[1]
    fend = Fend('.'.join(filename.split('.')[:-1] + ['.hdf5']), 'w')
    fend.load_fends(filename)
    fend.fends.close()
