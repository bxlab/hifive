#!/usr/bin/env python

"""A class for handling 5C fragment information."""

import os
import sys

import numpy
import h5py


class Fragment(object):

    """
    This class handles restriction enzyme digest-generated fragment data for 5C experiments.

    The Fragment class contains all of the genomic information for a 5C experiment, including fragment locations and orientations, chromosome name mapping, and region locations and indices.

    .. note::
      This class is also available as hifive.Fragment

    When initialized, this class creates an h5dict in which to store all data associated with this object.

    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`Fragment <hifive.fragment.Fragment>` class object.

    :Attributes: * **file** (*str.*) - A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcome.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a Fragment object."""
        self.file = filename
        self.silent = silent
        self.history = ""
        self.filetype = "fragment"
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
        Save fragment data to h5dict.

        :returns: None
        """
        self.history.replace("'None'", "None")
        fragfile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['file', 'silent']:
                continue
            elif isinstance(self[key], numpy.ndarray):
                fragfile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                fragfile.attrs[key] = self[key]
        fragfile.close()
        return None

    def load(self):
        """
        Load fragment data from h5dict specified at object creation.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :returns: None
        """
        fragfile = h5py.File(self.file, 'r')
        for key in fragfile.keys():
            self[key] = numpy.copy(fragfile[key])
        for key in fragfile['/'].attrs.keys():
            self[key] = fragfile['/'].attrs[key]
        fragfile.close()
        return None

    def load_fragments(self, filename, genome_name=None, re_name=None, regions=[],
                       minregionspacing=1000000):
        """
        Parse and store fragment data from a bed for a 5C assay file into an h5dict.

        :param filename: A file name to read restriction fragment data from. This should be a BED file containing fragment boundaries for all probed fragments and primer names that match those used for read mapping.
        :type filename: str.
        :param genome_name: The name of the species and build. Optional.
        :type genome_name: str.
        :param re_name: The name of the restriction enzyme used to produce the fragment set. Optional.
        :type re_name: str.
        :param regions: User-defined partitioning of fragments into different regions. This argument should be a list of lists containing the chromosome, start, and stop coordinates for each region.
        :type regions: list
        :param minregionspacing: If 'regions' is not defined, this is used to parse regions by inserting breaks where fragments are spaced apart greater than this value.
        :type minregionspacing: int.
        :returns: None

        :Attributes: * **chromosomes** (*ndarray*) - A numpy array containing chromosome names as strings. The position of the chromosome name in this array is referred to as the chromosome index.
                     * **fragments** (*ndarray*) - A numpy array of length N where N is the number of fragments and containing the fields 'chr', 'start', 'stop', 'mid', 'strand', 'region', and 'name'. With the exception of the 'name' field which is of type string, all of these are of type int32. The 'chr' and 'region' fields contain the indices of the chromosome and region, respectively. If the bed file used to create the Fragment object contains additional columns, these features are also included as fields with names corresponding to the bed header names. These additional fields are of type float32. Fragments are sorted by chromosome (the order in the 'chromosomes' array) and then by coordinates.
                     * **chr_indices** (*ndarray*) - A numpy array with a length of the number of chromosomes in 'chromosomes' + 1. This array contains the first position in 'fragments' for the chromosome in the corresponding position in the 'chromosomes' array. The last position in the array contains the total number of fragments.
                     * **regions** (*ndarray*) - A numpy array of length equal to the number of regions a containing the fields 'index', 'chromosome', 'start_frag', 'stop_frag', 'start' and 'stop'. Except for 'chromosome' which is a string, all fields are of type int32.
        """
        self.history += "Fragment.load_fragments(filename='%s', genome_name='%s', re_name='%s', regions=%s, minregionspacing=%i) - " % (filename, genome_name, re_name, str(regions), minregionspacing)
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.") % (filename),
            self.history += "Error: '%s' not found\n" % filename
            return None
        chromosomes = []
        chr2int = {}
        data = []
        input = open(filename, 'r')
        feature_names = []
        for line in input:
            temp = line.strip('\n').split('\t')
            try:
                start = int(temp[1])
            except:
                for i in range(6, len(temp)):
                    feature_names.append(temp[i])
                continue
            chrom = temp[0]
            start = int(temp[1])
            stop = int(temp[2])
            name = temp[3]
            if temp[5] == '+':
                strand = 0
            else:
                strand = 1
            features = []
            for i in range(6, 6 + len(feature_names)):
                features.append(float(temp[i]))
            if chrom not in chr2int:
                chr2int[chrom] = len(chromosomes)
                chromosomes.append(chrom)
            data.append([chr2int[chrom], start, stop, (start + stop) / 2, strand, name] + features)
        input.close()
        data.sort()
        dtypes = [('chr', numpy.int32), ('start', numpy.int32), ('stop', numpy.int32),
                  ('mid', numpy.int32), ('strand', numpy.int32), ('region', numpy.int32),
                  ('name', 'S32')]
        for feature in feature_names:
            dtypes.append((feature, numpy.float32))
        fragments = numpy.empty(len(data), dtype=numpy.dtype(dtypes))
        for i in range(len(data)):
            # make an entry for each fragment
            fragments['chr'][i] = data[i][0]
            fragments['start'][i] = data[i][1]
            fragments['stop'][i] = data[i][2]
            fragments['mid'][i] = data[i][3]
            fragments['strand'][i] = data[i][4]
            fragments['name'][i] = data[i][5]
            for j in range(len(feature_names)):
                fragments[feature_names[j]][i] = data[i][6 + j]
        chromosomes = numpy.array(chromosomes)
        # make note of chromosome positions in fend array
        chr_indices = numpy.zeros(chromosomes.shape[0] + 1, dtype=numpy.int32)
        chr_indices[1:] = numpy.bincount(fragments['chr'])
        for i in range(1, chr_indices.shape[0]):
            chr_indices[i] += chr_indices[i - 1]
        # if given by user, find region bounds
        region_array = []
        if len(regions) > 0:
            for region in regions:
                chrint = chr2int[region[0]]
                start_frag = chr_indices[chrint]
                while start_frag < chr_indices[chrint + 1] and fragments['mid'][start_frag] < region[1]:
                    start_frag += 1
                stop_frag = start_frag
                while stop_frag < chr_indices[chrint + 1] and fragments['mid'][stop_frag] < region[2]:
                    stop_frag += 1
                region_array.append([start_frag, stop_frag, chrint, region[1], region[2]])
        else:
            # if regions not given, parse regions
            for i in range(chr_indices.shape[0] - 1):
                start_frag = chr_indices[i]
                stop_frag = start_frag
                while start_frag < chr_indices[i + 1]:
                    while (stop_frag < chr_indices[i + 1] and (stop_frag == start_frag or
                           fragments['mid'][stop_frag] - fragments['mid'][stop_frag - 1] < minregionspacing)):
                        stop_frag += 1
                    region_array.append([i, start_frag, stop_frag, fragments['start'][start_frag],
                                        fragments['stop'][stop_frag - 1]])
                    start_frag = stop_frag
        # convert region list into array
        region_array.sort()
        regions = numpy.empty(len(region_array), dtype=numpy.dtype([
                              ('index', numpy.int32), ('start_frag', numpy.int32),
                              ('stop_frag', numpy.int32), ('chromosome', chromosomes.dtype),
                              ('start', numpy.int32), ('stop', numpy.int32)]))
        for i in range(len(region_array)):
            regions['index'][i] = i
            regions['chromosome'][i] = chromosomes[region_array[i][0]]
            regions['start_frag'][i] = region_array[i][1]
            regions['stop_frag'][i] = region_array[i][2]
            regions['start'][i] = region_array[i][3]
            regions['stop'][i] = region_array[i][4]
            fragments['region'][region_array[i][1]:region_array[i][2]] = i
        # write all data to h5dict
        self.fragments = fragments 
        if not re_name is None:
            self.re_name = re_name
        if not genome_name is None:
            self.genome_name = genome_name
        self.chr_indices = chr_indices
        self.chromosomes = chromosomes
        self.regions = regions
        self.history += 'Success\n'
        return None

if __name__ == '__main__':
    filename = sys.argv[1]
    fragment = Fragment('.'.join(filename.split('.')[:-1] + ['.frags']), 'w')
    fragment.load_fragments(filename)
    fragment.save()
