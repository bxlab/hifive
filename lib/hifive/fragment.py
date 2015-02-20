#!/usr/bin/env python

import os
import sys

import numpy
import h5py


class Fragment(object):
    """
    This class handles restriction enzyme digest-generated fragment data for 5C experiments.

    This class stores a list of chromosomes, a dictionary for converting from chromosome label to integer and back, fragment starts, stops, and chromosome number in an h5dict.

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
    """

    def __init__(self, filename, mode='r', silent=False):
        """
        Create a Fragment object.
        """
        self.fragments = h5py.File(filename, mode)
        self.silent = silent
        return None

    def save(self):
        """
        Save fragment data to h5dict.

        :returns: None
        """
        self.fragments.close()
        return None

    def load_fragments(self, filename, genome_name=None, re_name=None, regions=[], minregionspacing=1000000):
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
        """
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.") % (filename),
            return None
        chromosomes = []
        chr2int = {}
        data = []
        input = open(filename, 'r')
        fragment_index = 1
        chromosome_index = 2
        coordinate_index = 3
        length_index = 5
        for line in input:
            temp = line.strip('\n').split('\t')
            if temp[0] == 'chr':
                continue
            chrom = temp[0].strip('chr')
            start = int(temp[1])
            stop = int(temp[2])
            name = temp[3]
            if temp[5] == '+':
                strand = 0
            else:
                strand = 1
            if chrom not in chr2int:
                chr2int[chrom] = len(chromosomes)
                chromosomes.append(chrom)
            data.append([chr2int[chrom], start, stop, (start + stop) / 2, strand, name])
        input.close()
        data.sort()
        fragments = numpy.empty(len(data), dtype=numpy.dtype([('chr', numpy.int32), ('start', numpy.int32),
                ('stop', numpy.int32), ('mid', numpy.int32), ('strand', numpy.int32), ('name', 'S32')]))
        for i in range(len(data)):
            # make an entry for each fragment
            fragments['chr'][i] = data[i][0]
            fragments['start'][i] = data[i][1]
            fragments['stop'][i] = data[i][2]
            fragments['mid'][i] = data[i][3]
            fragments['strand'][i] = data[i][4]
            fragments['name'][i] = data[i][5]
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
                chrint = chr2int[region[0].strip('chr')]
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
        # write all data to h5dict
        self.fragments.create_dataset(name='fragments', data=fragments)
        if not re_name is None:
            self.fragments.attrs['re_name'] = re_name
        if not genome_name is None:
            self.fragments.attrs['genome_name'] = genome_name
        self.fragments.create_dataset(name='chr_indices', data=chr_indices)
        self.fragments.create_dataset(name='chromosomes', data=chromosomes)
        self.fragments.create_dataset(name='regions', data=regions)
        return None

if __name__ == '__main__':
    filename = sys.argv[1]
    fragment = Fragment('.'.join(filename.split('.')[:-1] + ['.hdf5']), 'w')
    fragment.load_fragments(filename)
    fragment.fragments.close()
