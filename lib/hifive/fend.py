#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

"""
This is a module class for handling restriction fragment data. The base
class "Fend" can load and save restriction fragment data for use with
"hifive" data storage classes "HiCData" and "FiveCData" and analysis classes
"HiCAnalysis" and "FiveCAnalysis".

Input data
----------

This class loads data from a text file formatted to work with "hicpipe".
The fragment occuring before the first and after the last restriction site
of each chromosome are discarded. Fragment ends (fends) are numbered starting
at zero, such that the fragment number for any given fend is simply the fend
number divided by two.

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


class Fend(object):
    """Base class for handling restriction fragment dataset.

    This class stores a list of chromosomes, a dictionary for converting from
    chromosome label to integer and back, fragment starts, stops, and
    chromosome number in an h5dict."""

    def __init__(self, filename, mode='r'):
        """
        __init__ method

        Initialize fend dataset and create an h5dict.

        Parameters
        ----------
        filename : string
            A filename to store h5dict in.
        mode : string, optional
            Specifies how to open the h5dict, depending on whether data is to
            be written or read.
        """
        self.fends = h5py.File(filename, mode)
        return None

    def load_fends(self, filename, genome_name=None, re_name=None):
        """
        load_fends method

        Parse and store fend data in h5dict.

        Parameters
        ----------
        filename : string
            A filename to read restriction fragment data from.
        genome_name : string, optional
            An indicator of genome species and build.
        re_name : string, optional
            An indicator of restriction enzyme.
        """
        if not os.path.exists(filename):
            print >> sys.stderr, ("Could not find %s. No data loaded.") % (filename),
            return None
        # if no genome name given, determine from filename
        if genome_name is None:
            genome_name = filename.split('/')[-1].split('_')[0]
        # if no re name given, determine from filename
        if re_name is None:
            re_name = filename.split('_')[-1].split('.')[0]
        chromosomes = []
        chr2int = {}
        data = {}
        input = open(filename, 'r')
        fragment_index = 1
        chromosome_index = 2
        coordinate_index = 3
        length_index = 5
        for line in input:
            temp = line.strip('\n').split('\t')
            if temp[0] == 'fend':
                fragment_index = temp.index('frag')
                chromosome_index = temp.index('chr')
                coordinate_index = temp.index('coord')
                length_index = temp.index('frag_len')
                continue
            chrom = temp[chromosome_index].strip('chr')
            frag = int(temp[fragment_index])
            if chrom not in chr2int:
                chr2int[chrom] = len(chromosomes)
                chromosomes.append(chrom)
            if frag not in data:
                # keep first instance of each fragment
                data[frag] = [0, chr2int[chrom], int(temp[coordinate_index]), int(temp[length_index])]
            else:
                # make sure that fragment has both ends present in file
                data[frag][0] += 1
        input.close()
        keys = data.keys()
        for key in keys:
            if data[key][0] != 1:
                del data[key]
        keys = data.keys()
        keys.sort()
        fends = numpy.empty(len(keys) * 2, dtype=numpy.dtype([('chr', numpy.int32), ('start', numpy.int32),
                ('stop', numpy.int32), ('mid', numpy.int32)]))
        for i, key in enumerate(keys):
            # make 2 entries for each fragment, one for each half
            index = i * 2
            fends['chr'][index] = data[key][1]
            fends['start'][index] = data[key][2]
            fends['stop'][index] = data[key][2] + data[key][3] / 2
            fends['mid'][index] = (fends['start'][index] + fends['stop'][index]) / 2
            index = i * 2 + 1
            fends['chr'][index] = data[key][1]
            fends['start'][index] = data[key][2] + data[key][3] / 2
            fends['stop'][index] = data[key][2] + data[key][3]
            fends['mid'][index] = (fends['start'][index] + fends['stop'][index]) / 2
        chromosomes = numpy.array(chromosomes)
        # make note of chromosome positions in fend array
        chr_indices = numpy.zeros(chromosomes.shape[0] + 1, dtype=numpy.int32)
        chr_indices[1:] = numpy.bincount(fends['chr'])
        for i in range(1, chr_indices.shape[0]):
            chr_indices[i] += chr_indices[i - 1]
        # write all data to h5dict
        dset = self.fends.create_dataset(name='fends', data=fends)
        dset.attrs['re_name'] = re_name
        dset.attrs['genome_name'] = genome_name
        self.fends.create_dataset(name='chr_indices', data=chr_indices)
        self.fends.create_dataset(name='chromosomes', data=chromosomes)
        return None

if __name__ == '__main__':
    filename = sys.argv[1]
    fend = Fend('.'.join(filename.split('.')[:-1] + ['.hdf5']), 'w')
    fend.load_fends(filename)
    fend.fends.close()
