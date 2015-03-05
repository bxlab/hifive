#!/usr/bin/env python

import os
import sys

import numpy
import h5py


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
    """

    def __init__(self, filename, mode='r', silent=False):
        """
        Create a Fend object.
        """
        self.fends = h5py.File(filename, mode)
        self.silent = silent
        return None

    def save(self):
        """
        Save fend data to h5dict.

        :returns: None
        """
        self.fends.close()
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
        """
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.") % (filename),
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
            return None
        # make note of chromosome positions in fend array
        chr_indices = numpy.zeros(chromosomes.shape[0] + 1, dtype=numpy.int32)
        chr_indices[1:] = numpy.bincount(fends['chr'])
        for i in range(1, chr_indices.shape[0]):
            chr_indices[i] += chr_indices[i - 1]
        # write all data to h5dict
        self.fends.create_dataset(name='fends', data=fends)
        self.fends.attrs['re_name'] = re_name
        self.fends.attrs['genome_name'] = genome_name
        self.fends.create_dataset(name='chr_indices', data=chr_indices)
        self.fends.create_dataset(name='chromosomes', data=chromosomes)
        return None

    def _load_from_fend(self, fname):
        chromosomes = []
        chr2int = {}
        data = {}
        gcscore = {}
        mapscore = {}
        input = open(fname, 'r')
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
                if 'frag_gc' in temp:
                    gc_index = temp.index('frag_gc')
                else:
                    gc_index = None
                if 'map_score' in temp:
                    map_index = temp.index('map_score')
                else:
                    map_index = None
                continue
            chrom = temp[chromosome_index].strip('chr')
            frag = int(temp[fragment_index])
            if chrom not in chr2int:
                chr2int[chrom] = len(chromosomes)
                chromosomes.append(chrom)
            if frag not in data:
                # keep first instance of each fragment
                data[frag] = [0, chr2int[chrom], int(temp[coordinate_index]), int(temp[length_index])]
                if not gc_index is None:
                    gcscore[(frag, 0)] = float(temp[gc_index])
                if not map_index is None:
                    mapscore[(frag, 0)] = float(temp[map_index])
            else:
                # make sure that fragment has both ends present in file
                data[frag][0] += 1
                if not gc_index is None:
                    gcscore[(frag, 1)] = float(temp[gc_index])
                if not map_index is None:
                    mapscore[(frag, 1)] = float(temp[map_index])
        input.close()
        keys = data.keys()
        for key in keys:
            if data[key][0] != 1:
                del data[key]
        keys = data.keys()
        keys.sort()
        dtypes = [('chr', numpy.int32), ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32)]
        if len(gcscore) > 0:
            dtypes.append(('gc', numpy.float32))
        if len(mapscore) > 0:
            dtypes.append(('mappability', numpy.float32))
        fends = numpy.empty(len(keys) * 2, dtype=numpy.dtype(dtypes))
        for i, key in enumerate(keys):
            # make 2 entries for each fragment, one for each half
            index = i * 2
            fends['chr'][index] = data[key][1]
            fends['start'][index] = data[key][2]
            fends['stop'][index] = data[key][2] + data[key][3] / 2
            fends['mid'][index] = (fends['start'][index] + fends['stop'][index]) / 2
            if (key, 0) in gcscore:
                fends['gc'][index] = gcscore[(key, 0)]
            if (key, 0) in mapscore:
                fends['mappability'][index] = mapscore[(key, 0)]
            index = i * 2 + 1
            fends['chr'][index] = data[key][1]
            fends['start'][index] = data[key][2] + data[key][3] / 2
            fends['stop'][index] = data[key][2] + data[key][3]
            fends['mid'][index] = (fends['start'][index] + fends['stop'][index]) / 2
            if (key, 1) in gcscore:
                fends['gc'][index] = gcscore[(key, 1)]
            if (key, 1) in mapscore:
                fends['mappability'][index] = mapscore[(key, 1)]
        chromosomes = numpy.array(chromosomes)
        return [fends, chromosomes]

    def _load_from_bed(self, fname):
        chromosomes = []
        chr2int = {}
        data = {}
        input = open(fname, 'r')
        for line in input:
            temp = line[:-1].split('\t')
            if temp[0] == 'chr':
                continue
            chrom = temp[0].strip('chr')
            if chrom not in data:
                data[chrom] = []
            data[chrom].append([int(temp[1]), int(temp[2])])
        input.close()
        chromosomes = data.keys()
        chromosomes.sort()
        for i, chrom in enumerate(chromosomes):
            data[chrom] = numpy.array(data[chrom], dtype=numpy.int32)
            if numpy.mean(data[chrom][:, 1] - data[chrom][:, 0]) < 10.0:
                data[chrom] = (data[chrom][:, 0] + data[chrom][:, 1]) / 2
            else:
                data[chrom] = numpy.r_[data[chrom][0, 0], data[chrom][:, 1]].astype(numpy.int32)
            data[chrom] = data[chrom][numpy.argsort(data[chrom])]
        indices = numpy.zeros(len(chromosomes) + 1, dtype=numpy.int32)
        for i, chrom in enumerate(chromosomes):
            indices[i + 1] = indices[i] + (data[chrom].shape[0] - 1) * 2
        fends = numpy.empty(indices[-1], dtype=numpy.dtype([('chr', numpy.int32), ('start', numpy.int32),
                ('stop', numpy.int32), ('mid', numpy.int32)]))
        for i, chrom in enumerate(chromosomes):
            # make 2 entries for each fragment, one for each half
            fends['chr'][indices[i]:indices[i + 1]] = i
            fends['start'][indices[i]:indices[i + 1]:2] = data[chrom][:-1]
            fends['start'][(indices[i] + 1):indices[i + 1]:2] = (data[chrom][:-1] + data[chrom][1:]) / 2
            fends['stop'][indices[i]:indices[i + 1]:2] = (data[chrom][:-1] + data[chrom][1:]) / 2
            fends['stop'][(indices[i] + 1):indices[i + 1]:2] = data[chrom][1:]
            fends['mid'][indices[i]:indices[i + 1]] = (fends['start'][indices[i]:indices[i + 1]] +
                                                       fends['stop'][indices[i]:indices[i + 1]]) / 2
        chromosomes = numpy.array(chromosomes)
        return [fends, chromosomes]


if __name__ == '__main__':
    filename = sys.argv[1]
    fend = Fend('.'.join(filename.split('.')[:-1] + ['.hdf5']), 'w')
    fend.load_fends(filename)
    fend.fends.close()
