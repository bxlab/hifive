#!/usr/bin/env python

import os
import sys
from math import ceil

import h5py
import numpy
try:
    import pyx
except:
    pass
try:
    from mpi4py import MPI
except:
    pass

import hic_binning as h_binning
import fivec_binning as f_binning
import libraries._bi as _bi


class BI(object):
    """
    This class uses a HiC or FiveC object to generate BI scores. These scores can be saved to an h5dict. Once generated or loaded, BI scores can be used to find structural boundaries, create feature-centered BI profiles, or bound-centered feature profiles.

    .. note::
      This class is also available as hifive.BI

    The BI score is calculated by finding the absolute difference between a set of matched bins; the first set contains interactions between sequences falling within a range given by 'window' up- and downstream of the boundary point and binned into 'height'-sized bins and the group of sequences falling between zero and 'width' base pairs upstream of the boundary point; the second set contains interactions between sequences falling within a range given by 'window' up- and downstream of the boundary point and binned into 'height'-sized bins and the group of sequences falling between zero and 'width' base pairs downstream of the boundary point.

    :param width: A range, in base-pairs, specifying how large a neighborhood around a boundary point is used for binning in calculating the BI.
    :type width: int.
    :param height: The size to bin fragments or fends into that fall within the range specified by 'window' of the boundary point and interacting with the fragments or fends follwing within the 'width' range around the boundary point. If set to zero, fragments or fends are left unbinned in this dimension for BI calculations.
    :type height: int.
    :param window: A range, in base-pairs, specifying how far away from a boundary point interactions are included in calculating the BI.
    :type window: int.
    :param mincount: The minimum number of valid pairs of bins needed to calculate the boundary index for a given location.
    :type mincount: int.
    """

    def __init__(self, width=10000, window=1000000, height=0, mincount=10):
        """
        Initialize :class:`BI` object.
        """
        self.width = int(width)
        self.window = int(window)
        self.height = int(height)
        self.mincount = int(mincount)
        return None

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        self.__dict__[key] = value
        return None

    def load(self, filename):
        """
        Load BI data and parameters from an h5dict.

        :param filename: This specifies the file name of the :class:`BI` object to load from.
        :type filename: str.
        """
        # ensure data h5dict exists
        if not os.path.exists(filename):
            print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename),
            return None
        print >> sys.stderr, ("Loading data..."),
        input = h5py.File(filename, 'r')
        for key in input.keys():
            self[key] = input[key][...]
        for key in input['/'].attrs.keys():
            self[key] = input['/'].attrs[key]
        input.close()
        if 'chromosomes' in self.__dict__:
            self.chr2int = {}
            for i in range(self.chromosomes.shape[0]):
                self.chr2int[self.chromosomes[i]] = i
        print >> sys.stderr, ("Done\n"),
        return None

    def save(self, filename):
        """
        Save BI data and parameters to an h5dict.

        :param filename: This specifies the file name to which to save the :class:`BI` object.
        :type filename: str.
        """
        # see if h5dict exists
        if os.path.exists(filename):
            print >> sys.stderr, ("%s already exists. Overwriting.\n") % (filename),
        print >> sys.stderr, ("Saving data..."),
        output = h5py.File(filename, 'w')
        for key in self.__dict__.keys():
            if isinstance(self[key], (basestring, int, float)):
                output.attrs[key] = self[key]
            elif not isinstance(self[key], dict):
                output.create_dataset(name=key, data=self[key])
        output.close()
        print >> sys.stderr, ("Done\n"),
        return None

    def find_bi_from_hic(self, hic, datatype='fend', chroms=[]):
        """
        Calculate the boundary index for HiC data. This function is MPI compatible.

        :param hic: A :class:`HiC` class object containing fend and count data.
        :type hic: :class:`HiC`
        :param chroms: A list of chromosome names to limit BI calculations to. If empty, all chromosomes in 'hic' are used.
        :type chroms: list
        """
        if 'mpi4py' in sys.modules.keys():
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()
        else:
            comm = None
            rank = 0
            num_procs = 1
        if rank == 0:
            print >> sys.stderr, ("Calculating BI..."),
        # save relevant parameters from hic object
        self.chromosomes = numpy.zeros(len(hic.chr2int), dtype=numpy.dtype('S20'))
        for chrom, chrint in hic.chr2int.iteritems():
            self.chromosomes[chrint] = chrom
        self.chr2int = dict(hic.chr2int)
        self.hicfilename = hic.file
        # determine which chromosomes to consider
        if len(chroms) == 0:
            chroms = list(hic.fends['chromosomes'])
        if len(chroms) == 0:
            chroms = numpy.copy(self.chromosomes)
        if rank == 0:
            worker_size = max(1, int(ceil(len(chroms) / float(num_procs))))
            for i in range(1, min(num_procs, len(chroms))):
                comm.send(chroms[(worker_size * (i - 1)):(worker_size * i)], dest=i, tag=11)
            for i in range(min(num_procs, len(chroms)), num_procs):
                comm.send([], dest=i, tag=11)
            worker_chroms = chroms[(min(num_procs - 1, len(chroms) - 1) * worker_size):]
        else:
            worker_chroms = comm.recv(source=0, tag=11)
        # for each chromosome, find BI
        BI = numpy.zeros(0, dtype=numpy.dtype([('chromosome', numpy.int32), ('start', numpy.int32),
                                               ('stop', numpy.int32), ('mid', numpy.int32), ('score', numpy.float32)]))
        temp_BI = {}
        for chrom in worker_chroms:
            # pull relevant data from hic object
            temp = h_binning.unbinned_cis_signal(hic, chrom, datatype=datatype, arraytype='compact',
                                                   maxdistance=(self.window + self.width), skipfiltered=True,
                                                   returnmapping=True)
            if temp is None or temp[1].shape[0] < 3:
                continue
            data, mapping = temp
            temp_BI[chrom] = numpy.zeros(mapping.shape[0] - 1, dtype=numpy.dtype([('chromosome', numpy.int32),
                                 ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32),
                                 ('score', numpy.float32)]))
            temp_BI[chrom]['chromosome'][:] = self.chr2int[chrom]
            temp_BI[chrom]['mid'][:] = (hic.fends['fends']['stop'][mapping[:-1]] +
                                 hic.fends['fends']['start'][mapping[1:]]) / 2
            temp_BI[chrom]['start'][:] = hic.fends['fends']['mid'][mapping[:-1]]
            temp_BI[chrom]['stop'][:] = hic.fends['fends']['mid'][mapping[1:]]
            if self.height == 0:
                temp = numpy.zeros((data.shape[1] * 2, 4), dtype=numpy.float32)
            else:
                temp = numpy.zeros((self.window / self.height * 2 + 2, 4), dtype=numpy.float32)
            scores = numpy.zeros(mapping.shape[0] - 1, dtype=numpy.float32)
            if self.height == 0:
                _bi.find_bi(data, hic.fends['fends']['mid'][mapping], temp_BI[chrom]['mid'], scores, temp, self.width,
                            self.window, self.mincount)
            else:
                _bi.find_bi_height(data, hic.fends['fends']['mid'][mapping], temp_BI[chrom]['mid'], scores, temp,
                                   self.width, self.window, self.height, self.mincount)
            temp_BI[chrom]['score'][:] = scores
            valid = numpy.where(scores != numpy.inf)[0]
            temp_BI[chrom] = temp_BI[chrom][valid]
        if rank == 0:
            for i in range(1, num_procs):
                temp_BI.update(comm.recv(source=i, tag=11))
            for chrom in chroms:
                BI = numpy.hstack((BI, temp_BI[chrom]))
            order = numpy.lexsort((BI['mid'], BI['chromosome']))
            BI = BI[order]
            for i in range(1, num_procs):
                comm.send(BI, dest=i, tag=11)
        else:
            comm.send(temp_BI, dest=0, tag=11)
            BI = comm.recv(source=0, tag=11)
        self.BI = BI
        # create an array containing the first index for each chromosome
        self.chr_indices = numpy.zeros(self.chromosomes.shape[0] + 1, dtype=numpy.int32)
        self.chr_indices[1:] = numpy.bincount(self.BI['chromosome'][:], minlength=self.chromosomes.shape[0])
        for i in range(1, self.chr_indices.shape[0]):
            self.chr_indices[i] += self.chr_indices[i - 1]
        if rank == 0:
            print >> sys.stderr, ("Done\n"),
        return None

    def find_bi_from_fivec(self, fivec, datatype='fragment', regions=[]):
        """
        Calculate the boundary index for 5C data.

        :param fivec: A :class:`FiveC` class object containing fend and count data.
        :type fivec: :class:`FiveC`
        :param regions: A list of region numbers to limit BI calculations to. If empty, all regions in 'fivec' are used.
        :type regions: list
        """
        print >> sys.stderr, ("Calculating BI..."),
        # save relevant parameters from hic object
        self.chromosomes = numpy.zeros(len(fivec.chr2int), dtype=numpy.dtype('S20'))
        for chrom, chrint in fivec.chr2int.iteritems():
            self.chromosomes[chrint] = chrom
        self.chr2int = dict(fivec.chr2int)
        self.fivecfilename = fivec.file
        # determine which chromosomes to consider
        if len(regions) == 0:
            regions = range(fivec.frags['regions'].shape[0])
        # for each region, find BI
        BI = numpy.zeros(0, dtype=numpy.dtype([('chromosome', numpy.int32), ('start', numpy.int32),
                                               ('stop', numpy.int32), ('mid', numpy.int32), ('score', numpy.float32)]))
        for region in regions:
            # pull relevant data from fivec object
            temp = f_binning.unbinned_cis_signal(fivec, region, datatype=datatype, arraytype='full',
                                                     skipfiltered=True, returnmapping=True)
            if temp is None or temp[1].shape[0] < 3:
                continue
            start_frag = fivec.frags['regions']['start_frag'][region]
            stop_frag = fivec.frags['regions']['stop_frag'][region]
            distances = (fivec.frags['fragments']['mid'][start_frag:stop_frag].reshape(1, -1) - 
                         fivec.frags['fragments']['mid'][start_frag:stop_frag].reshape(-1, 1))
            mapping = temp[1]
            max_bin = 0
            for i in range(temp[0].shape[0] - 1):
                where = numpy.where(distances[i, (i + 1):] > self.window)[0]
                if where.shape[0] > 0:
                    max_bin = max(max_bin, where[0])
                    temp[0][i, (i + 1 + where[0]):, :] = 0
            data = numpy.zeros((temp[0].shape[0], max_bin, 2), dtype=numpy.float32)
            for i in range(temp[0].shape[0] - 1):
                data[i, :min(data.shape[0] - i - 1, max_bin), :] = temp[0][i, (i + 1):(i + 1 + min(data.shape[0] -
                                                                                       i - 1, max_bin)), :]
            temp_BI = numpy.zeros(mapping.shape[0] - 1, dtype=numpy.dtype([('chromosome', numpy.int32),
                                 ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32),
                                 ('score', numpy.float32)]))
            temp_BI['chromosome'][:] = self.chr2int[chrom]
            temp_BI['mid'][:] = (fivec.frags['fragments']['stop'][mapping[:-1]] +
                                 fivec.frags['fragments']['start'][mapping[1:]]) / 2
            temp_BI['start'][:] = fivec.frags['fragments']['mid'][mapping[:-1]]
            temp_BI['stop'][:] = fivec.frags['fragments']['mid'][mapping[1:]]
            if self.height == 0:
                temp = numpy.zeros((data.shape[1] * 2, 4), dtype=numpy.float32)
            else:
                temp = numpy.zeros((self.window / self.height * 2 + 2, 4), dtype=numpy.float32)
            scores = numpy.zeros(mapping.shape[0] - 1, dtype=numpy.float32)
            if self.height == 0:
                _bi.find_bi(data, fivec.frags['fragments']['mid'][mapping], temp_BI['mid'], scores, temp,
                           self.width, self.mincount)
            else:
                _bi.find_bi_height(data, fivec.frags['fragments']['mid'][mapping], temp_BI['mid'], scores, temp,
                                   self.width, self.window, self.height, self.mincount)
            temp_BI['score'][:] = -scores
            valid = numpy.where(scores != numpy.inf)[0]
            if valid.shape[0] > 0:
                BI = numpy.hstack((BI, temp_BI[valid]))
        order = numpy.lexsort((BI['mid'], BI['chromosome']))
        BI = BI[order]
        self.BI = BI
        # create an array containing the first index for each chromosome
        self.chr_indices = numpy.zeros(self.chromosomes.shape[0] + 1, dtype=numpy.int32)
        self.chr_indices[1:] = numpy.bincount(self.BI['chromosome'][:], minlength=self.chromosomes.shape[0])
        for i in range(1, self.chr_indices.shape[0]):
            self.chr_indices[i] += self.chr_indices[i - 1]
        print >> sys.stderr, ("Done\n"),
        return None

    def smooth_bi(self, width):
        """
        Using a Gaussian weighting approach, smooth BI values. Weights extend out to 2.5 standard deviatons. This function is MPI compatible.

        :param width: The distance, in base pairs, equivalent to one standard deviation for the Gaussian weights.
        :type width: int.
        """
        if 'mpi4py' in sys.modules.keys():
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()
        else:
            comm = None
            rank = 0
            num_procs = 1
            self.smoothing = width
        if rank == 0:
            print >> sys.stderr, ("Finding smoothed BIs..."),
        chroms = list(self.chromosomes)
        # if needed, get list of chromosomes
        if len(chroms) == 0:
            chroms = numpy.copy(self.chromosomes)
        adjusted = {}
        if rank == 0:
            worker_size = max(1, int(ceil(len(chroms) / float(num_procs))))
            for i in range(1, min(num_procs, len(chroms))):
                comm.send(chroms[(worker_size * (i - 1)):(worker_size * i)], dest=i, tag=11)
            for i in range(min(num_procs, len(chroms)), num_procs):
                comm.send([], dest=i, tag=11)
            worker_chroms = chroms[(min(num_procs - 1, len(chroms) - 1) * worker_size):]
        else:
            worker_chroms = comm.recv(source=0, tag=11)
        for chrom in worker_chroms:
            chrint = self.chr2int[chrom]
            if self.chr_indices[chrint] + 3 >= self.chr_indices[chrint + 1]:
                continue
            original = numpy.copy(self.BI['score'][self.chr_indices[chrint]:self.chr_indices[chrint + 1]])
            mids = self.BI['mid'][self.chr_indices[chrint]:self.chr_indices[chrint + 1]]
            smoothed = numpy.copy(original)
            _bi.gaussian_smoothing(mids, original, smoothed, width)
            adjusted[chrom] = smoothed
        if rank == 0:
            for i in range(1, num_procs):
                adjusted.update(comm.recv(source=i, tag=11))
            BI = numpy.zeros(self.BI.shape[0], dtype=numpy.dtype([('chromosome', numpy.int32), ('start', numpy.int32),
                 ('stop', numpy.int32), ('mid', numpy.int32), ('score', numpy.float32), ('original', numpy.float32)]))
            BI['chromosome'] = self.BI['chromosome'][:]
            BI['mid'] = self.BI['mid'][:]
            BI['start'] = self.BI['start'][:]
            BI['stop'] = self.BI['stop'][:]
            BI['original'] = self.BI['score'][:]
            for chrom in adjusted:
                chrint = self.chr2int[chrom]
                BI['score'][self.chr_indices[chrint]:self.chr_indices[chrint + 1]] = adjusted[chrom][:]
            for i in range(1, num_procs):
                comm.send(BI, dest=i, tag=11)
        else:
            comm.send(adjusted, dest=0, tag=11)
            BI = comm.recv(source=0, tag=11)
        self.BI = BI
        if rank == 0:
            print >> sys.stderr, ("Done\n"),
        return None

    def plot_bi(self, width, height, chromosome, start=None, stop=None):
        """
        Plot a line representation of the BI score for a specified region.

        :param width: The width of the image to return in whatever units :class:`Pyx` is using.
        :type width: float
        :param height: The height of the image to return in whatever units :class:`Pyx` is using.
        :type height: float
        :param chromosome: The name of the chromosome to pull BI information from.
        :type chromosome: str.
        :param start: The coordinates corresponding to the left edge of the image. If not specified, the first BI midpoint is used.
        :type start: int.
        :param stop: The coordinates corresponding to the right edge of the image. If not specified, the first BI midpoint is used.
        :type stop: int.
        """
        if 'pyx' not in sys.modules.keys():
            print >> sys.stderr, ("The pyx module must be installed to use this function.")
            return None
        if 'BI' not in self.__dict__:
            print >> sys.stderr, ("No BI scores present, no plot generated.\n"),
            return None
        if chromosome not in self.chr2int:
            print >> sys.stderr, ("Unrecognized chromosome name, no plot generated.\n"),
            return None
        chrint = self.chr2int[chromosome]
        if start == None:
            start = self.BI['mid'][self.chr_indices[chrint]]
        if stop == None:
            stop = self.BI['mid'][self.chr_indices[chrint + 1] - 1]
        print >> sys.stderr, ("Plotting BIs for %s:%i-%i...") % (chromosome, start, stop),
        c = pyx.canvas.canvas()
        start_index = numpy.searchsorted(self.BI['mid'][self.chr_indices[chrint]:self.chr_indices[chrint + 1]], start)
        start_index += self.chr_indices[chrint]
        stop_index = numpy.searchsorted(self.BI['mid'][self.chr_indices[chrint]:self.chr_indices[chrint + 1]], stop)
        stop_index += self.chr_indices[chrint]
        X = (self.BI['mid'][start_index:stop_index] - start) / float(stop - start) * width
        Y = self.BI['score'][start_index:stop_index]
        Y /= numpy.amax(numpy.abs(Y))
        Y = (Y + 1) * height / 2.0
        line_path = pyx.path.path(pyx.path.moveto(0, height / 2.0))
        for i in range(X.shape[0]):
            line_path.append(pyx.path.lineto(X[i], Y[i]))
        line_path.append(pyx.path.lineto(width, height / 2.0))
        c.fill(line_path)
        print >> sys.stderr, ("Done\n"),
        return c

    def find_bi_profile(self, coordinates, binsize=2000, numbins=51):
        """
        Find the mean BI signal in a number of bins centered around a set of user-specified positions.

        :param coordinates: A dictionary keyed with chromosome names and containing 1d arrays of coordinates for each chromosome.
        :type coordinates: dict.
        :param binsize: Width, in base pairs, of each bin to aggregate signal in.
        :type binsize: int.
        :param numbins: The number of bins, centered around each coordinate, to aggregate signal in.
        :type numbins: int.
        """
        print >> sys.stderr, ("Finding BI profile..."),
        signal = numpy.zeros((numbins, 2), dtype=numpy.float32)
        for chrom in coordinates:
            if chrom not in self.chr2int:
                continue
            chrint = self.chr2int[chrom]
            sorted_coords = numpy.copy(coordinates[chrom])
            sorted_coords.sort()
            _bi.find_bi_profile(self.BI['start'], self.BI['stop'], self.BI['score'], self.chr_indices,
                                sorted_coords, signal, binsize, chrint)
        where = numpy.where(signal[:, 1] > 0.0)
        signal[where, 0] /= signal[where, 1]
        print >> sys.stderr, ("Done\n"),
        return signal[:, 0]

    def find_bi_bounds(self, cutoff=2.0, chroms=[]):
        """
        Find BI scores that are higher than adjacent scores and are more than the cutoff * standard deviation above the minimum adjacent trough score.

        :param cutoff: Minimum height ratio between peak and lower adjacent trough.
        :type cutoff: float
        :param chroms: If specified, find bounds for only chromosomes in this list. Otherwise use all chromosomes present in the :class:`BI` object.
        :type chroms: list
        """
        print >> sys.stderr, ("Finding BI bounds..."),
        # if needed, get list of chromosomes
        if len(chroms) == 0:
            chroms = numpy.copy(self.chromosomes)
        adjusted = numpy.copy(self.BI['score'])
        for chrom in chroms:
            chrint = self.chr2int[chrom]
            if self.chr_indices[chrint] + 3 >= self.chr_indices[chrint + 1]:
                continue
            # if needed, smooth BI line
            original = numpy.copy(self.BI['score'][self.chr_indices[chrint]:self.chr_indices[chrint + 1]])
            original -= numpy.mean(original)
            adjusted[self.chr_indices[chrint]:self.chr_indices[chrint + 1]] = original
        adjusted /= numpy.std(adjusted)
        bounds = numpy.zeros(0, dtype=numpy.dtype([('chr', numpy.int32), ('coord', numpy.int32),
                             ('score', numpy.float32)]))
        for chrom in chroms:
            chrint = self.chr2int[chrom]
            if self.chr_indices[chrint] + 3 >= self.chr_indices[chrint + 1]:
                continue
            scores = adjusted[self.chr_indices[chrint]:self.chr_indices[chrint +1]]
            coords = self.BI['mid'][self.chr_indices[chrint]:self.chr_indices[chrint + 1]]
            peaks = numpy.where((scores[1:-1] > scores[2:]) * (scores[1:-1] > scores[:-2]))[0] + 1
            valid = []
            if scores[peaks[0]] - numpy.amin(scores[:peaks[0]]) >= cutoff:
                valid.append(0)
            for i in range(1, peaks.shape[0] - 1):
                if scores[peaks[i]] - numpy.amin(scores[(peaks[i - 1] + 1):peaks[i + 1]]) >= cutoff:
                    valid.append(i)
            if scores[peaks[-1]] - numpy.amin(scores[(peaks[-1] + 1):]) >= cutoff:
                valid.append(peaks.shape[0] - 1)
            peaks = peaks[valid]
            new_bounds = numpy.zeros(peaks.shape[0], dtype=numpy.dtype([('chr', numpy.int32),
                                     ('coord', numpy.int32), ('score', numpy.float32)]))
            new_bounds['chr'][:] = chrint
            new_bounds['coord'][:] = self.BI['mid'][self.chr_indices[chrint] + peaks]
            new_bounds['score'][:] = scores[peaks]
            print "Chr %s   %i out of %i" % (chrom, new_bounds.shape[0], scores.shape[0])
            bounds = numpy.hstack((bounds, new_bounds))
        print >> sys.stderr, ("Done\n"),
        return bounds
