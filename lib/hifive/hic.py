#!/usr/bin/env python

import os
import sys
from math import ceil, floor, log, exp

import numpy
import h5py
try:
    from mpi4py import MPI
except:
    pass
try:
    import mlpy
except:
    pass

from fend import Fend
from hic_data import HiCData
from libraries._hic_binning import find_fend_coverage, dynamically_bin_unbinned_upper, dynamically_bin_upper_from_upper, remap_counts
import hic_binning
import libraries._hic_distance as _distance
import plotting


class HiC(object):
    """
    This is the class for handling HiC analysis.

    This class relies on :class:`Fend <hifive.fend.Fend>` and :class:`HiCData <hifive.hic_data.HiCData>` for genomic position and interaction count data. Use this class to perform filtering of fends based on coverage, model fend bias and distance dependence, and downstream analysis and manipulation. This includes binning of data, plotting of data, modeling of data, and statistical analysis.

        .. note::
          This class is also available as hifive.HiC

    When initialized, this class creates an h5dict in which to store all data associated with this object.

    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`HiC <hifive.hic.HiC>` class object.
    """

    def __init__(self, filename, mode='r', silent=False):
        """
        Create a HiC object.
        """

        self.file = os.path.abspath(filename)
        if 'mpi4py' in sys.modules.keys():
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.num_procs = self.comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.num_procs = 1
        if self.rank == 0:
            self.silent = silent
        else:
            self.silent = True
        if mode != 'w':
            self.load()
        return None

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        self.__dict__[key] = value
        return None

    def load_data(self, filename):
        """
        Load fend-pair counts and fend object from :class:`HiCData <hifive.hic_data.HiCData>` object.

        :param filename: Specifies the file name of the :class:`HiCData <hifive.hic_data.HiCData>` object to associate with this analysis.
        :type filename: str.
        :returns: None
        """
        if self.rank > 0:
            return None
        filename = os.path.abspath(filename)
        # ensure data h5dict exists
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename),
            return None
        self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                       os.path.dirname(self.file)), os.path.basename(filename))
        self.data = h5py.File(filename, 'r')
        fendfilename = self.data['/'].attrs['fendfilename']
        if fendfilename[:2] == './':
            fendfilename = fendfilename[2:]
        parent_count = fendfilename.count('../')
        fendfilename = '/'.join(os.path.abspath(filename).split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(fendfilename),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        # ensure fend h5dict exists
        if not os.path.exists(fendfilename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s.\n") % (fendfilename),
            return None
        self.fends = h5py.File(fendfilename, 'r')
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.fends['chromosomes']):
            self.chr2int[chrom] = i
        # create arrays
        self.filter = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.int32)
        self.corrections = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.float32)
        return None

    def save(self, out_fname=None):
        """
        Save analysis parameters to h5dict.

        :param filename: Specifies the file name of the :class:`HiC <hifive.hic.HiC>` object to save this analysis to.
        :type filename: str.
        :returns: None
        """
        if self.rank > 0:
            return None
        if not out_fname is None:
            original_file = os.path.abspath(self.file)
            self.file = out_fname
            if 'datafilename' in self.__dict__:
                datafilename = self.datafilename
                if datafilename[:2] == './':
                    datafilename = datafilename[2:]
                parent_count = datafilename.count('../')
                datafilename = '/'.join(original_file.split('/')[:-(1 + parent_count)] +
                                        datafilename.lstrip('/').split('/')[parent_count:])
                self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(datafilename)),
                                               os.path.dirname(self.file)), os.path.basename(datafilename))
            if 'fendfilename' in self.__dict__:
                fendfilename = self.fendfilename
                if fendfilename[:2] == './':
                    fendfilename = fendfilename[2:]
                parent_count = fendfilename.count('../')
                fendfilename = '/'.join(original_file.split('/')[:-(1 + parent_count)] +
                                        fendfilename.lstrip('/').split('/')[parent_count:])
                self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                               os.path.dirname(self.file)), os.path.basename(fendfilename))
        datafile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['data', 'fends', 'file', 'chr2int', 'comm', 'rank', 'num_procs', 'silent']:
                continue
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                datafile.attrs[key] = self[key]
        datafile.close()
        return None

    def load(self):
        """
        Load analysis parameters from h5dict specified at object creation and open h5dicts for associated :class:`HiCData <hifive.hic_data.HiCData>` and :class:`Fend <hifive.fend.Fend>` objects.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :returns: None
        """
        datafile = h5py.File(self.file, 'r')
        for key in datafile.keys():
            self[key] = numpy.copy(datafile[key])
        for key in datafile['/'].attrs.keys():
            self[key] = datafile['/'].attrs[key]
        # ensure data h5dict exists
        if 'datafilename' in self.__dict__:
            datafilename = self.datafilename
            if datafilename[:2] == './':
                datafilename = datafilename[2:]
            parent_count = datafilename.count('../')
            datafilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                datafilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(datafilename):
                if not self.silent:
                    print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (datafilename),
            else:
                self.data = h5py.File(datafilename, 'r')
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

    def reset_filter(self):
        """
        Return all fends to a valid filter state.

        :returns: None
        """
        self.filter.fill(1)
        return None

    def filter_fends(self, mininteractions=10, mindistance=0, maxdistance=0):
        """
        Iterate over the dataset and remove fends that do not have 'minobservations' within 'maxdistance' of themselves using only unfiltered fends.

        In order to create a set of fends that all have the necessary number of interactions, after each round of filtering, fend interactions are retallied using only interactions that have unfiltered fends at both ends.

        :param mininteractions: The required number of interactions for keeping a fend in analysis.
        :type mininteractions: int.
        :param mindistance: The minimum inter-fend distance used to count fend interactions.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fend distance used to count fend interactions. A value of 0 indicates no maximum should be used.
        :type maxdistance: int.
        :returns: None
        """
        if self.rank > 0:
            return None
        if not self.silent:
            print >> sys.stderr, ("Filtering fends..."),
        self.mininteractions = mininteractions
        self.mindistance = mindistance
        self.maxdistance = maxdistance
        original_count = numpy.sum(self.filter)
        previous_valid = original_count + 1
        current_valid = original_count
        coverage = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        # determine maximum ranges of valid interactions for each fend
        max_fend = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        min_fend = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        _distance.find_max_fend(max_fend,
                                self.fends['fends']['mid'][...],
                                self.fends['fends']['chr'][...],
                                self.fends['chr_indices'][...],
                                0,
                                maxdistance)
        _distance.find_min_fend(min_fend,
                        self.fends['fends']['mid'][...],
                        self.fends['fends']['chr'][...],
                        self.fends['chr_indices'][...],
                        mindistance)
         # copy needed arrays
        data = self.data['cis_data'][...]
        indices = self.data['cis_indices'][...]
        # repeat until all remaining fends have mininteraction valid interactions
        while current_valid < previous_valid:
            previous_valid = current_valid
            coverage.fill(0)
            find_fend_coverage(data,
                               indices,
                               self.filter,
                               min_fend,
                               max_fend,
                               coverage,
                               mininteractions)
            current_valid = numpy.sum(self.filter)
        if not self.silent:
            print >> sys.stderr, ("Removed %i of %i fends\n") % (original_count - current_valid, original_count),
        return None

    def find_distance_parameters(self, numbins=90, minsize=200, maxsize=0, corrected=False):
        """
        Count reads and possible interactions from valid fend pairs in each distance bin to find mean bin signals. This function is MPI compatible.

        This partitions the range of interaction distances (measured from mipoints of the involved fends) from the 'minsize' to 'maxsize' into a number of partitions equal to 'numbins'. The first bin contains all distances less than or equal to 'minsize'. The remaining bins are defined such that their log ranges are equal to one another. The curve defined by the mean interaction value of each bin can be smoothed using a  triangular smoothing operation.

        :param numbins: The number of bins to divide the distance range into. The first bin extends from zero to 'minsize', while the remaining bins are divided into evenly-spaced log-sized bins from 'minsize' to 'maxsize' or the maximum inter-fend distance, whichever is greater.
        :type numbins: int.
        :param minsize: The upper size limit of the smallest distance bin.
        :type minsize: int.
        :param maxsize: If this value is larger than the largest included chromosome, it will extend bins out to maxsize. If this value is smaller, it is ignored.
        :type maxsize: int.
        :param corrected: If True, correction values are applied to counts prior to summing.
        :type corrected: bool.
        :returns: None
        """
        if not self.silent:
            print >> sys.stderr, ('Finding distance arrays...'),
        # determine max distance range of data
        chr_indices = self.fends['chr_indices'][...]
        max_dist = 0
        valid_chroms = []
        for i in range(chr_indices.shape[0] - 1):
            start_fend = chr_indices[i]
            stop_fend = chr_indices[i + 1]
            while start_fend < stop_fend and self.filter[start_fend] == 0:
                start_fend += 1
            while stop_fend > start_fend and self.filter[stop_fend - 1] == 0:
                stop_fend -= 1
            if stop_fend - 1 > start_fend:
                valid_chroms.append(i)
                max_dist = max(max_dist,
                               self.fends['fends']['mid'][stop_fend - 1] - self.fends['fends']['mid'][start_fend])
        valid_chroms = numpy.asarray(valid_chroms, dtype=numpy.int32)
        num_chroms = valid_chroms.shape[0]
        # create cutoff values evenly spaced in log space
        cutoffs = numpy.linspace(numpy.log(max(minsize, 1.0)), numpy.log(max(maxsize, max_dist)),
                                 numbins).astype(numpy.float32)
        cutoffs[-1] += 1.0
        bin_size = numpy.zeros((num_chroms + 1, numbins), dtype=numpy.int64)
        count_sum = numpy.zeros((num_chroms + 1, numbins), dtype=numpy.float64)
        logdistance_sum = numpy.zeros((num_chroms + 1, numbins), dtype=numpy.float64)
        if 'binned' in self.data.attrs.keys() and self.data.attrs['binned'] == True:
            binned = 1
        else:
            binned = 0
        # for each chromosome, find counts, possible interactions, and distance sums for each bin
        for h, i in enumerate(valid_chroms):
            start_fend = chr_indices[i]
            stop_fend = chr_indices[i + 1]
            rev_mapping = numpy.where(self.filter[start_fend:stop_fend] == 1)[0].astype(numpy.int32)
            num_valid = rev_mapping.shape[0]
            if not self.silent:
                print >> sys.stderr, ('\r%s\rFinding distances for chromosome %s...') % \
                                     (' ' * 80, self.fends['chromosomes'][i]),
            # partition total possible interactions into roughly even-sized groups to spread across nodes
            total_pairs = num_valid * (num_valid - 1) / 2
            node_cutoffs = numpy.linspace(0, total_pairs, self.num_procs + 1).astype(numpy.float64)
            pair_sums = numpy.r_[0, num_valid - numpy.arange(1, num_valid + 1)].astype(numpy.int64)
            for j in range(1, pair_sums.shape[0]):
                pair_sums[j] += pair_sums[j - 1]
            node_ranges = numpy.searchsorted(pair_sums, node_cutoffs, side='left')
            node_start = node_ranges[self.rank]
            if (node_start > 0 and node_cutoffs[self.rank] - pair_sums[node_start - 1] <
                                   pair_sums[node_start] - node_cutoffs[self.rank]):
                node_start -= 1
            node_stop = node_ranges[self.rank + 1]
            if (node_stop > 0 and node_cutoffs[self.rank + 1] - pair_sums[node_stop - 1] <
                                   pair_sums[node_stop] - node_cutoffs[self.rank + 1]):
                node_stop -= 1
            mapping = numpy.zeros(stop_fend - start_fend, dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(num_valid, dtype=numpy.int32)
            mids = self.fends['fends']['mid'][rev_mapping + start_fend]
            # pull relevant data
            start_index = self.data['cis_indices'][rev_mapping[node_start] + start_fend]
            stop_index = self.data['cis_indices'][max(rev_mapping[node_stop - 1] + 1,
                                                      rev_mapping[node_start]) + start_fend]
            indices = self.data['cis_data'][start_index:stop_index, :]
            counts = indices[:, 2].astype(numpy.float64)
            if corrected:
                counts /= (self.corrections[indices[:, 0]] * self.corrections[indices[:, 1]])
            indices = indices[:, :2]
            indices -= start_fend
            # find bin sums
            _distance.find_distance_bin_sums(mapping,
                                             rev_mapping,
                                             cutoffs,
                                             mids,
                                             counts,
                                             indices,
                                             bin_size,
                                             count_sum,
                                             logdistance_sum,
                                             node_start,
                                             node_stop,
                                             h,
                                             binned)
        if self.rank == 0:
            # exchange arrays
            if not self.silent:
                print >> sys.stderr, ('\r%s\rExchanging distance arrays...') % (' ' * 80),
            for i in range(1, self.num_procs):
                bin_size += self.comm.recv(source=i, tag=11)
                count_sum += self.comm.recv(source=i, tag=11)
                logdistance_sum += self.comm.recv(source=i, tag=11)
            # find chromosome means
            bin_size[-1, :] = numpy.sum(bin_size[:-1, :], axis=0)
            count_sum[-1, :] = numpy.sum(count_sum[:-1, :], axis=0)
            logdistance_sum[-1, :] = numpy.sum(logdistance_sum[:-1, :], axis=0)
            bin_size = bin_size.astype(numpy.float64)
            count_sum = count_sum.astype(numpy.float64)
            valid = numpy.where(count_sum[-1, :] > 0.0)[0]
            count_means = numpy.log(count_sum[-1, valid] / bin_size[-1, valid])
            distance_means = logdistance_sum[-1, valid] / bin_size[-1, valid]
            # find distance line parameters, cutoffs, slopes and intercepts
            distance_parameters = numpy.zeros((valid.shape[0] - 1, 3), dtype=numpy.float32)
            distance_parameters[:-1, 0] = distance_means[1:-1]
            distance_parameters[-1, 0] = numpy.inf
            distance_parameters[:, 1] = ((count_means[1:] - count_means[:-1]) /
                                         (distance_means[1:] - distance_means[:-1]))
            distance_parameters[:, 2] = (count_means[1:] - distance_parameters[:, 1] * distance_means[1:])
            # distribute distance parameters and chromosome means to all nodes
            for i in range(1, self.num_procs):
                self.comm.send(distance_parameters, dest=i, tag=11)
        else:
            self.comm.send(bin_size, dest=0, tag=11)
            self.comm.send(count_sum, dest=0, tag=11)
            self.comm.send(logdistance_sum, dest=0, tag=11)
            distance_parameters = self.comm.recv(source=0, tag=11)
        self.distance_parameters = distance_parameters
        if not corrected:
            chrom_sums = self._find_chromosome_means()
            if self.rank == 0:
                self.chromosome_means = numpy.zeros(chr_indices.shape[0] - 1, dtype=numpy.float32)
                chrom_sums[valid_chroms] /= numpy.sum(bin_size[:-1, :], axis=1)
                chrom_sums[-1] /= numpy.sum(bin_size[-1, :])
                chrom_sums[valid_chroms] = numpy.log(chrom_sums[valid_chroms])
                chrom_sums[-1] = numpy.log(chrom_sums[-1])
                self.chromosome_means[valid_chroms] = (chrom_sums[-1] - chrom_sums[valid_chroms]).astype(numpy.float32)
                for i in range(1, self.num_procs):
                    self.comm.send(self.chromosome_means, dest=i, tag=11)
            else:
                self.chromosome_means = self.comm.recv(source=0, tag=11)
        if not self.silent:
            print >> sys.stderr, ('\r%s\rFinding distance curve... Done\n') % (' ' * 80),
        return None

    def _find_chromosome_means(self):
        chr_indices = self.fends['chr_indices'][...]
        self.chromosome_means = numpy.zeros(chr_indices.shape[0] - 1, dtype=numpy.float32)
        chrom_sums = numpy.zeros(chr_indices.shape[0], dtype=numpy.float64)
        for i in range(chr_indices.shape[0] - 1):
            start_fend = chr_indices[i]
            stop_fend = chr_indices[i + 1]
            num_valid = numpy.sum(self.filter[start_fend:stop_fend])
            if num_valid == 0:
                continue
            rev_mapping = numpy.where(self.filter[start_fend:stop_fend] == 1)[0].astype(numpy.int32)
            if not self.silent:
                print >> sys.stderr, ('\r%s\rFinding mean for chromosome %s...') % \
                                     (' ' * 80, self.fends['chromosomes'][i]),
            mapping = numpy.zeros(stop_fend - start_fend, dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(num_valid, dtype=numpy.int32)
            mids = self.fends['fends']['mid'][rev_mapping + start_fend]
            # pull relevant data
            start_index = self.data['cis_indices'][start_fend]
            stop_index = self.data['cis_indices'][stop_fend]
            node_ranges = numpy.round(numpy.linspace(start_index, stop_index, self.num_procs + 1)).astype(numpy.int64)
            start_index = node_ranges[self.rank]
            stop_index = node_ranges[self.rank + 1]
            data = self.data['cis_data'][start_index:stop_index, :]
            data[:, 0] = mapping[data[:, 0] - start_fend]
            data[:, 1] = mapping[data[:, 1] - start_fend]
            valid = numpy.where((data[:, 0] != -1) * (data[:, 1] != -1))[0]
            data = data[valid, :]
            _distance.find_chromosome_sums(mids,
                                           data,
                                           self.distance_parameters,
                                           chrom_sums,
                                           i)
        if self.rank == 0:
            for i in range(1, self.num_procs):
                chrom_sums += self.comm.recv(source=i, tag=11)
            chrom_sums[-1] = numpy.sum(chrom_sums[:-1])
            for i in range(1, self.num_procs):
                self.comm.Send(chrom_sums, dest=i, tag=13)
        else:
            self.comm.send(chrom_sums, dest=0, tag=11)
            self.comm.Recv(chrom_sums, source=0, tag=13)
        return chrom_sums

    def _smooth_distance_means(self, smoothed):
        """Smooth distance curve in log-space using triangular smoothing of
        smoothed+1 degree."""
        if not self.silent:
            print >> sys.stderr, ("Smoothing distance bin means..."),
        where = numpy.where(self.distance_means > 0)[0]
        log_means = numpy.log(self.distance_means[where])
        new_means = numpy.copy(log_means)
        for i in range(1, smoothed):
            factor = 1
            for j in range(1, i + 1):
                scaling = (i + 1 - j) / float(i + 1)
                new_means[i] += (log_means[i - j] + log_means[i + j]) * scaling
                new_means[-i] += (log_means[-i - j] + log_means[new_means.shape[0] - 1 - i + j]) * scaling
            new_means[i] /= i + 1
            new_means[-i] /= i + 1
        for i in range(1, smoothed + 1):
            new_means[smoothed:(-smoothed)] += ((log_means[(smoothed - i):(-smoothed - i)] +
                                                log_means[(smoothed + i):(new_means.shape[0] - smoothed + i)]) *
                                                (smoothed + 1 - i) / (smoothed + 1))
        new_means[smoothed:(-smoothed)] /= smoothed + 1
        self.distance_means[where] = numpy.exp(new_means)
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None

    def find_fend_corrections(self, mindistance=0, maxdistance=0, minchange=0.0001, burnin_iterations=10000,
                              annealing_iterations=10000, learningrate=0.1, display=0, chroms=None,
                              precalculate=True):
        """
        Using gradient descent, learn correction values for each valid fend based on a Poisson distribution of observations. This function is MPI compatible.

        :param mindistance: The minimum inter-fend distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fend distance to be included in modeling.
        :type maxdistance: int.
        :param minchange: The minimum mean change in fend correction parameter values needed to keep running past 'burnin_iterations' number of iterations during burn-in phase.
        :type minchange: float
        :param burnin_iterations: The number of iterations to use with constant learning rate in gradient descent for learning fend corrections.
        :type burnin_iterations: int.
        :param annealing_iterations: The number of iterations to use with a linearly-decreasing learning rate in gradient descent for learning fend corrections.
        :type annealing_iterations: int.
        :param learningrate: The gradient scaling factor for parameter updates.
        :type learningrate: float
        :param display: Specifies how many iterations between when cost is calculated and displayed as model is learned. If 'display' is zero, the cost is not calculated of displayed.
        :type display: int.
        :param chroms: A list of chromosomes to calculate corrections for. If set as None, all chromosome corrections are found.
        :type chroms: list
        :param precalculate: Specifies whether the correction values should be initialized at the fend means.
        :type precalculate: bool.
        :returns: None
        """
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
            chroms.sort()
        if 'binned' in self.data.attrs.keys() and self.data.attrs['binned'] == True:
            binned = 1
        else:
            binned = 0
        for chrom in chroms:
            chrint = self.chr2int[chrom]
            if self.rank == 0:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rFinding fend correction arrays for chromosome %s...") %\
                        (' ' * 80, chrom),
                start_fend = self.fends['chr_indices'][chrint]
                stop_fend = self.fends['chr_indices'][chrint + 1]
                while start_fend < stop_fend and self.filter[start_fend] == 0:
                    start_fend += 1
                while stop_fend > start_fend and self.filter[stop_fend - 1] == 0:
                    stop_fend -= 1
                if stop_fend > start_fend:
                    for i in range(1, self.num_procs):
                        self.comm.send([start_fend, stop_fend], dest=i, tag=11)
                else:
                    for i in range(1, self.num_procs):
                        self.comm.send([-1, -1], dest=i, tag=11)
                    if not self.silent:
                        print >> sys.stderr, ("insufficient data\n"),
                    continue
            else:
                start_fend, stop_fend = self.comm.recv(source=0, tag=11)
                if start_fend == -1:
                    continue
            num_fends = stop_fend - start_fend
            rev_mapping = numpy.where(self.filter[start_fend:stop_fend] == 1)[0].astype(numpy.int32)
            num_valid = rev_mapping.shape[0]
            node_ranges = numpy.round(numpy.linspace(0, num_valid, self.num_procs + 1)).astype(numpy.int32)
            mapping = numpy.zeros(num_fends, dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(num_valid).astype(numpy.int32)
            fend_ranges = numpy.zeros((num_valid, 3), dtype=numpy.int64)
            mids = self.fends['fends']['mid'][rev_mapping + start_fend]
            if maxdistance == 0:
                maxdistance = mids[-1] - mids[0]
            # find number of downstream interactions for each fend using distance limits
            _distance.find_fend_ranges(rev_mapping,
                                       mids,
                                       fend_ranges,
                                       mindistance,
                                       maxdistance,
                                       node_ranges[self.rank],
                                       node_ranges[self.rank + 1],
                                       binned,
                                       start_fend)
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    fend_ranges[node_ranges[i]:node_ranges[i + 1], :] = self.comm.recv(source=i, tag=11)
                for i in range(1, self.num_procs):
                    self.comm.Send(fend_ranges, dest=i, tag=13)
            else:
                self.comm.send(fend_ranges[node_ranges[self.rank]:node_ranges[self.rank + 1], :], dest=0, tag=11)
                self.comm.Recv(fend_ranges, source=0, tag=13)
            # find range of fends for each node, creating as even spacing as possible
            total_pairs = numpy.sum(fend_ranges[:, 0])
            node_ranges = numpy.round(numpy.linspace(0, total_pairs, self.num_procs + 1)).astype(numpy.int64)
            temp = numpy.r_[0, fend_ranges[:, 0]]
            for i in range(1, temp.shape[0]):
                temp[i] += temp[i - 1]
            start = numpy.searchsorted(temp, node_ranges[self.rank])
            if start > 0:
                if node_ranges[self.rank] - temp[start - 1] < temp[start] - node_ranges[self.rank]:
                    start -= 1
            stop = numpy.searchsorted(temp, node_ranges[self.rank + 1])
            if stop > 0:
                if node_ranges[self.rank + 1] - temp[stop - 1] < temp[stop] - node_ranges[self.rank + 1]:
                    stop -= 1
            num_pairs = temp[stop] - temp[start]
            del temp
            # pull needed data for each node and determine number of nonzero pairs
            start_index = self.data['cis_indices'][start_fend + rev_mapping[start]]
            stop_index = max(start_index, self.data['cis_indices'][start_fend + rev_mapping[stop - 1] + 1])
            temp_data = self.data['cis_data'][start_index:stop_index, :]
            temp_data = temp_data[numpy.where(temp_data[:, 1] < stop_fend)[0], :]
            temp_data[:, 0] = mapping[temp_data[:, 0] - start_fend]
            temp_data[:, 1] = mapping[temp_data[:, 1] - start_fend]
            nonzero_pairs = _distance.find_nonzeros_in_range(fend_ranges, temp_data)
            # allocate nonzero index and count arrays and fill
            nonzero_indices0 = numpy.zeros(nonzero_pairs, dtype=numpy.int32)
            nonzero_indices1 = numpy.zeros(nonzero_pairs, dtype=numpy.int32)
            counts = numpy.zeros(nonzero_pairs, dtype=numpy.int32)
            _distance.find_nonzero_node_indices(fend_ranges,
                                                nonzero_indices0,
                                                nonzero_indices1,
                                                counts,
                                                temp_data)
            del temp_data
            # allocate zero index arrays and fill
            zero_indices0 = numpy.zeros(num_pairs - nonzero_pairs, dtype=numpy.int32)
            zero_indices1 = numpy.zeros(num_pairs - nonzero_pairs, dtype=numpy.int32)
            _distance.find_zero_node_indices(rev_mapping,
                                             fend_ranges,
                                             nonzero_indices0,
                                             nonzero_indices1,
                                             zero_indices0,
                                             zero_indices1,
                                             start,
                                             stop,
                                             binned,
                                             start_fend)
            # count total interactions per fend and sent to root node
            interactions = numpy.bincount(nonzero_indices0, minlength=rev_mapping.shape[0])
            interactions += numpy.bincount(nonzero_indices1, minlength=rev_mapping.shape[0])
            interactions += numpy.bincount(zero_indices0, minlength=rev_mapping.shape[0])
            interactions += numpy.bincount(zero_indices1, minlength=rev_mapping.shape[0])
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    interactions += self.comm.recv(source=i, tag=11)
                temp = numpy.zeros(rev_mapping.shape[0], dtype=numpy.float64)
            else:
                self.comm.send(interactions, dest=0, tag=11)
                temp = None
            # calculate predicted distance means for all interactions
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding distance means for chromosome %s...") % (' ' * 80, chrom),
            nonzero_means = numpy.zeros(nonzero_indices0.shape[0], dtype=numpy.float32)
            _distance.find_remapped_distance_means(nonzero_indices0,
                                                   nonzero_indices1,
                                                   mids,
                                                   nonzero_means,
                                                   self.distance_parameters,
                                                   self.chromosome_means[chrint])
            zero_means = numpy.zeros(zero_indices0.shape[0], dtype=numpy.float32)
            _distance.find_remapped_distance_means(zero_indices0,
                                                   zero_indices1,
                                                   mids,
                                                   zero_means,
                                                   self.distance_parameters,
                                                   self.chromosome_means[chrint])
            del mids
            # if requested, find fend means after adjustment by distance means as correction starting point
            if precalculate:
                enrichments = counts / nonzero_means
                count_sums = numpy.bincount(nonzero_indices0, weights=enrichments, minlength=rev_mapping.shape[0])
                count_sums += numpy.bincount(nonzero_indices1, weights=enrichments, minlength=rev_mapping.shape[0])
                if self.rank == 0:
                    for i in range(1, self.num_procs):
                        count_sums += self.comm.recv(source=i, tag=11)
                    corrections = ((count_sums / interactions.astype(numpy.float32)) ** 0.5).astype(numpy.float32)
                    for i in range(1, self.num_procs):
                        self.comm.send(corrections, dest=i, tag=11)
                else:
                    self.comm.send(count_sums, dest=0, tag=11)
                    corrections = self.comm.recv(source=0, tag=11)
            else:
                corrections = self.corrections[numpy.where(mapping >= 0)[0] + start_fend]
            gradients = numpy.zeros(corrections.shape[0], dtype=numpy.float64)
            cont = True
            # calculate correction gradients
            learningstep = learningrate / max(1, annealing_iterations)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLearning corrections...") % (' ' * 80),
            cost = 0.0
            for phase in ['burnin', 'annealing']:
                if phase == 'burnin':
                    iterations = burnin_iterations
                else:
                    iterations = annealing_iterations
                if iterations == 0:
                    cont = False
                else:
                    cont = True
                iteration = 1
                while cont:
                    gradients.fill(0.0)
                    if display > 0 and iteration%display == 0:
                        findcost = True
                    else:
                        findcost = False
                    _distance.calculate_gradients(zero_indices0,
                                                  zero_indices1,
                                                  nonzero_indices0,
                                                  nonzero_indices1,
                                                  nonzero_means,
                                                  zero_means,
                                                  counts,
                                                  corrections,
                                                  gradients)
                    if findcost:
                        cost = _distance.calculate_cost(zero_indices0,
                                                        zero_indices1,
                                                        nonzero_indices0,
                                                        nonzero_indices1,
                                                        nonzero_means,
                                                        zero_means,
                                                        counts,
                                                        corrections)
                        if self.rank == 0:
                            for i in range(1, self.num_procs):
                                cost += self.comm.recv(source=i, tag=11)
                        else:
                            self.comm.send(cost, dest=0, tag=11)
                    # if using multiple cores, pass gradients to root
                    self._exchange_gradients(gradients, temp)
                    change = self._exchange_corrections(corrections,
                                                         gradients,
                                                         interactions,
                                                         learningrate)
                    if not self.silent and findcost > 0:
                        print >> sys.stderr, ("\r%s phase:%s iteration:%i  cost:%f  change:%f %s") %\
                                             ('Learning corrections...', phase, iteration, cost, change, ' ' * 20),
                    if phase == 'annealing':
                        learningrate -= learningstep
                        if iteration >= annealing_iterations:
                            cont = False
                    elif iteration >= burnin_iterations and change <= minchange:
                        cont = False
                    iteration += 1
            self.corrections[rev_mapping + start_fend] = corrections
            cost = _distance.calculate_cost(zero_indices0,
                                            zero_indices1,
                                            nonzero_indices0,
                                            nonzero_indices1,
                                            nonzero_means,
                                            zero_means,
                                            counts,
                                            corrections)
            if self.rank == 0:
                for i in range(1, self.num_procs):
                    cost += self.comm.recv(source=i, tag=11)
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rLearning corrections... chromosome %s  Final Cost:%f  Done\n") % \
                        (' ' * 80, chrom, cost),
            else:
                self.comm.send(cost, dest=0, tag=11)
        if not self.silent:
            print >> sys.stderr, ("\rLearning corrections... Done%s\n") % (' ' * 80),
        return None

    def _exchange_gradients(self, gradients, temp):
        if self.rank == 0:
            for i in range(1, self.num_procs):
                self.comm.Recv(temp, source=i, tag=13)
                gradients += temp
        else:
            self.comm.Send(gradients, dest=0, tag=13)
        return None

    def _exchange_corrections(self, corrections, gradients, interactions, learningrate):
        if self.rank == 0:
            gradients /= interactions
            old_corrections = numpy.copy(corrections)
            corrections[:] = numpy.minimum(100.0, numpy.maximum(0.01,
                                           corrections - learningrate * gradients))
            change = numpy.amax(numpy.abs(corrections / old_corrections - 1.0))
            for i in range(1, self.num_procs):
                self.comm.Send(corrections, dest=i, tag=13)
            for i in range(1, self.num_procs):
                self.comm.send(change, dest=i, tag=11)
        else:
            self.comm.Recv(corrections, source=0, tag=13)
            change = self.comm.recv(source=0, tag=11)
        return change

    def find_express_fend_corrections(self, iterations=100, mindistance=0, maxdistance=0, remove_distance=True, 
                                      usereads='cis', mininteractions=None, chroms=[]):
        """
        Using iterative approximation, learn correction values for each valid fend.

        :param iterations: The number of iterations to use for learning fend corrections.
        :type iterations: int.
        :param mindistance: This is the minimum distance between fend midpoints needed to be included in the analysis. All possible and observed interactions with a distance shorter than this are ignored. If 'usereads' is set to 'trans', this value is ignored.
        :param maxdistance: The maximum inter-fend distance to be included in modeling. If 'usereads' is set to 'trans', this value is ignored.
        :type maxdistance: int.
        :type mindistance: int.
        :param remove_distance: Specifies whether the estimated distance-dependent portion of the signal is removed prior to learning fend corrections.
        :type remove_distance: bool.
        :param usereads: Specifies which set of interactions to use, 'cis', 'trans', or 'all'.
        :type usereads: str.
        :param mininteractions: If a non-zero 'mindistance' is specified or only 'trans' interactions are used, fend filtering will be performed again to ensure that the data being used is sufficient for analyzed fends. This parameter may specify how many interactions are needed for valid fends. If not given, the value used for the last call to :func:`filter_fends` is used or, barring that, one.
        :type mininteractions: int.
        :param chroms: A list of chromosomes to calculate corrections for. If set as None, all chromosome corrections are found.
        :type chroms: list
        :returns: None
        """
        if self.rank > 0:
            return None
        if mininteractions is None:
            if 'mininteractions' in self.__dict__.keys():
                mininteractions = self.mininteractions
            else:
                mininteractions = 1
        if not self.silent:
            print >> sys.stderr, ("Creating needed arrays for fast correction..."),
        # make sure usereads has a valid value
        read_int = {'cis':0, 'all':1, 'trans':2}
        if usereads not in read_int:
            if not self.silent:
                print >> sys.stderr, ("usereads does not have a valid value.\n"),
            return None
        useread_int = read_int[usereads]
        if 'binned' in self.data.attrs.keys() and self.data.attrs['binned'] == True:
            binned = 1
        else:
            binned = 0
        # create needed arrays
        fend_means = numpy.zeros(self.filter.shape[0], dtype=numpy.float64)
        interactions = numpy.zeros(self.filter.shape[0], dtype=numpy.int64)
        chr_indices = self.fends['chr_indices'][...]
        mids = self.fends['fends']['mid'][...]
        chrs = self.fends['fends']['chr'][...]
        filt = numpy.copy(self.filter)
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
            chroms.sort()
        elif not chroms is None and not isinstance(chroms, list):
            chroms = [chroms]
        for chrm, i in self.chr2int.iteritems():
            if chrm not in chroms:
                filt[chr_indices[i]:chr_indices[i + 1]] = 0
        all_valid = numpy.sum(filt)
        if not self.silent:
            print >> sys.stderr, ("Done\nCopy data for fast correction..."),
        # copy needed arrays from h5dict
        if useread_int < 2:
            data = self.data['cis_data'][...]
            distances = mids[data[:, 1]] - mids[data[:, 0]]
            if maxdistance == 0:
                maxdistance = numpy.amax(distances) + 1
            valid = numpy.where(filt[data[:, 0]] * filt[data[:, 1]] *
                                (distances >= mindistance) * (distances < maxdistance))[0]
            data = data[valid, :]
            data_indices = numpy.zeros(filt.shape[0] + 1, dtype=numpy.int64)
            data_indices[1:] += numpy.bincount(data[:, 0], minlength=filt.shape[0])
            for i in range(1, data_indices.shape[0]):
                data_indices[i] += data_indices[i - 1]
        else:
            data = None
        if useread_int > 0:
            trans_data = self.data['trans_data'][...]
            valid = numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0]
            trans_data = trans_data[valid, :]
        else:
            trans_data = None
        if not self.silent:
            print >> sys.stderr, ("Done\nCount interactions for fast correction..."),
        # double check that, given the type of reads being used for learning, there are enough for each fend
        # to meet the mininteraction criteria
        observed_interactions = numpy.zeros(filt.shape[0], dtype=numpy.int32)
        if not trans_data is None:
            observed_interactions += numpy.bincount(trans_data[:, 0], minlength=filt.shape[0])
            observed_interactions += numpy.bincount(trans_data[:, 1], minlength=filt.shape[0])
        if not data is None:
            observed_interactions += numpy.bincount(data[:, 0], minlength=filt.shape[0])
            observed_interactions += numpy.bincount(data[:, 1], minlength=filt.shape[0])
        if numpy.amin(observed_interactions[numpy.where(filt)]) < mininteractions:
            if not self.silent:
                print >> sys.stderr, ("\nInsufficient interactions for one or more fends.\n"),
                print >> sys.stderr, ("Try resetting and refiltering fends or expanding distance range.\n"),
            return None
        _distance.find_distancebound_interactions(interactions,
                                                  chr_indices,
                                                  mids,
                                                  filt,
                                                  useread_int,
                                                  binned,
                                                  mindistance,
                                                  maxdistance)
        if not self.silent:
            print >> sys.stderr, ("Done\nPrecalculate Distances for fast correction..."),
        # precalculate interaction distance means for all included interactions
        if not remove_distance or data is None:
            distance_means = None
            if trans_data is None:
                mu = (2.0 * numpy.sum(data[:, 2])) / numpy.sum(interactions)
                trans_mu = 1.0
            elif data is None:
                mu = 1.0
                trans_mu = (2.0 * numpy.sum(trans_data[:, 2])) / numpy.sum(interactions)
            else:
                mu = (2.0 * (numpy.sum(data[:, 2]) + numpy.sum(trans_data[:, 2]))) / numpy.sum(interactions)
                trans_mu = mu
        else:
            mu = 1.0
            distance_means = numpy.zeros(data.shape[0], dtype=numpy.float32)
            _distance.find_data_distance_means(distance_means,
                                               filt,
                                               data,
                                               data_indices,
                                               mids,
                                               chrs,
                                               self.distance_parameters,
                                               self.chromosome_means)
            if not trans_data is None:
                total_possible = numpy.sum(filt).astype(numpy.int64) ** 2
                for i in range(self.fends['chr_indices'].shape[0] - 1):
                    start = self.fends['chr_indices'][i]
                    stop = self.fends['chr_indices'][i + 1]
                    total_possible -= numpy.sum(filt[start:stop].astype(numpy.int64)) ** 2
                trans_mu = (2.0 * numpy.sum(trans_data[:, 2])) / total_possible
            else:
                trans_mu = 1.0
        if not self.silent:
            print >> sys.stderr, ("Done\nFinding fend corrections...\n"),
        # calculate corrections
        for iteration in range(iterations):
            cost = _distance.find_fend_means(distance_means,
                                             interactions,
                                             fend_means,
                                             data,
                                             trans_data,
                                             filt,
                                             self.corrections,
                                             mu,
                                             trans_mu)
            if not self.silent:
                print >> sys.stderr, ("\rIteration: %i  Cost: %f    ") % (iteration, cost),
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinal cost: %f\n") % (' ' * 80, cost),
        return None

    def find_trans_means(self):
        """
        Calculate the mean signals across all valid fend-pair trans interactions for each chromosome pair.

        :returns: None
        """
        if not self.silent:
            print >> sys.stderr, ("Finding mean signals across trans interactions..."),
        chr_indices = self.fends['chr_indices'][...]
        num_chroms = chr_indices.shape[0] - 1
        possible = numpy.zeros(num_chroms * (num_chroms - 1) / 2, dtype=numpy.int32)
        pos = 0
        for i in range(chr_indices.shape[0] - 2):
            valid1 = numpy.sum(self.filter[chr_indices[i]:chr_indices[i + 1]])
            for j in range(i + 1, chr_indices.shape[0] - 1):
                valid2 = numpy.sum(self.filter[chr_indices[j]:chr_indices[j + 1]])
                possible[pos] = valid1 * valid2
                pos += 1
        trans_data = self.data['trans_data'][...]
        valid = numpy.where(self.filter[trans_data[:, 0]] * self.filter[trans_data[:, 1]])[0]
        trans_data = trans_data[valid, :]
        del valid
        chrom = self.fends['fends']['chr'][trans_data[:, 0]]
        indices = chrom * (num_chroms - 1) - (chrom * (chrom + 1) / 2) - 1
        indices += self.fends['fends']['chr'][trans_data[:, 1]]
        del chrom
        counts = trans_data[:, 2] / (self.corrections[trans_data[:, 0]] * self.corrections[trans_data[:, 1]])
        del trans_data
        actual = numpy.bincount(indices, weights=counts, minlength=possible.shape[0])
        self.trans_means = actual / numpy.maximum(1.0, possible.astype(numpy.float32))
        if not self.silent:
            print >> sys.stderr, ('Done\n'),
        return None

    def learn_fend_3D_model(self, chrom, minobservations=10):
        """
        Learn coordinates for a 3D model of data using an approximate PCA dimensional reduction.

        This function makes use of the :mod:`mlpy` function :func:`PCAFast` to reduce the data to a set of three coordinates per fend. Cis data for all unfiltered fends for the specified chromosome are dynamically binned to yield a complete distance matrix. The diagonal is set equal to the highest valid enrichment value after dynamic binning. This N x N matrix is passed to :func:`PCAFast` and reduced to an N x 3 matrix.

        :param chrom: The chromosome to learn the model for.
        :type chrom: str.
        :param minobservations: The minimum number of observed reads needed to cease bin expansion in the dynamic binning phase.
        :type minobservations: int.
        :returns: Array containing a row for each valid fend and columns containing X coordinate, Y coordinate, Z coordinate, and sequence coordinate (fend midpoint).
        """
        if 'mlpy' not in sys.modules.keys():
            if not self.silent:
                print >> sys.stderr, ("The mlpy module must be installed to use this function.")
            return None
        if not self.silent:
            print >> sys.stderr, ("Learning fend-resolution 3D model..."),
        unbinned, mapping = hic_binning.unbinned_cis_signal(self, chrom, datatype='fend', arraytype='upper',
                                                            skipfiltered=True, returnmapping=True)
        mids = self.fends['fends']['mid'][mapping]
        if not self.silent:
            print >> sys.stderr, ("Dynamically binning data..."),
        dynamic = numpy.copy(unbinned)
        dynamically_bin_unbinned_upper(unbinned, mids, dynamic, minobservations)
        dynamic = numpy.log(dynamic[:, 0], dynamic[:, 1])
        data_mat = numpy.zeros((mids.shape[0], mids.shape[0]), dtype=numpy.float32)
        indices = numpy.triu_indices(mids.shape[0], 1)
        data_mat[indices] = dynamic
        data_mat[indices[1], indices[0]] = dynamic
        del dynamic
        del indices
        if not self.silent:
            print >> sys.stderr, ("Done\nModeling chromosome..."),
        pca_fast = mlpy.PCAFast(k=3)
        pca_fast.learn(data_mat)
        coordinates = pca_fast.transform(data_mat)
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return numpy.hstack((coordinates, mids.reshape(-1, 1)))

    def cis_heatmap(self, chrom, start=None, stop=None, startfend=None, stopfend=None, binsize=0, binbounds=None,
                    datatype='enrichment', arraytype='compact', maxdistance=0, skipfiltered=False, returnmapping=False,
                    dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                    removefailed=False, image_file=None, **kwargs):
        """
        Return a heatmap of cis data of the type and shape specified by the passed arguments.

        This function returns a heatmap for a single chromosome region, bounded by either 'start' and 'stop' or 'startfend' and 'stopfend' ('start' and 'stop' take precedence), or if given, the outer coordinates of the array passed by 'binbounds'. If none of these are specified, data for the complete chromosome is used. The data in the array is determined by the 'datatype', being raw, fend-corrected, distance-corrected, enrichment, or expected data. The array shape is given by 'arraytype' and can be compact, upper, or full. See :mod:`hic_binning <hifive.hic_binning>` for further explanation of 'datatype' and 'arraytype'. The returned data will include interactions ranging from zero to 'maxdistance' apart. If maxdistance is zero, all interactions within the requested bounds are returned. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param chrom: The name of a chromosome to obtain data from.
        :type chrom: str.
        :param start: The smallest coordinate to include in the array, measured from fend midpoints. If both 'start' and 'startfend' are given, 'start' will override 'startfend'. If unspecified, this will be set to the midpoint of the first fend for 'chrom'. Optional.
        :type start: int.
        :param stop: The largest coordinate to include in the array, measured from fend midpoints. If both 'stop' and 'stopfend' are given, 'stop' will override 'stopfend'. If unspecified, this will be set to the midpoint of the last fend plus one for 'chrom'. Optional.
        :type stop: int.
        :param startfend: The first fend to include in the array. If unspecified and 'start' is not given, this is set to the first fend in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given preference. Optional
        :type startfend: int.
        :param stopfend: The first fend not to include in the array. If unspecified and 'stop' is not given, this is set to the last fend in 'chrom' plus one. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given preference. Optional.
        :type stopfend: str.
        :param binsize: This is the coordinate width of each bin. If 'binsize' is zero, unbinned data is returned. If 'binbounds' is not None, this value is ignored.
        :type binsize: int.
        :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fend not falling in a bin is ignored. Optional.
        :type binbounds: numpy array
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of fends or bins, M is the maximum number of steps between included fend pairs or bin pairs and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
        :type arraytype: str.
        :param maxdistance: This specifies the maximum coordinate distance between bins that will be included in the array. If set to zero, all distances are included.
        :type maxdistance: str.
        :param skipfiltered: If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
        :type skipfiltered: bool.
        :param returnmapping: If 'True', a list containing the data array and a 1d array containing fend numbers included in the data array if unbinned or a 2d array of N x 4 containing the first fend and last fend plus one included in each bin and first and last coordinates if binned is return. Otherwise only the data array is returned.
        :type returnmapping: bool.
        :param dynamically_binned: If 'True', return dynamically binned data.
        :type dynamically_binned: bool.
        :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
        :type minobservations: int.
        :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
        :type searchdistance: int.
        :param expansion_binsize: The size of bins to use for data to pull from when expanding dynamic bins. If set to zero, unbinned data is used.
        :type expansion_binsize: int.
        :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
        :type removefailed: bool.
        :param image_file: If a filename is specified, a PNG image file is written containing the heatmap data. Arguments for the appearance of the image can be passed as additional keyword arguments.
        :type image_file: str.
        :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'. If returnmapping is True, a list is returned containined the requested data array and an array of associated positions (dependent on the binning options selected).

        """
        # check that all values are acceptable
        datatypes = {'raw': 0, 'fend': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
        if datatype not in datatypes:
            if not self.silent:
                print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        else:
            datatype_int = datatypes[datatype]
        if ((dynamically_binned == True and arraytype not in ['compact', 'upper']) or
            (arraytype not in ['full', 'compact', 'upper'])):
            if not self.silent:
                print >> sys.stderr, ("Unrecognized or inappropriate array type. No data returned.\n"),
            return None
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            # determine if data is to be unbinned or binned
            if binbounds is None and binsize == 0:
                # data should be unbinned
                data = hic_binning.unbinned_cis_signal(self, chrom, start=start, stop=stop, startfend=startfend,
                                                       stopfend=stopfend, datatype=datatype, arraytype=arraytype,
                                                       maxdistance=maxdistance, skipfiltered=skipfiltered,
                                                       returnmapping=returnmapping, silent=self.silent)
            else:
                # data should be binned
                data = hic_binning.bin_cis_signal(self, chrom, start=start, stop=stop, startfend=startfend,
                                                  stopfend=stopfend, binsize=binsize, binbounds=binbounds,
                                                  datatype=datatype, arraytype=arraytype, maxdistance=maxdistance,
                                                  returnmapping=returnmapping, silent=self.silent)
        else:
            if expansion_binsize == 0:
                # data should be dynamically binned with unbinned expansion data
                expansion, fends = hic_binning.unbinned_cis_signal(self, chrom, start=start, stop=stop,
                                                                   startfend=startfend, stopfend=stopfend,
                                                                   datatype=datatype, arraytype=arraytype,
                                                                   maxdistance=maxdistance, skipfiltered=True,
                                                                   returnmapping=True, silent=self.silent)
                mids = self.fends['fends']['mid'][fends]
                binned, mapping = hic_binning.bin_cis_array(self, expansion, fends, start=start, stop=stop,
                                                            binsize=binsize, binbounds=binbounds, arraytype=arraytype,
                                                            returnmapping=True, silent=self.silent)
            else:
                # data should be dynamically binned with binned expansion data
                expansion, mapping = hic_binning.bin_cis_signal(self, chrom, start=start, stop=stop,
                                                                startfend=startfend, stopfend=stopfend,
                                                                binsize=expansion_binsize, binbounds=None,
                                                                datatype=datatype, arraytype=arraytype,
                                                                maxdistance=maxdistance, returnmapping=True,
                                                                silent=self.silent)
                mids = (mapping[:, 2] + mapping[:, 3]) / 2
                binned, mapping = hic_binning.bin_cis_signal(self, chrom, start=start, stop=stop, startfend=startfend,
                                                             stopfend=stopfend, binsize=binsize, binbounds=binbounds,
                                                             datatype=datatype, arraytype=arraytype,
                                                             maxdistance=maxdistance, returnmapping=True,
                                                             silent=self.silent)
            hic_binning.dynamically_bin_cis_array(expansion, mids, binned, mapping[:, 2:],
                                                  minobservations=minobservations, searchdistance=searchdistance,
                                                  removefailed=removefailed, silent=self.silent)
            if returnmapping:
                data = [binned, mapping]
            else:
                data = binned
        if not image_file is None:
            if 'symmetricscaling' not in kwargs:
                if datatype == 'enrichment':
                    kwargs['symmetricscaling'] = True
                else:
                    kwargs['symmetricscaling'] = False
            if isinstance(data, list):
                binned = data[0]
            else:
                binned = data
            if arraytype == 'compact':
                img = plotting.plot_compact_array(binned, silent=self.silent, **kwargs)
            elif arraytype == 'full':
                img = plotting.plot_full_array(binned, silent=self.silent, **kwargs)
            else:
                img = plotting.plot_upper_array(binned, silent=self.silent, **kwargs)
            img.save(image_file, format='png')
        return data

    def trans_heatmap(self, chrom1, chrom2, start1=None, stop1=None, startfend1=None, stopfend1=None, binbounds1=None,
                      start2=None, stop2=None, startfend2=None, stopfend2=None, binbounds2=None, binsize=1000000,
                      datatype='enrichment', returnmapping=False, dynamically_binned=False, minobservations=0,
                      searchdistance=0, expansion_binsize=0, removefailed=False, image_file=None, **kwargs):
        """
        Return a heatmap of trans data of the type and shape specified by the passed arguments.

        This function returns a heatmap for trans interactions between two chromosomes within a region, bounded by either 'start1', 'stop1', 'start2' and 'stop2' or 'startfend1', 'stopfend1', 'startfend2', and 'stopfend2' ('start' and 'stop' take precedence), or if given, the outer coordinates of the arrays passed by 'binbounds1' and 'binbounds2'. The data in the array is determined by the 'datatype', being raw, fend-corrected, distance-corrected, enrichment, or expected data. The array shape is always rectangular. See :mod:`hic_binning <hifive.hic_binning>` for further explanation of 'datatype'. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param chrom1: The name of the first chromosome to obtain data from.
        :type chrom1: str.
        :param chrom2: The name of the second chromosome to obtain data from.
        :type chrom2: str.
        :param start1: The coordinate at the beginning of the smallest bin from 'chrom1'. If unspecified, 'start1' will be the first multiple of 'binsize' below the 'startfend1' mid. If there is a conflict between 'start1' and 'startfend1', 'start1' is given preference. Optional.
        :type start1: int.
        :param stop1: The largest coordinate to include in the array from 'chrom1', measured from fend midpoints. If both 'stop1' and 'stopfend1' are given, 'stop1' will override 'stopfend1'. 'stop1' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop1: int.
        :param startfend1: The first fend from 'chrom1' to include in the array. If unspecified and 'start1' is not given, this is set to the first valid fend in 'chrom1'. In cases where 'start1' is specified and conflicts with 'startfend1', 'start1' is given preference. Optional
        :type startfend1: int.
        :param stopfend1: The first fend not to include in the array from 'chrom1'. If unspecified and 'stop1' is not given, this is set to the last valid fend in 'chrom1' + 1. In cases where 'stop1' is specified and conflicts with 'stopfend1', 'stop1' is given preference. Optional.
        :param binbounds1: An array containing start and stop coordinates for a set of user-defined bins to use for partitioning 'chrom1'. Any fend not falling in a bin is ignored.
        :type binbounds1: numpy array
        :param start2: The coordinate at the beginning of the smallest bin from 'chrom2'. If unspecified, 'start2' will be the first multiple of 'binsize' below the 'startfend2' mid. If there is a conflict between 'start2' and 'startfend2', 'start2' is given preference. Optional.
        :type start2: int.
        :param stop2: The largest coordinate to include in the array from 'chrom2', measured from fend midpoints. If both 'stop2' and 'stopfend2' are given, 'stop2' will override 'stopfend2'. 'stop2' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop2: int.
        :param startfend2: The first fend from 'chrom2' to include in the array. If unspecified and 'start2' is not given, this is set to the first valid fend in 'chrom2'. In cases where 'start2' is specified and conflicts with 'startfend2', 'start2' is given preference. Optional
        :type startfend2: int.
        :param stopfend2: The first fend not to include in the array from 'chrom2'. If unspecified and 'stop2' is not given, this is set to the last valid fend in 'chrom2' + 1. In cases where 'stop2' is specified and conflicts with 'stopfend2', 'stop1' is given preference. Optional.
        :type stopfend2: str.
        :param binbounds2: An array containing start and stop coordinates for a set of user-defined bins to use for partitioning 'chrom2'. Any fend not falling in a bin is ignored.
        :type binbounds2: numpy array
        :param binsize: This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
        :type binsize: int.
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param returnmapping: If 'True', a list containing the data array and two 2d arrays of N x 4 containing the first fend and last fend plus one included in each bin and first and last coordinates for the first and second chromosomes is returned. Otherwise only the data array is returned.
        :type returnmapping: bool.
        :param dynamically_binned: If 'True', return dynamically binned data.
        :type dynamically_binned: bool.
        :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
        :type minobservations: int.
        :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
        :type searchdistance: int.
        :param expansion_binsize: The size of bins to use for data to pull from when expanding dynamic bins. If set to zero, unbinned data is used.
        :type expansion_binsize: int.
        :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
        :type removefailed: bool.
        :param image_file: If a filename is specified, a PNG image file is written containing the heatmap data. Arguments for the appearance of the image can be passed as additional keyword arguments.
        :type image_file: str.
        :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'. If returnmapping is True, a list is returned containined the requested data array and an array of associated positions (dependent on the binning options selected).
        """
        # check that all values are acceptable
        datatypes = {'raw': 0, 'fend': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
        if datatype not in datatypes:
            if not self.silent:
                print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        else:
            datatype_int = datatypes[datatype]
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            data = hic_binning.bin_trans_signal(self, chrom1, chrom2, start1=start1, stop1=stop1,
                                                startfend1=startfend1, stopfend1=stopfend1, binbounds1=binbounds1,
                                                start2=start2, stop2=stop2, startfend2=startfend2, stopfend2=stopfend2,
                                                binbounds2=binbounds2, binsize=binsize, datatype=datatype,
                                                returnmapping=returnmapping, silent=self.silent)
        else:
            expansion, mapping1, mapping2 = hic_binning.bin_trans_signal(self, chrom1, chrom2, start1=start1,
                                                                         stop1=stop1, startfend1=startfend1,
                                                                         stopfend1=stopfend1, binbounds1=binbounds1,
                                                                         start2=start2, stop2=stop2,
                                                                         startfend2=startfend2, stopfend2=stopfend2,
                                                                         binbounds2=binbounds2,
                                                                         binsize=expansion_binsize, datatype=datatype,
                                                                         returnmapping=True, silent=self.silent)
            mids1 = (mapping1[:, 2] + mapping1[:, 3]) / 2
            mids2 = (mapping2[:, 2] + mapping2[:, 3]) / 2
            binned, mapping1, mapping2 = hic_binning.bin_trans_signal(self, chrom1, chrom2, start1=start1, stop1=stop1,
                                                                      startfend1=startfend1, stopfend1=stopfend1,
                                                                      binbounds1=binbounds1, start2=start2,
                                                                      stop2=stop2, startfend2=startfend2,
                                                                      stopfend2=stopfend2, binbounds2=binbounds2,
                                                                      binsize=binsize, datatype=datatype,
                                                                      returnmapping=True, silent=self.silent)
            hic_binning.dynamically_bin_trans_array(expansion, mids1, mids2, binned, mapping1[:, 2:], mapping2[:, 2:],
                                                    minobservations=minobservations, searchdistance=searchdistance,
                                                    removefailed=removefailed, silent=self.silent)
            if returnmapping:
                data = [binned, mapping1, mapping2]
            else:
                data = binned
        if not image_file is None:
            if 'symmetricscaling' not in kwargs:
                if datatype == 'enrichment':
                    kwargs['symmetricscaling'] = True
                else:
                    kwargs['symmetricscaling'] = False
            if isinstance(data, list):
                binned = data[0]
            else:
                binned = data
            img = plotting.plot_full_array(binned, silent=self.silent, **kwargs)
            img.save(image_file, format='png')
        return data

    def write_heatmap_dict(self, filename, binsize, includetrans=True, datatype='enrichment', chroms=[]):
        """
        Create an h5dict file containing binned interaction arrays, bin positions, and an index of included chromosomes. This function is MPI compatible.

        :param filename: Location to write h5dict object to.
        :type filename: str.
        :param binsize: Size of bins for interaction arrays.
        :type binsize: int.
        :param includetrans: Indicates whether trans interaction arrays should be calculated and saved.
        :type includetrans: bool.
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param chroms: A list of chromosome names indicating which chromosomes should be included. If left empty, all chromosomes are included. Optional.
        :type chroms: list
        :returns: None
        """
        if (chroms is None or
                (isinstance(chroms, list) and
                (len(chroms) == 0 or
                (len(chroms) == 1 and chroms[0] == ''))) or
                chroms == ''):
            chroms = self.chr2int.keys()
            chroms.sort()
        hic_binning.write_heatmap_dict(self, filename, binsize, includetrans=includetrans,
                                       datatype=datatype, chroms=chroms, silent=self.silent)
        return None
