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
    :returns: :class:`HiC <hifive.hic.HiC>` class object.
    """

    def __init__(self, filename, mode='r'):
        """
        Create a HiC object.
        """

        self.file = os.path.abspath(filename)
        if mode != 'w':
            self.load()
        if 'mpi4py' in sys.modules.keys():
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.num_procs = self.comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.num_procs = 1
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
        filename = os.path.abspath(filename)
        # ensure data h5dict exists
        if not os.path.exists(filename):
            print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename),
            return None
        self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                       os.path.dirname(self.file)), os.path.basename(filename))
        self.data = HiCData(filename).data
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
            print >> sys.stderr, ("Could not find %s.\n") % (fendfilename),
            return None
        self.fends = Fend(fendfilename).fends
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.fends['chromosomes']):
            self.chr2int[chrom] = i
        # create arrays
        self.filter = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.int32)
        self.corrections = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.float32)
        return None

    def save(self):
        """
        Save analysis parameters to h5dict.

        :returns: None
        """
        if self.rank > 0:
            return None
        datafile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['data', 'fends', 'file', 'chr2int', 'comm', 'rank', 'num_procs']:
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
                print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (datafilename),
            else:
                self.data = HiCData(datafilename).data
        # ensure fend h5dict exists
        if 'fendfilename' in self.__dict__:
            fendfilename = self.fendfilename
            if fendfilename[:2] == './':
                fendfilename = fendfilename[2:]
            parent_count = fendfilename.count('../')
            fendfilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(fendfilename):
                print >> sys.stderr, ("Could not find %s. No fends loaded.\n") % (fendfilename),
            else:
                self.fends = Fend(fendfilename).fends
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.fends['chromosomes']):
            self.chr2int[chrom] = i
        datafile.close()
        return None

    def filter_fends(self, mininteractions=10, maxdistance=0):
        """
        Iterate over the dataset and remove fends that do not have 'minobservations' within 'maxdistance' of themselves using only unfiltered fends.

        In order to create a set of fends that all have the necessary number of interactions, after each round of filtering, fend interactions are retallied using only interactions that have unfiltered fends at both ends.

        :param mininteractions: The required number of interactions for keeping a fend in analysis.
        :type mininteractions: int.
        :param maxdistance: The maximum inter-fend distance used to count fend interactions. A value of 0 indicates all cis-data should be used.
        :type maxdistance: int.
        :returns: None
        """
        print >> sys.stderr, ("Filtering fends..."),
        self.mininteractions = mininteractions
        self.maxdistance = maxdistance
        original_count = numpy.sum(self.filter)
        previous_valid = original_count + 1
        current_valid = original_count
        coverage = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        # determine maximum ranges of valid interactions for each fend
        max_fend = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        _distance.find_max_fend(max_fend,
                                self.fends['fends']['mid'][...],
                                self.fends['fends']['chr'][...],
                                self.fends['chr_indices'][...],
                                0,
                                maxdistance)
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
                               max_fend,
                               coverage,
                               mininteractions)
            current_valid = numpy.sum(self.filter)
        print >> sys.stderr, ("Removed %i of %i fends\n") % (original_count - current_valid, original_count),
        return None

    def find_distance_means(self, numbins=90, minsize=200, maxsize=0, smoothed=0, startfend=0, stopfend=None,
                             corrected=False):
        """
        Count reads and possible interactions from valid fend pairs in each distance bin to find mean bin signals. This function is MPI compatible.

        This partitions the range of interaction distances (measured from mipoints of the involved fends) from the 'minsize' to 'maxsize' into a number of partitions equal to 'numbins'. The first bin contains all distances less than or equal to 'minsize'. The remaining bins are defined such that their log ranges are equal to one another. The curve defined by the mean interaction value of each bin can be smoothed using a  triangular smoothing operation.

        :param numbins: The number of bins to divide the distance range into. The first bin extends from zero to 'minsize', while the remaining bins are divided into evenly-spaced log-sized bins from 'minsize' to 'maxsize' or the maximum inter-fend distance, whichever is greater.
        :type numbins: int.
        :param minsize: The upper size limit of the smallest distance bin.
        :type minsize: int.
        :param maxsize: If this value is larger than the largest included chromosome, it will extend bins out to maxsize. If this value is smaller, it is ignored.
        :type maxsize: int.
        :param smoothed: Indicates the degree of smoothing to the distance curve for noise reduction. A value of less than 1 indicates no smoothing.
        :type smoothed: int.
        :param startfend: The first fend to include interactions from in calculating bin means.
        :type startfend: int.
        :param stopfend: The first fend to exclude interactions from in calculating bin means. If stopfend is None, there is no upper bound on included fends.
        :type stopfend: int.
        :param corrected: If True, correction values are applied to counts prior to summing.
        :type corrected: bool.
        :returns: None
        """
        if stopfend is None:
            stopfend = self.filter.shape[0]
        self.distance_smoothing = smoothed
        if self.rank == 0:
            print >> sys.stderr, ("Find distance means..."),
            # determine maximum valid inter-fend distance
            max_distance = 0
            for i in range(self.fends['fends']['chr'][startfend],
                           self.fends['fends']['chr'][stopfend - 1] + 1):
                start_fend = max(self.fends['chr_indices'][i], startfend)
                while start_fend < min(stopfend, self.fends['chr_indices'][i + 1]) and self.filter[start_fend] == 0:
                    start_fend += 1
                stop_fend = min(stopfend, self.fends['chr_indices'][i + 1])
                while stop_fend > max(start_fend, self.fends['chr_indices'][i]) and self.filter[stop_fend - 1] == 0:
                    stop_fend -= 1
                if start_fend < stop_fend:
                    max_distance = max(max_distance, self.fends['fends']['mid'][stop_fend - 1] -
                                       self.fends['fends']['mid'][start_fend])
            max_distance = max(maxsize, max_distance)
            # find bin cutoffs
            self.distance_bins = numpy.empty(numbins, dtype=numpy.int32)
            self.distance_bins[0] = minsize
            if numbins > 1:
                binsize = (log(max_distance) - log(minsize)) / (numbins - 1)
                self.distance_bins[1:] = numpy.ceil(numpy.exp(numpy.arange(1, numbins) * binsize +
                                                    log(minsize))).astype(numpy.int32)
            self.distance_mids = numpy.zeros(numbins, dtype=numpy.float32)
            self.distance_mids[:] = self.distance_bins.astype(numpy.float32) / exp(binsize * (1.0 - 0.5 ** 0.5))
            # send bins to other workers
            for i in range(1, self.num_procs):
                self.comm.send([max_distance, self.distance_bins, self.distance_mids], dest=i, tag=11)
            # divide fends amongst workers
            needed = (numpy.where(self.filter[startfend:stopfend] > 0)[0] + startfend).astype(numpy.int32)
            numpy.random.shuffle(needed)
            worker_size = int(ceil(needed.shape[0] / float(self.num_procs)))
            for i in range(1, self.num_procs):
                self.comm.send(needed[((i - 1) * worker_size):(i * worker_size)], dest=i, tag=11)
            node_needed = needed[((self.num_procs - 1) * worker_size):]
        else:
            max_distance, self.distance_bins, self.distance_mids = self.comm.recv(source=0, tag=11)
            node_needed = self.comm.recv(source=0, tag=11)
        # allocate arrays
        bin_sums = numpy.zeros(self.distance_bins.shape[0], dtype=numpy.float64)
        bin_counts = numpy.zeros(self.distance_bins.shape[0], dtype=numpy.int64)
        data = self.data['cis_data'][...]
        data_indices = self.data['cis_indices'][...]
        mids = self.fends['fends']['mid'][...]
        chroms = self.fends['fends']['chr'][...]
        chr_indices = self.fends['chr_indices']
        # update arrays for each fend
        for fend in node_needed:
            stop_fend = min(stopfend, chr_indices[chroms[fend] + 1])
            _distance.find_fend_distance_bin_values(self.filter,
                                                    self.corrections,
                                                    data,
                                                    data_indices,
                                                    mids,
                                                    self.distance_bins,
                                                    bin_sums,
                                                    bin_counts,
                                                    fend,
                                                    stop_fend,
                                                    int(corrected))
        if self.rank == 0:
            for i in range(1, self.num_procs):
                temp = self.comm.recv(source=i, tag=11)
                bin_sums += temp[0]
                bin_counts += temp[1]
            where = numpy.where(bin_counts > 0)[0]
            bin_sums[where] /= bin_counts[where]
            self.distance_means = bin_sums.astype(numpy.float32)[:]
            print >> sys.stderr, ("Done\n"),
            # if requested, smooth distance curve
            if smoothed > 0:
                self._smooth_distance_means(smoothed)
            # send results back to workers
            for i in range(1, self.num_procs):
                self.comm.send(self.distance_means, dest=i, tag=11)
        else:
            self.comm.send([bin_sums, bin_counts], dest=0, tag=11)
            self.distance_means = self.comm.recv(source=0, tag=11)
        return None

    def _smooth_distance_means(self, smoothed):
        """Smooth distance curve in log-space using triangular smoothing of
        smoothed+1 degree."""
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
        print >> sys.stderr, ("Done\n"),
        return None

    def find_fend_corrections(self, maxdistance=0, burnin_iterations=10000, annealing_iterations=10000,
                              learningrate=0.01, recalculate_distance=0, display=0):
        """
        Using gradient descent, learn correction values for each valid fend based on a Poisson distribution of observations. This function is MPI compatible.

        :param maxdistance: The maximum inter-fend distance to be included in modeling.
        :type maxdistance: int.
        :param burnin_iterations: The number of iterations to use with constant learning rate in gradient descent for learning fend corrections.
        :type burnin_iterations: int.
        :param annealing_iterations: The number of iterations to use with a linearly-decreasing learning rate in gradient descent for learning fend corrections.
        :type annealing_iterations: int.
        :param learningrate: The gradient scaling factor for parameter updates.
        :type learningrate: float
        :param recalculate_distance: Number of iterations that should pass before recalculating the distance bin means to account for the updated fend corrections. If set to zero, no recalculation is performed.
        :type recalculate_distance: int.
        :param display: Specifies how many iterations between when cost is calculated and displayed as model is learned. If 'display' is zero, the cost is not calculated of displayed.
        :type display: int.
        :returns: None
        """
        num_bins = self.distance_bins.shape[0]
        minsize = self.distance_bins[0]
        smoothed = int(self.distance_smoothing)
        # find range of fends for worker. starts and stops must always include whole fragments
        valid = numpy.where(self.filter)[0]
        worker_size = int(ceil(valid.shape[0] / float(self.num_procs)))
        if self.rank == 0:
            print >> sys.stderr, ("Find fend correction arrays..."),
            if self.num_procs > 1:
                start_fend = valid[worker_size * (self.num_procs - 1)]
                start_fend -= start_fend % 2
                stop_fend = valid[-1] + 1
                stop_fend += stop_fend % 2
            else:
                start_fend = 0
                stop_fend = self.filter.shape[0]
        else:
            start_fend = valid[worker_size * (self.rank - 1)]
            start_fend -= start_fend % 2
            stop_fend = valid[worker_size * self.rank]
            stop_fend -= stop_fend % 2
        # create arrays needed for root and find possible interactions
        if self.rank == 0:
            all_gradients = numpy.zeros(self.filter.shape[0], dtype=numpy.float32)
            all_corrections = numpy.copy(self.corrections).astype(numpy.float32)
            interactions = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
            _distance.find_possible_interactions(interactions,
                                                 self.filter,
                                                 self.fends['fends']['mid'][:],
                                                 self.fends['chr_indices'][...],
                                                 maxdistance)
            fend_ranges = numpy.zeros((self.num_procs, 3), dtype=numpy.int32)
        else:
            all_gradients = None
            all_corrections = None
            fend_ranges = None
        # create arrays common to all workers
        num_fends = stop_fend - start_fend
        max_fend = numpy.zeros(num_fends, dtype=numpy.int32)
        _distance.find_max_fend(max_fend,
                                self.fends['fends']['mid'][start_fend:],
                                self.fends['fends']['chr'][start_fend:stop_fend],
                                self.fends['chr_indices'][...],
                                start_fend,
                                maxdistance)
        num_fends2 = max_fend[num_fends - 1]
        stop_fend2 = num_fends2 + start_fend
        gradients = numpy.zeros(num_fends2, dtype=numpy.float32)
        corrections = self.corrections[start_fend:stop_fend2].astype(numpy.float32)
        # if using multiple cores, send worker ranges to root and allocate arrays for receiving gradient data
        if self.num_procs > 1:
            if self.rank != 0:
                self.comm.send([start_fend, stop_fend, stop_fend2], dest=0, tag=11)
                gradient_arrays = None
            else:
                gradient_arrays = {}
                fend_ranges[0, 0] = start_fend
                fend_ranges[0, 1] = stop_fend
                fend_ranges[0, 2] = stop_fend2
                for i in range(1, self.num_procs):
                    temp = self.comm.recv(source=i, tag=11)
                    fend_ranges[i, 0] = temp[0]
                    fend_ranges[i, 1] = temp[1]
                    fend_ranges[i, 2] = temp[2]
                    gradient_arrays[i] = numpy.empty(temp[2] - temp[0], dtype=numpy.float32)
        # pull only needed data, adjusting as necessary
        data_indices = self.data['cis_indices'][start_fend:(stop_fend + 1)]
        data = self.data['cis_data'][data_indices[0]:data_indices[-1], :]
        data_indices -= data_indices[0]
        data[:, :2] -= start_fend
        filter = self.filter[start_fend:stop_fend2]
        original_chr_indices = self.fends['chr_indices'][...]
        chr_indices = original_chr_indices - start_fend
        # precalculate interaction distance means for all included interactions, if needed
        max_bin = numpy.amax(max_fend - numpy.arange(num_fends))
        interaction_means = numpy.zeros((num_fends, max_bin), dtype=numpy.float32)
        distance_mid_logs = numpy.log(self.distance_mids).astype(numpy.float32)
        distance_mean_logs = numpy.log(self.distance_means).astype(numpy.float32)
        _distance.find_interaction_distance_means(interaction_means,
                                                  filter,
                                                  self.fends['fends']['mid'][start_fend:stop_fend2],
                                                  self.distance_means,
                                                  self.distance_mean_logs,
                                                  self.distance_mids,
                                                  distance_mid_logs,
                                                  distance_mid_logs[1:] - distance_mid_logs[:-1],
                                                  max_fend)
      # calculate correction gradients
        if self.rank == 0:
            print >> sys.stderr, ("Done\n"),
            print >> sys.stderr, ("Learning corrections..."),
        learningstep = learningrate / max(1, annealing_iterations)
        for phase in ['burnin', 'annealing']:
            if phase == 'burnin':
                iterations = burnin_iterations
            else:
                iterations = annealing_iterations
            for iteration in range(iterations):
                gradients.fill(0.0)
                if display > 0 and iteration%display == 0:
                    findcost = 1
                else:
                    findcost = 0
                cost = _distance.calculate_gradients(data,
                                                     data_indices,
                                                     filter,
                                                     interaction_means,
                                                     corrections,
                                                     gradients,
                                                     max_fend,
                                                     findcost)
                # if using multiple cores, pass gradients to root
                if self.num_procs > 1:
                    cost = self._exchange_gradients(gradients, all_gradients, fend_ranges, gradient_arrays, cost)
                else:
                    all_gradients = gradients
                if self.rank == 0:
                    if findcost > 0:
                        print >> sys.stderr, ("\r%s phase:%s iteration:%i  cost:%f ") %\
                                             ('Learning corrections...', phase, iteration, cost),
                    # update corrections on root only
                    _distance.update_corrections(all_corrections,
                                                 all_gradients,
                                                 interactions,
                                                 learningrate)
                # if using multiple cores, distribute needed correction values to workers
                if self.num_procs > 1:
                    self._exchange_corrections(corrections, all_corrections, fend_ranges)
                else:
                    corrections = all_corrections[start_fend:stop_fend2]
                if phase == 'annealing':
                    learningrate -= learningstep
                if recalculate_distance > 0 and (iteration + 1) % recalculate_distance == 0:
                    if self.rank == 0:
                        self.corrections = all_corrections.astype(numpy.float32)
                        if self.num_procs > 1:
                            for i in range(1, self.num_procs):
                                self.comm.Send(self.corrections, dest=i, tag=13)
                    else:
                        self.comm.Recv(self.corrections, source=0, tag=13)
                        self.corrections = self.corrections.astype(numpy.float32)
                    self.find_distance_means(numbins=num_bins, minsize=minsize, smoothed=smoothed,
                                             corrected=True)
                    interaction_means[:, :] = 0.0
                    distance_mid_logs = numpy.log(self.distance_mids).astype(numpy.float32)
                    distance_mean_logs = numpy.log(self.distance_means).astype(numpy.float32)
                    _distance.find_interaction_distance_means(interaction_means,
                                                  filter,
                                                  self.fends['fends']['mid'][start_fend:stop_fend2],
                                                  self.distance_means,
                                                  self.distance_mean_logs,
                                                  self.distance_mids,
                                                  distance_mid_logs,
                                                  distance_mid_logs[1:] - distance_mid_logs[:-1],
                                                  max_fend)
        if self.rank == 0:
            self.corrections = all_corrections.astype(numpy.float32)
            if self.num_procs > 1:
                for i in range(1, self.num_procs):
                    self.comm.Send(self.corrections, dest=i, tag=13)
            print >> sys.stderr, ("\rLearning corrections... Done%s\n") % (' ' * 60),
        elif self.num_procs > 1:
            self.comm.Recv(self.corrections, source=0, tag=13)
            self.corrections = self.corrections.astype(numpy.float32)
        return None

    def find_fend_corrections2(self, mindistance=0, maxdistance=0, minchange=0.0, burnin_iterations=10000,
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
        if chroms is None:
            chroms = self.chr2int.keys()
            chroms.sort()
        for chrom in chroms:
            chrint = self.chr2int[chrom]
            if self.rank == 0:
                print >> sys.stderr, ("%s\rFinding fend correction arrays for chromosome %s...") % (' ' * 80, chrom),
                start_fend = self.chr_indices[chrint]
                stop_fend = self.chr_indices[chrint + 1]
                while start_fend < stop_fend and self.filter[start_fend] == 0:
                    start_fend += 1
                while stop_fend > start_fend and self.filter[stop_fend - 1] == 0:
                    stop_fend -= 1
                num_valid = numpy.sum(self.filter[start_fend:stop_fend])
                indices0, indices1 = numpy.triu_indices(num_valid, 1)
                if self.num_procs > 1:
                    node_ranges = numpy.round(numpy.linspace(0, indices.shape[0], self.num_procs)).astype(numpy.int32)
                    for i in range(1, self.num_procs):
                        self.comm.send([start_fend, stop_fend], dest=i, tag=11)
                        self.comm.send(indices0[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                        self.comm.send(indices1[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    indices0 = indices0[:node_ranges[1]]
                    indices1 = indices1[:node_ranges[1]]
            else:
                start_fend, stop_fend = self.comm.recv(source=0, tag=11)
                indices0 = self.comm.recv(source=0, tag=11)
                indices1 = self.comm.recv(source=0, tag=11)
                num_valid = numpy.sum(self.filter[start_fend:stop_fend])
            num_fends = stop_fend - start_fend
            mapping = numpy.zeros(num_fends, dtype=numpy.int32) - 1
            rev_mapping = numpy.where(self.filter[start_fend:stop_fend] == 1)[0].astype(numpy.int32)
            mapping[rev_mapping] = numpy.arange(num_valid)
            distances = numpy.zeros(indices0.shape[0], dtype=numpy.float32)
            distance_mid_logs = numpy.log(self.distance_mids).astype(numpy.float32)
            distance_mean_logs = numpy.log(self.distance_means).astype(numpy.float32)
            _distance.find_remapped_distance_means(indices0, indices1, distances,
                                                   self.fends['fends']['mid'][rev_mapping + start_fend],
                                                   self.distance_mids, self.distance_means, distance_mean_logs,
                                                   distance_mid_logs, distance_mid_logs[1:] - distance_mid_logs[:-1],
                                                   mindistance, maxdistance)
            valid = numpy.where(distances > 0.0)[0]
            indices0 = indices0[valid]
            indices1 = indices1[valid]
            distances = distances[valid]
            if self.rank == 0:
                if self.num_procs > 1:
                    for i in range(1, self.num_procs):
                        indices0 = numpy.hstack((indices0, self.comm.recv(source=i, tag=11)))
                    node_ranges = numpy.round(numpy.linspace(0, indices0.shape[0], self.num_procs)).astype(numpy.int32)
                    for i in range(1, self.num_procs):
                        self.comm.send(indices0[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    indices0 = indices0[:node_ranges[1]]
                    for i in range(1, self.num_procs):
                        indices1 = numpy.hstack((indices1, self.comm.recv(source=i, tag=11)))
                    for i in range(1, self.num_procs):
                        self.comm.send(indices1[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    indices1 = indices1[:node_ranges[1]]
                    for i in range(1, self.num_procs):
                        distances = numpy.hstack((distances, self.comm.recv(source=i, tag=11)))
                    for i in range(1, self.num_procs):
                        self.comm.send(distances[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    distances = distances[:node_ranges[1]]
                interactions = numpy.bincount(indices0, minlength=rev_mapping.shape[0])
                interactions += numpy.bincount(indices1, minlength=rev_mapping.shape[0])
                if self.num_procs > 1:
                    for i in range(1, self.num_procs):
                        interactions += self.comm.recv(source=i, tag=11)
                    temp = numpy.zeros(rev_mapping.shape[0], dtype=numpy.float32)
            else:
                self.comm.send(indices0, dest=0, tag=11)
                indices0 = self.comm.recv(source=0, tag=11)
                self.comm.send(indices1, dest=0, tag=11)
                indices1 = self.comm.recv(source=0, tag=11)
                self.comm.send(distances, dest=0, tag=11)
                distances = self.comm.recv(source=0, tag=11)
                interactions = numpy.bincount(indices0, minlength=rev_mapping.shape[0])
                interactions += numpy.bincount(indices1, minlength=rev_mapping.shape[0])
                self.comm.send(interactions, dest=0, tag=11)
                temp = None
            data = numpy.zeros(indices0.shape[0], dtype=numpy.int32)
            start_index = self.data['cis_indices'][start_fend + indices0[0]]
            stop_index = self.data['cis_indices'][start_fend + indices0[-1] + 1]
            temp_data = self.data['cis_data'][start_index:stop_index, :]
            remap_counts(indices0, indices1, mapping, data, temp_data, start_fend)
            if precalculate:
                count_sums = numpy.bincount(indices0, weights=data, minlength=rev_mapping.shape[0])
                count_sums += numpy.bincount(indices1, weights=data, minlength=rev_mapping.shape[0])
                if self.rank == 0:
                    for i in range(1, self.num_procs):
                        count_sums += self.comm.recv(source=i, tag=11)
                    corrections = (count_sums / interactions.astype(numpy.float32)).astype(numpy.float32)
                    for i in range(1, self.num_procs):
                        self.comm.send(corrections, dest=i, tag=11)
            else:
                corrections = self.corrections[numpy.where(mapping >= 0)[0] + start_fend]
            prev_corrections = numpy.copy(corrections)
            gradients = numpy.zeros(corrections.shape[0], dtype=numpy.float32)
            cont = True
            node_ranges = numpy.round(numpy.linspace(0, corrections.shape[0], self.num_procs)).astype(numpy.int32)
            start_index = node_ranges[self.rank]
            stop_index = node_ranges[self.rank + 1]
            # calculate correction gradients
            learningstep = learningrate / max(1, annealing_iterations)
            if self.rank == 0:
                print >> sys.stderr, ("\r%s") % (' ' * 80),
            for phase in ['burnin', 'annealing']:
                if phase == 'burnin':
                    iterations = burnin_iterations
                else:
                    iterations = annealing_iterations
                iteration = 1
                cont = True
                while cont:
                    gradients.fill(0.0)
                    if display > 0 and iteration%display == 0:
                        findcost = 1
                    else:
                        findcost = 0
                    cost = _distance.calculate_gradients2(indices0,
                                                          indices1,
                                                          data,
                                                          distances,
                                                          corrections,
                                                          gradients,
                                                          findcost)
                    # if using multiple cores, pass gradients to root
                    if self.num_procs > 1:
                        cost = self._exchange_gradients2(gradients, temp, cost)
                    # update corrections
                    change = self._exchange_corrections2(corrections,
                                                         previous_corrections,
                                                         gradients,
                                                         interactions,
                                                         learningrate,
                                                         start_index,
                                                         stop_index)
                    previous_corrections = numpy.copy(corrections)
                    if self.rank == 0 and findcost > 0:
                        print >> sys.stderr, ("\r%s phase:%s iteration:%i  cost:%f  change:%f ") %\
                                             ('Learning corrections...', phase, iteration, cost, change),
                    if phase == 'annealing':
                        learningrate -= learningstep
                        if iteration == annealing_iterations:
                            cont = False
                    else:
                        if iteration >= burnin_iterations and change <= minchange:
                            cont = False
                    iteration += 1
            self.corrections[rev_mapping + start_fend] = corrections
        if self.rank == 0:
            print >> sys.stderr, ("\rLearning corrections... Done%s\n") % (' ' * 60),
        return None

    def _exchange_gradients(self, gradients, all_gradients, fend_ranges, gradient_arrays, cost):
        """Compile all gradients on root"""
        if self.rank == 0:
            all_gradients.fill(0.0)
            all_gradients[fend_ranges[0, 0]:fend_ranges[0, 2]] += gradients
            for i in range(1, fend_ranges.shape[0]):
                cost += self.comm.recv(source=i, tag=11)
                self.comm.Recv(gradient_arrays[i], source=i, tag=13)
                all_gradients[fend_ranges[i, 0]:fend_ranges[i, 2]] += gradient_arrays[i]
        else:
            self.comm.send(cost, dest=0, tag=11)
            self.comm.Send(gradients, dest=0, tag=13)
        return cost

    def _exchance_gradients2(self, gradients, temp, cost):
        if self.rank == 0:
            for i in range(1, self.num_procs):
                self.comm.Recv(temp, source=i, tag=13)
                gradients += temp
                cost += self.comm.recv(source=i, tag=11)
        else:
            self.comm.Send(gradients, dest=0, tag=13)
            self.comm.send(cost, dest=0, tag=11)
        return None

    def _exchange_corrections(self, corrections, all_corrections, fend_ranges):
        """Distribute correction values needed to each worker"""
        if self.rank == 0:
            for i in range(1, fend_ranges.shape[0]):
                self.comm.Send(all_corrections[fend_ranges[i, 0]:fend_ranges[i, 2]], dest=i, tag=13)
            corrections[:] = all_corrections[fend_ranges[0, 0]:fend_ranges[0, 2]]
        else:
            self.comm.Recv(corrections, source=0, tag=13)
        return None

    def _exchange_corrections2(self, corrections, previous_corrections, gradients, interactions, learningrate,
                               start_index, stop_index):
        if self.rank == 0:
            corrections -= learningrate * gradients
            for i in range(1, self.num_procs):
                self.comm.Send(corrections, dest=i, tag=13)
            change = numpy.sum(numpy.abs(previous_corrections[start_index:stop_index] /
                                         corrections[start_index:stop_index] - 1.0))
            for i in range(1, self.num_procs):
                change += self.comm.recv(source=i, tag=11)
            change /= corrections.shape[0]
            for i in range(1, self.num_procs):
                self.comm.send(change, dest=i, tag=11)
        else:
            self.comm.Recv(corrections, source=i, tag=13)
            change = numpy.sum(numpy.abs(previous_corrections[start_index:stop_index] /
                                         corrections[start_index:stop_index] - 1.0))
            self.comm.send(change, dest=0, tag=11)
            change = self.comm.recv(source=0, tag=11)
        return change

    def find_express_fend_corrections(self, iterations=100, mindistance=0, remove_distance=True, usereads='cis',
                                      recalculate_distance=0, mininteractions=None):
        """
        Using iterative approximation, learn correction values for each valid fend. This function is MPI compatible.

        :param iterations: The number of iterations to use for learning fend corrections.
        :type iterations: int.
        :param mindistance: This is the minimum distance between fend midpoints needed to be included in the analysis. All possible and observed interactions with a distance shorter than this are ignored. If 'usereads' is set to 'trans', this value is ignored.
        :type mindistance: int.
        :param remove_distance: Specifies whether the estimated distance-dependent portion of the signal is removed prior to learning fend corrections.
        :type remove_distance: bool.
        :param usereads: Specifies which set of interactions to use, 'cis', 'trans', or 'all'.
        :type usereads: str.
        :param recalculate_distance: Number of iterations that should pass before recalculating the distance bin means to account for the updated fend corrections. If set to zero, no recalculation is performed.
        :type recalculate_distance: int.
        :param mininteractions: If a non-zero 'mindistance' is specified or only 'trans' interactions are used, fend filtering will be performed again to ensure that the data being used is sufficient for analyzed fends. This parameter may specify how many interactions are needed for valid fends. If not given, the value used for the last call to :func:`filter_fends` is used or, barring that, one.
        :type mininteractions: int.
        :returns: None
        """
        num_bins = self.distance_bins.shape[0]
        minsize = self.distance_bins[0]
        smoothed = int(self.distance_smoothing)
        if mininteractions is None:
            if 'mininteractions' in self.__dict__.keys():
                mininteractions = self.mininteractions
            else:
                mininteractions = 1
        if self.rank == 0:
            print >> sys.stderr, ("Creating needed arrays for fast correction..."),
            # make sure usereads has a valid value
            read_int = {'cis':0, 'all':1, 'trans':2}
            if usereads not in read_int:
                print >> sys.stderr, ("usereads does not have a valid value.\n"),
                return None
            useread_int = read_int[usereads]
            # create needed arrays
            fend_means = numpy.zeros(self.filter.shape[0], dtype=numpy.float64)
            interactions = numpy.zeros(self.filter.shape[0], dtype=numpy.int64)
            mids = self.fends['fends']['mid'][:]
            chr_indices = self.fends['chr_indices'][:]
            all_valid = numpy.sum(self.filter)
            min_fend = numpy.zeros((self.filter.shape[0], 2), dtype=numpy.int32)
            mids = self.fends['fends']['mid'][:]
            print >> sys.stderr, ("Done\nCopy data for fast correction..."),
            # copy needed arrays from h5dict
            if useread_int < 2:
                data = self.data['cis_data'][...]
                valid = numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]] *
                                    (mids[data[:, 1]] - mids[data[:, 0]] >= mindistance))[0]
                data = data[valid, :]
            else:
                data = None
            if useread_int > 0:
                trans_data = self.data['trans_data'][:, :]
                valid = numpy.where(self.filter[trans_data[:, 0]] * self.filter[trans_data[:, 1]])[0]
                trans_data = trans_data[valid, :]
            else:
                trans_data = None
            print >> sys.stderr, ("Done\nCount interactions for fast correction..."),
            # filter any fends with too few observed interactions
            current_fends = numpy.sum(self.filter)
            previous_fends = current_fends + 1
            while current_fends < previous_fends:
                if data is None:
                    observed_interactions = numpy.bincount(numpy.r_[trans_data[:, 0], trans_data[:, 1]],
                                                           minlength=self.filter.shape[0])
                elif trans_data is None:
                    observed_interactions = numpy.bincount(numpy.r_[data[:, 0], data[:, 1]],
                                                           minlength=self.filter.shape[0])
                else:
                    observed_interactions = numpy.bincount(numpy.r_[data[:, 0], data[:, 1], trans_data[:, 0],
                                                           trans_data[:, 1]], minlength=self.filter.shape[0])
                self.filter[numpy.where(observed_interactions < mininteractions)] = 0
                previous_fends = current_fends
                current_fends = numpy.sum(self.filter)
                if not data is None:
                    data = data[numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]])[0], :]
                if not trans_data is None:
                    trans_data = trans_data[numpy.where(self.filter[trans_data[:, 0]] *
                                                        self.filter[trans_data[:, 1]])[0], :]
            if not data is None:
                data_indices = numpy.r_[0, numpy.bincount(data[:, 0], minlength=self.filter.shape[0])].astype(numpy.int32)
                for i in range(1, data_indices.shape[0]):
                    data_indices[i] += data_indices[i - 1]
            print >> sys.stderr, ("Done\nFind MinDistance for fast correction..."),
            for i in range(chr_indices.shape[0]-1):
                # first fend outside of mindistance range
                start = chr_indices[i]
                stop = chr_indices[i]
                for j in range(chr_indices[i], chr_indices[i + 1]):
                    while start < j and mids[j] - mids[start] >= mindistance:
                        start += 1
                    min_fend[j, 0] = start
                    stop = max(stop, j + 1)
                    while stop < chr_indices[i + 1] and mids[stop] - mids[j] < mindistance:
                        stop += 1
                    min_fend[j, 1] = stop
            _distance.find_mindistance_interactions(interactions,
                                                    chr_indices,
                                                    min_fend,
                                                    self.filter,
                                                    useread_int)
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
                distance_mid_logs = numpy.log(self.distance_mids).astype(numpy.float32)
                distance_mean_logs = numpy.log(self.distance_means).astype(numpy.float32)
                _distance.find_data_distance_means(distance_means,
                                                   self.filter,
                                                   data,
                                                   data_indices,
                                                   mids,
                                                   self.distance_means,
                                                   self.distance_mean_logs,
                                                   self.distance_mids,
                                                   distance_mid_logs,
                                                   distance_mid_logs[1:] - distance_mid_logs[:-1],
                                                   0)
                if not trans_data is None:
                    total_possible = numpy.sum(self.filter) ** 2
                    for i in range(self.fends['chr_indices'].shape[0] - 1):
                        start = self.fends['chr_indices'][i]
                        stop = self.fends['chr_indices'][i + 1]
                        total_possible -= numpy.sum(self.filter[start:stop]) ** 2
                    trans_mu = (2.0 * numpy.sum(trans_data[:, 2])) / total_possible
                else:
                    trans_mu = 1.0
            print >> sys.stderr, ("Done\nFinding fend corrections...\n"),
            # calculate corrections
            for iteration in range(iterations):
                cost = _distance.find_fend_means(distance_means,
                                                 interactions,
                                                 fend_means,
                                                 data,
                                                 trans_data,
                                                 self.filter,
                                                 self.corrections,
                                                 mu,
                                                 trans_mu)
                print >> sys.stderr, ("\rIteration: %i  Cost: %f    ") % (iteration, cost),
                if (recalculate_distance > 0 and (iteration + 1) % recalculate_distance == 0 and
                        not data is None and remove_distance):
                    for i in range(1, self.num_procs):
                        self.comm.send(1, dest=i, tag=11)
                    self.find_distance_means(numbins=num_bins, minsize=minsize, smoothed=smoothed,
                                             corrected=True)
                    distance_mid_logs = numpy.log(self.distance_mids).astype(numpy.float32)
                    distance_mean_logs = numpy.log(self.distance_means).astype(numpy.float32)
                    _distance.find_data_distance_means(distance_means,
                                                       self.filter,
                                                       data,
                                                       data_indices,
                                                       mids,
                                                       self.distance_means,
                                                       self.distance_mean_logs,
                                                       self.distance_mids,
                                                       distance_mid_logs,
                                                       distance_mid_logs[1:] - distance_mid_logs[:-1],
                                                       0)
            print >> sys.stderr, ("\r%s\rFinal cost: %f\n") % (' ' * 80, cost),
            for i in range(1, self.num_procs):
                self.comm.send(0, dest=i, tag=11)      
        else:
            task = self.comm.recv(source=0, tag=11)
            while task == 1:
                self.find_distance_means(numbins=num_bins, minsize=minsize, smoothed=smoothed,
                                         corrected=True)
                task = self.comm.recv(source=0, tag=11)
        return None

    def find_trans_mean(self):
        """
        Calculate the mean signal across all valid fend-pair trans interactions.

        :returns: None
        """
        print >> sys.stderr, ("Finding mean signal across trans interactions..."),
        possible = 0
        chr_indices = self.fends['chr_indices'][...]
        for i in range(chr_indices.shape[0] - 2):
            valid1 = numpy.sum(self.filter[chr_indices[i]:chr_indices[i + 1]])
            for j in range(i + 1, chr_indices.shape[0] - 1):
                valid2 = numpy.sum(self.filter[chr_indices[j]:chr_indices[j + 1]])
                possible += valid1 * valid2
        trans_data = self.data['trans_data'][...]
        actual = numpy.sum(self.filter[trans_data[:, 0]] * self.filter[trans_data[:, 1]] * trans_data[:, 2])
        self.trans_mean = actual / float(possible)
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
            print >> sys.stderr, ("The mlpy module must be installed to use this function.")
            return None
        print >> sys.stderr, ("Learning fend-resolution 3D model..."),
        unbinned, mapping = hic_binning.unbinned_cis_signal(self, chrom, datatype='fend', arraytype='upper',
                                                            skipfiltered=True, returnmapping=True)
        mids = self.fends['fends']['mid'][mapping]
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
        print >> sys.stderr, ("Done\nModeling chromosome..."),
        pca_fast = mlpy.PCAFast(k=3)
        pca_fast.learn(data_mat)
        coordinates = pca_fast.transform(data_mat)
        print >> sys.stderr, ("Done\n"),
        return numpy.hstack((coordinates, mids.reshape(-1, 1)))

    def cis_heatmap(self, chrom, start=0, stop=None, startfend=None, stopfend=None, binsize=0, binbounds=None,
                    datatype='enrichment', arraytype='compact', maxdistance=0, skipfiltered=False, returnmapping=False,
                    dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                    removefailed=False):
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
        :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'. If returnmapping is True, a list is returned containined the requested data array and an array of associated positions (dependent on the binning options selected).
        """
        # check that all values are acceptable
        datatypes = {'raw': 0, 'fend': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
        if datatype not in datatypes:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        else:
            datatype_int = datatypes[datatype]
        if ((dynamically_binned == True and arraytype not in ['compact', 'upper']) or
            (arraytype not in ['full', 'compact', 'upper'])):
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
                                                       returnmapping=returnmapping)
            else:
                # data should be binned
                data = hic_binning.bin_cis_signal(self, chrom, start=start, stop=stop, startfend=startfend,
                                                  stopfend=stopfend, binsize=binsize, binbounds=binbounds,
                                                  datatype=datatype, arraytype=arraytype, maxdistance=maxdistance,
                                                  returnmapping=returnmapping)
        else:
            if expansion_binsize == 0:
                # data should be dynamically binned with unbinned expansion data
                expansion, fends = hic_binning.unbinned_cis_signal(self, chrom, start=start, stop=stop,
                                                                   startfend=startfend, stopfend=stopfend,
                                                                   datatype=datatype, arraytype=arraytype,
                                                                   maxdistance=maxdistance, skipfiltered=True,
                                                                   returnmapping=True)
                mids = self.fends['fends']['mid'][fends]
                binned, mapping = hic_binning.bin_cis_array(self, expansion, fends, start=start, stop=stop,
                                                            binsize=binsize, binbounds=binbounds, arraytype=arraytype,
                                                            returnmapping=True)
            else:
                # data should be dynamically binned with binned expansion data
                expansion, mapping = hic_binning.bin_cis_signal(self, chrom, start=start, stop=stop,
                                                                startfend=startfend, stopfend=stopfend,
                                                                binsize=expansion_binsize, binbounds=None,
                                                                datatype=datatype, arraytype=arraytype,
                                                                maxdistance=maxdistance, returnmapping=True)
                mids = (mapping[:, 2] + mapping[:, 3]) / 2
                binned, mapping = hic_binning.bin_cis_signal(self, chrom, start=start, stop=stop, startfend=startfend,
                                                             stopfend=stopfend, binsize=binsize, binbounds=binbounds,
                                                             datatype=datatype, arraytype=arraytype,
                                                             maxdistance=maxdistance, returnmapping=True)
            hic_binning.dynamically_bin_cis_array(expansion, mids, binned, mapping[:, 2:],
                                                  minobservations=minobservations, searchdistance=searchdistance,
                                                  removefailed=removefailed)
            if returnmapping:
                data = [binned, mapping]
            else:
                data = binned
        return data

    def trans_heatmap(self, chrom1, chrom2, start1=0, stop1=None, startfend1=None, stopfend1=None, binbounds1=None,
                      start2=0, stop2=None, startfend2=None, stopfend2=None, binbounds2=None, binsize=1000000,
                      datatype='enrichment', returnmapping=False, dynamically_binned=False, minobservations=0,
                      searchdistance=0, expansion_binsize=0, removefailed=False):
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
        :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'. If returnmapping is True, a list is returned containined the requested data array and an array of associated positions (dependent on the binning options selected).
        """
        # check that all values are acceptable
        datatypes = {'raw': 0, 'fend': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
        if datatype not in datatypes:
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
                                                returnmapping=returnmapping)
        else:
            expansion, mapping1, mapping2 = hic_binning.bin_trans_signal(self, chrom1, chrom2, start1=start1,
                                                                         stop1=stop1, startfend1=startfend1,
                                                                         stopfend1=stopfend1, binbounds1=binbounds1,
                                                                         start2=start2, stop2=stop2,
                                                                         startfend2=startfend2, stopfend2=stopfend2,
                                                                         binbounds2=binbounds2,
                                                                         binsize=expansion_binsize, datatype=datatype,
                                                                         returnmapping=True)
            mids1 = (mapping1[:, 2] + mapping1[:, 3]) / 2
            mids2 = (mapping2[:, 2] + mapping2[:, 3]) / 2
            binned, mapping1, mapping2 = hic_binning.bin_trans_signal(self, chrom1, chrom2, start1=start1, stop1=stop1,
                                                                      startfend1=startfend1, stopfend1=stopfend1,
                                                                      binbounds1=binbounds1, start2=start2,
                                                                      stop2=stop2, startfend2=startfend2,
                                                                      stopfend2=stopfend2, binbounds2=binbounds2,
                                                                      binsize=binsize, datatype=datatype,
                                                                      returnmapping=True)
            hic_binning.dynamically_bin_trans_array(expansion, mids1, mids2, binned, mapping1[:, 2:], mapping2[:, 2:],
                                                    minobservations=minobservations, searchdistance=searchdistance,
                                                    removefailed=removefailed)
            if returnmapping:
                data = [binned, mapping1, mapping2]
            else:
                data = binned
        return data

    def write_heatmap_dict(self, filename, binsize, includetrans=True, remove_distance=False, chroms=[]):
        """
        Create an h5dict file containing binned interaction arrays, bin positions, and an index of included chromosomes. This function is MPI compatible.

        :param filename: Location to write h5dict object to.
        :type filename: str.
        :param binsize: Size of bins for interaction arrays.
        :type binsize: int.
        :param includetrans: Indicates whether trans interaction arrays should be calculated and saved.
        :type includetrans: bool.
        :param remove_distance: If 'True', the expected value is calculated including the expected distance mean. Otherwise, only fend corrections are used.
        :type remove_distance: bool.
        :param chroms: A list of chromosome names indicating which chromosomes should be included. If left empty, all chromosomes are included. Optional.
        :type chroms: list
        :returns: None
        """
        hic_binning.write_heatmap_dict(self, filename, binsize, includetrans=includetrans,
                                       remove_distance=remove_distance, chroms=chroms)
        return None
