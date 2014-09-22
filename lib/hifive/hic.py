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
from libraries._hic_binning import find_fend_coverage, dynamically_bin_unbinned_upper, dynamically_bin_upper_from_upper
from hic_binning import unbinned_cis_signal
import libraries._hic_distance as _distance


class HiC(object):
    """
    This is the class for handling HiC analysis.

    This class relies on :class:`Fend` and :class:`HiCData` for genomic position and interaction count data. Use this class to perform filtering of fends based on coverage, model fend bias and distance dependence, and downstream analysis and manipulation. This includes binning of data, plotting of data, modeling of data, and statistical analysis.

        .. note::
          This class is also available as hifive.HiC

    When initialized, this class creates an h5dict in which to store all data associated with this object.

    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    """

    def __init__(self, filename, mode='r'):
        """
        Create a HiC object.
        """

        self.file = os.path.abspath(filename)
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
        Load fend-pair counts and fend object from :class:`HiCData` object.

        :param filename: Specifies the file name of the HiCData object to associate with this analysis.
        :type filename: str.
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
        """
        if 'mpi4py' in sys.modules and MPI.COMM_WORLD.Get_rank() > 0:
            return None
        datafile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['data', 'fends', 'file', 'chr2int']:
                continue
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                datafile.attrs[key] = self[key]
        datafile.close()
        return None

    def load(self):
        """
        Load analysis parameters from h5dict specified at object creation and open h5dicts for associated :class:`HiCData` and :class:`Fend` objects.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.
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
        """
        # check if MPI is available
        if 'mpi4py' in sys.modules.keys():
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()
        else:
            comm = None
            rank = 0
            num_procs = 1
        if stopfend is None:
            stopfend = self.filter.shape[0]
        self.distance_smoothing = smoothed
        if rank == 0:
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
            for i in range(1, num_procs):
                comm.send([max_distance, self.distance_bins, self.distance_mids], dest=i, tag=11)
            # divide fends amongst workers
            needed = (numpy.where(self.filter[startfend:stopfend] > 0)[0] + startfend).astype(numpy.int32)
            numpy.random.shuffle(needed)
            worker_size = int(ceil(needed.shape[0] / float(num_procs)))
            for i in range(1, num_procs):
                comm.send(needed[((i - 1) * worker_size):(i * worker_size)], dest=i, tag=11)
            node_needed = needed[((num_procs - 1) * worker_size):]
        else:
            max_distance, self.distance_bins, self.distance_mids = comm.recv(source=0, tag=11)
            node_needed = comm.recv(source=0, tag=11)
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
        if rank == 0:
            for i in range(1, num_procs):
                temp = comm.recv(source=i, tag=11)
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
            for i in range(1, num_procs):
                comm.send(self.distance_means, dest=i, tag=11)
        else:
            comm.send([bin_sums, bin_counts], dest=0, tag=11)
            self.distance_means = comm.recv(source=0, tag=11)
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
        """
        # check if MPI is available
        if 'mpi4py' in sys.modules.keys():
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()
        else:
            comm = None
            rank = 0
            num_procs = 1
        num_bins = self.distance_bins.shape[0]
        minsize = self.distance_bins[0]
        smoothed = int(self.distance_smoothing)
        # find range of fends for worker. starts and stops must always include whole fragments
        valid = numpy.where(self.filter)[0]
        worker_size = int(ceil(valid.shape[0] / float(num_procs)))
        if rank == 0:
            print >> sys.stderr, ("Find fend correction arrays..."),
            if num_procs > 1:
                start_fend = valid[worker_size * (num_procs - 1)]
                start_fend -= start_fend % 2
                stop_fend = valid[-1] + 1
                stop_fend += stop_fend % 2
            else:
                start_fend = 0
                stop_fend = self.filter.shape[0]
        else:
            start_fend = valid[worker_size * (rank - 1)]
            start_fend -= start_fend % 2
            stop_fend = valid[worker_size * rank]
            stop_fend -= stop_fend % 2
        # create arrays needed for root and find possible interactions
        if rank == 0:
            all_gradients = numpy.zeros(self.filter.shape[0], dtype=numpy.float32)
            all_corrections = numpy.copy(self.corrections).astype(numpy.float32)
            interactions = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
            _distance.find_possible_interactions(interactions,
                                                 self.filter,
                                                 self.fends['fends']['mid'][:],
                                                 self.fends['chr_indices'][...],
                                                 maxdistance)
            fend_ranges = numpy.zeros((num_procs, 3), dtype=numpy.int32)
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
        if num_procs > 1:
            if rank != 0:
                comm.send([start_fend, stop_fend, stop_fend2], dest=0, tag=11)
                gradient_arrays = None
            else:
                gradient_arrays = {}
                fend_ranges[0, 0] = start_fend
                fend_ranges[0, 1] = stop_fend
                fend_ranges[0, 2] = stop_fend2
                for i in range(1, num_procs):
                    temp = comm.recv(source=i, tag=11)
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
        _distance.find_interaction_distance_means(interaction_means,
                                                  filter,
                                                  self.fends['fends']['mid'][start_fend:stop_fend2],
                                                  self.distance_means,
                                                  self.distance_mids,
                                                  max_fend)
      # calculate correction gradients
        if rank == 0:
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
                if num_procs > 1:
                    cost = self._exchange_gradients(rank, comm, gradients, all_gradients,
                                                    fend_ranges, gradient_arrays, cost)
                else:
                    all_gradients = gradients
                if rank == 0:
                    if findcost > 0:
                        print >> sys.stderr, ("\r%s phase:%s iteration:%i  cost:%f ") %\
                                             ('Learning corrections...', phase, iteration, cost),
                    # update corrections on root only
                    _distance.update_corrections(all_corrections,
                                                 all_gradients,
                                                 interactions,
                                                 learningrate)
                # if using multiple cores, distribute needed correction values to workers
                if num_procs > 1:
                    self._exchange_corrections(rank, comm, corrections, all_corrections, fend_ranges)
                else:
                    corrections = all_corrections[start_fend:stop_fend2]
                if phase == 'annealing':
                    learningrate -= learningstep
                if recalculate_distance > 0 and (iteration + 1) % recalculate_distance == 0:
                    if rank == 0:
                        self.corrections = all_corrections.astype(numpy.float32)
                        if num_procs > 1:
                            for i in range(1, num_procs):
                                comm.Send(self.corrections, dest=i, tag=13)
                    else:
                        comm.Recv(self.corrections, source=0, tag=13)
                        self.corrections = self.corrections.astype(numpy.float32)
                    self.find_distance_means(numbins=num_bins, minsize=minsize, smoothed=smoothed,
                                             corrected=True)
                    interaction_means[:, :] = 0.0
                    _distance.find_interaction_distance_means(interaction_means,
                                                  filter,
                                                  self.fends['fends']['mid'][start_fend:stop_fend2],
                                                  self.distance_means,
                                                  self.distance_mids,
                                                  max_fend)
        if rank == 0:
            self.corrections = all_corrections.astype(numpy.float32)
            if num_procs > 1:
                for i in range(1, num_procs):
                    comm.Send(self.corrections, dest=i, tag=13)
            print >> sys.stderr, ("\rLearning corrections... Done%s\n") % (' ' * 60),
        elif num_procs > 1:
            comm.Recv(self.corrections, source=0, tag=13)
            self.corrections = self.corrections.astype(numpy.float32)
        return None

    def _exchange_gradients(self, rank, comm, gradients, all_gradients, fend_ranges, gradient_arrays, cost):
        """Compile all gradients on root"""
        if rank == 0:
            all_gradients.fill(0.0)
            all_gradients[fend_ranges[0, 0]:fend_ranges[0, 2]] += gradients
            for i in range(1, fend_ranges.shape[0]):
                cost += comm.recv(source=i, tag=11)
                comm.Recv(gradient_arrays[i], source=i, tag=13)
                all_gradients[fend_ranges[i, 0]:fend_ranges[i, 2]] += gradient_arrays[i]
        else:
            comm.send(cost, dest=0, tag=11)
            comm.Send(gradients, dest=0, tag=13)
        return cost

    def _exchange_corrections(self, rank, comm, corrections, all_corrections, fend_ranges):
        """Distribute correction values needed to each worker"""
        if rank == 0:
            for i in range(1, fend_ranges.shape[0]):
                comm.Send(all_corrections[fend_ranges[i, 0]:fend_ranges[i, 2]], dest=i, tag=13)
            corrections[:] = all_corrections[fend_ranges[0, 0]:fend_ranges[0, 2]]
        else:
            comm.Recv(corrections, source=0, tag=13)
        return None

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
        """
        num_bins = self.distance_bins.shape[0]
        minsize = self.distance_bins[0]
        smoothed = int(self.distance_smoothing)
        if mininteractions is None:
            if 'mininteractions' in self.__dict__.keys():
                mininteractions = self.mininteractions
            else:
                mininteractions = 1
        # check if MPI is available
        if 'mpi4py' in sys.modules.keys():
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()
        else:
            comm = None
            rank = 0
            num_procs = 1
        if rank == 0:
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
                _distance.find_data_distance_means(distance_means,
                                                   self.filter,
                                                   data,
                                                   data_indices,
                                                   mids,
                                                   self.distance_means,
                                                   self.distance_mids,
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
                    for i in range(1, num_procs):
                        comm.send(1, dest=i, tag=11)
                    self.find_distance_means(numbins=num_bins, minsize=minsize, smoothed=smoothed,
                                             corrected=True)
                    _distance.find_data_distance_means(distance_means,
                                                       self.filter,
                                                       data,
                                                       data_indices,
                                                       mids,
                                                       self.distance_means,
                                                       self.distance_mids,
                                                       0)
            print >> sys.stderr, ("\r%s\rFinal cost: %f\n") % (' ' * 80, cost),
            for i in range(1, num_procs):
                comm.send(0, dest=i, tag=11)      
        else:
            task = comm.recv(source=0, tag=11)
            while task == 1:
                self.find_distance_means(numbins=num_bins, minsize=minsize, smoothed=smoothed,
                                         corrected=True)
                task = comm.recv(source=0, tag=11)
        return None

    def find_trans_mean(self):
        """
        Calculate the mean signal across all valid fend-pair trans interactions.
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
        """
        if 'mlpy' not in sys.modules.keys():
            print >> sys.stderr, ("The mlpy module must be installed to use this function.")
            return None
        print >> sys.stderr, ("Learning fend-resolution 3D model..."),
        unbinned, mapping = unbinned_cis_signal(self, chrom, datatype='fend', arraytype='upper',
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
        distances = numpy.sum((coordinates.reshape(1, -1, 3) - coordinates.reshape(-1, 1, 3))**2, axis=2)**0.5
        distances = -numpy.log(distances[numpy.triu_indices(distances.shape[0], 1)])
        distances -= numpy.mean(distances)
        distances /= numpy.std(distances)
        data_mat = data_mat[numpy.triu_indices(mids.shape[0], 1)]
        data_mat -= numpy.mean(data_mat)
        data_mat /= numpy.std(data_mat)
        print numpy.corrcoef(distances, data_mat)[0, 1]
        return numpy.hstack((coordinates, mids.reshape(-1, 1)))
