#!/usr/bin/env python

import os
import sys
from math import ceil, floor, log, exp

import numpy
from scipy.stats import linregress
import h5py

from fragment import Fragment
from fivec_data import FiveCData
from libraries._fivec_binning import find_fragment_coverage
import libraries._fivec_distance as _distance


class FiveC(object):
    """
    This is the class for handling 5C analysis.

    This class relies on :class:`Fragment <hifive.fragment.Fragment>` and :class:`FiveCData <hifive.fivec_data.FiveCData>` for genomic position and interaction count data. Use this class to perform filtering of fragments based on coverage, model fragment bias and distance dependence, and downstream analysis and manipulation. This includes binning of data, plotting of data, and statistical analysis.

    .. note::
      This class is also available as hifive.FiveC

    When initialized, this class creates an h5dict in which to store all data associated with this object.
    
    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :returns: :class:`FiveC <hifive.fivec.FiveC>` class object.
    """

    def __init__(self, filename, mode='r'):
        """
        Create a FiveC object.
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
        Load fragment-pair counts and fragment object from :class:`FiveCData <hifive.fivec_data.FiveCData>` object.

        :param filename: Specifies the file name of the :class:`FiveCData <hifive.fivec_data.FiveCData>` object to associate with this analysis.
        :type filename: str.
        :returns: None
        """
        # ensure data h5dict exists
        if not os.path.exists(filename):
            print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename.split('/')[-1]),
            return None
        self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                       os.path.dirname(self.file)), os.path.basename(filename))
        self.data = HiCData(filename).data
        fragfilename = self.data['/'].attrs['fragfilename']
        if fragfilename[:2] == './':
            fragfilename = fragfilename[2:]
        parent_count = fragfilename.count('../')
        fragfilename = '/'.join(os.path.abspath(filename).split('/')[:-(1 + parent_count)] +
                                fragfilename.lstrip('/').split('/')[parent_count:])
        self.fragfilename = "%s/%s" % (os.path.relpath(os.path.dirname(fragfilename),
                                       os.path.dirname(self.file)), os.path.basename(fragfilename))
        # ensure fend h5dict exists
        if not os.path.exists(fragfilename):
            print >> sys.stderr, ("Could not find %s.\n") % (fragfilename),
            return None
        self.frags = Fragment(fragfilename).fragments
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.frags['chromosomes']):
            self.chr2int[chrom] = i
        # create arrays
        self.filter = numpy.ones(self.frags['fragments'].shape[0], dtype=numpy.int32)
        self.corrections = numpy.zeros(self.frags['fragments'].shape[0], dtype=numpy.float32)
        return None

    def save(self):
        """
        Save analysis parameters to h5dict.
        
        :returns: None
        """
        datafile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['data', 'frags', 'file', 'chr2int']:
                continue
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                datafile.attrs[key] = self[key]
        datafile.close()
        return None

    def load(self):
        """
        Load analysis parameters from h5dict specified at object creation and open h5dicts for associated :class:`FiveCData <hifive.fivec_data.FiveCData>` and :class:`Fragment <hifive.fragment.Fragment>` objects.

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
                self.data = FiveCData(datafilename).data
        # ensure fragment h5dict exists
        if 'fragfilename' in self.__dict__:
            fragfilename = self.fragfilename
            if fragfilename[:2] == './':
                fragfilename = fragfilename[2:]
            parent_count = fragfilename.count('../')
            fragfilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                fragfilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(fragfilename):
                print >> sys.stderr, ("Could not find %s. No fragments loaded.\n") % (fragfilename),
            else:
                self.frags = Fragment(fragfilename).fragments
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.frags['chromosomes']):
            self.chr2int[chrom] = i
        datafile.close()
        return None

    def filter_fragments(self, mininteractions=20):
        """
        Iterate over the dataset and remove fragments that do not have 'minobservations' using only unfiltered fragments.

        In order to create a set of fragments that all have the necessary number of interactions, after each round of filtering, fragment interactions are retallied using only interactions that have unfiltered fragments at both ends.

        :param mininteractions: The required number of interactions for keeping a fragment in analysis.
        :type mininteractions: int.
        :returns: None
        """
        print >> sys.stderr, ("Filtering fragments..."),
        original_count = numpy.sum(self.filter)
        previous_valid = original_count + 1
        current_valid = original_count
        coverage = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
         # copy needed arrays
        data = self.data['cis_data'][...]
        indices = self.data['cis_indices'][...]
        regions = self.frags['regions'][...]
        # repeat until all remaining fragments have minobservation valid observations
        while current_valid < previous_valid:
            previous_valid = current_valid
            coverage.fill(0)
            find_fragment_coverage(data,
                                   indices,
                                   self.filter,
                                   regions['start_frag'],
                                   regions['stop_frag'],
                                   coverage,
                                   mininteractions)
            current_valid = numpy.sum(self.filter)
        print >> sys.stderr, ("Removed %i of %i fragments\n") % (original_count - current_valid, original_count),
        return None

    def find_distance_parameters(self):
        """
        Regress log counts versus inter-fragment distances to find slope and intercept values and then find the standard deviation of corrected counts.

        :returns: None
        """
        print >> sys.stderr, ("Finding distance parameters..."),
        # copy needed arrays
        data = self.data['cis_data'][...]
        mids = self.frags['fragments']['mid'][...]
        # find which pairs are both unfiltered
        valid = numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]])[0]
        # find distances between fragment pairs
        log_distances = numpy.log(mids[data[valid, 1]] - mids[data[valid, 0]])
        # find regression line
        log_counts = numpy.log(data[valid, 2]) - self.corrections[data[valid, 0]] - self.corrections[data[valid, 1]]
        temp = linregress(log_distances, log_counts)[:2]
        self.gamma = -float(temp[0])
        self.mu = float(temp[1])
        self.sigma = float(numpy.std(log_counts - self.mu + self.gamma * log_distances))
        print >> sys.stderr, ("Done\n"),
        return None

    def find_fragment_corrections(self, maxdistance=0, burnin_iterations=5000, annealing_iterations=10000,
                                  learningrate=0.01, recalculate_distance=100, display=10):
        """
         Using gradient descent, learn correction values for each valid fragment based on a Log-Normal distribution of observations.

        :param maxdistance: The maximum inter-fragment distance to be included in modeling.
        :type maxdistance: int.
        :param burnin_iterations: The number of iterations to use with constant learning rate in gradient descent for learning fragment corrections.
        :type burnin_iterations: int.
        :param annealing_iterations: The number of iterations to use with a linearly-decreasing learning rate in gradient descent for learning fragment corrections.
        :type annealing_iterations: int.
        :param learningrate: The gradient scaling factor for parameter updates.
        :type learningrate: float
        :param recalculate_distance: Number of iterations that should pass before recalculating the distance function parameters to account for the updated fragment corrections. If set to zero, no recalculation is performed.
        :type recalculate_distance: int.
        :param display: Specifies how many iterations between when cost is calculated and displayed as model is learned. If 'display' is zero, the cost is not calculated of displayed.
        :type display: int.
        :returns: None
        """
        # determine if distance parameters have been calculated
        if 'gamma' not in self.__dict__.keys():
            self.find_distance_parameters()
        print >> sys.stderr, ("Learning corrections..."),
        # copy and calculate needed arrays
        data = self.data['cis_data'][...]
        valid = numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]])[0]
        log_counts_n = numpy.log(data[:, 2] - 0.5).astype(numpy.float32)
        log_counts = numpy.log(data[:, 2]).astype(numpy.float32)
        log_counts_p = numpy.log(data[:, 2] + 0.5).astype(numpy.float32)
        indices = self.data['cis_indices'][...]
        distances = numpy.log(self.frags['fragments']['mid'][data[:, 1]] - self.frags['fragments']['mid'][data[:, 0]])
        distance_signal = (self.mu - self.gamma * distances).astype(numpy.float32)
        valid_distances = distances[valid]
        # create empty arrays
        gradients = numpy.zeros(self.filter.shape[0], dtype=numpy.float32)
        sigma_gradient = numpy.zeros(1, dtype=numpy.float32)
        # find maximum valid fragment for specified maxdistance
        max_fragment = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        _distance.find_max_frag(max_fragment,
                                      self.frags['fragments']['mid'][...],
                                      self.frags['regions']['start_frag'][...],
                                      self.frags['regions']['stop_frag'][...],
                                      maxdistance)
        # find number of interactions for each fragment
        interactions = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        _distance.find_fragment_interactions(self.filter,
                                                   data,
                                                   indices,
                                                   interactions,
                                                   max_fragment)
        all_interactions = numpy.sum(interactions) / 2
        # cycle through learning phases
        find_variance = 0
        for phase in ['burnin', 'annealing']:
            if phase == 'burnin':
                iterations = burnin_iterations
            else:
                iterations = annealing_iterations
                learningstep = learningrate / max(1, annealing_iterations)
            for iteration in range(iterations):
                # if necessary, recalculate distance parameters
                if recalculate_distance > 0 and (iteration + 1) % recalculate_distance == 0:
                    valid_counts = (log_counts[valid] - self.corrections[data[valid, 0]] -
                                    self.corrections[data[valid, 1]])
                    distance_signal = (self.mu - self.gamma * distances).astype(numpy.float32)
                    self.sigma = float(numpy.std(valid_counts - distance_signal[valid]))
                # if requested and correct iteration, indicate cost is to be calculated
                if display > 0 and iteration%display == 0:
                    find_cost = 1
                else:
                    find_cost = 0
                # find gradients
                gradients.fill(0.0)
                cost = _distance.calculate_gradients(data,
                                                     log_counts_n,
                                                     log_counts,
                                                     log_counts_p,
                                                     indices,
                                                     self.filter,
                                                     distance_signal,
                                                     self.corrections,
                                                     gradients,
                                                     max_fragment,
                                                     self.sigma,
                                                     find_cost)
                # update gradients
                _distance.update_corrections(self.filter,
                                             self.corrections,
                                             gradients,
                                             interactions,
                                             learningrate)
                # if not burnin phase, update sigma
                if phase != 'burnin':
                    learningrate -= learningstep
                # if appropriate iteration and requested, update display
                if find_cost:
                    print >> sys.stderr, ("\rLearning corrections... phase:%s iteration:%i  cost:%f ") %\
                                         (phase, iteration, cost),
        print >> sys.stderr, ("\rLearning corrections... Done %f %s\n") % (cost, ' ' * 60),
        return None

    def find_express_fragment_corrections(self, iterations=1000, remove_distance=False, recalculate_distance=100):
        """
        Using iterative approximation, learn correction values for each valid fragment.

        :param iterations: The number of iterations to use for learning fragment corrections.
        :type iterations: int.
        :param remove_distance: Specifies whether the estimated distance-dependent portion of the signal is removed prior to learning fragment corrections.
        :type remove_distance: bool.
        :param recalculate_distance: Number of iterations that should pass before recalculating the distance bin means to account for the updated fragment corrections. If set to zero, no recalculation is performed.
        :type recalculate_distance: int.
        :returns: None
        """
        print >> sys.stderr, ("Learning corrections..."),
        # copy and calculate needed arrays
        data = self.data['cis_data'][...]
        valid = numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]])[0]
        data = data[valid, :]
        log_counts = numpy.log(data[:, 2]).astype(numpy.float32)
        corrections = numpy.copy(self.corrections)
        if remove_distance:

            distances = numpy.log(self.frags['fragments']['mid'][data[:, 1]] -
                                  self.frags['fragments']['mid'][data[:, 0]]).astype(numpy.float32)
            corrected_counts = log_counts - corrections[data[:, 0]] - corrections[data[:, 1]]
            self.gamma = -linregress(distances, corrected_counts)[0]
            distance_signal = -self.gamma * distances
        else:
            distances = None
            distance_signal = None
        # create empty arrays
        fragment_means = numpy.zeros(self.filter.shape[0], dtype=numpy.float64)
        interactions = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        # find number of interactions for each fragment
        for i in range(self.frags['regions'].shape[0]):
            interactions = (numpy.bincount(data[:, 0], minlength=interactions.shape[0]) +
                            numpy.bincount(data[:, 1], minlength=interactions.shape[0])).astype(numpy.int32)
        # learn corrections
        for iteration in range(iterations):
            # if necessary, recalculate distance parameters
            if not remove_distance and recalculate_distance > 0 and (iteration + 1) % recalculate_distance == 0:
                corrected_counts = log_counts - corrections[data[:, 0]] - corrections[data[:, 1]]
                self.gamma = -linregress(distances, corrected_counts)[0]
                distance_signal = -self.gamma * distances
            # update corrections
            cost = _distance.find_fragment_means(distance_signal,
                                                 interactions,
                                                 fragment_means,
                                                 data,
                                                 log_counts,
                                                 corrections)
            print >> sys.stderr, ("\rLearning corrections... iteration:%i  cost:%f ") % (iteration, cost),
        self.corrections[:] = corrections
        print >> sys.stderr, ("\rLearning corrections... Done %f %s\n") % (cost, ' ' * 60),
        return None

    def find_trans_mean(self):
        """
        Calculate the mean signal across all valid fragment-pair trans (inter-region) interactions.
        
        :returns: None
        """
        print >> sys.stderr, ("Finding mean signal across trans interactions..."),
        possible = 0
        for i in range(self.frags['regions'].shape[0] - 1):
            valid1 = numpy.sum(self.filter[self.frags['regions']['start_frag'][i]:
                                           self.frags['regions']['stop_frag'][i]])
            for j in range(i + 1, self.frags['regions'].shape[0]):
                valid2 = numpy.sum(self.filter[self.frags['regions']['start_frag'][j]:
                                               self.frags['regions']['stop_frag'][j]])
                possible += valid1 * valid2
        trans_data = self.data['trans_data'][...]
        actual = numpy.sum(self.filter[trans_data[:, 0]] * self.filter[trans_data[:, 1]] * trans_data[:, 2])
        self.trans_mean = actual / float(possible)
        print >> sys.stderr, ('Done\n'),
        return None
