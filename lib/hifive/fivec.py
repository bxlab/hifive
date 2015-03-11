#!/usr/bin/env python

import os
import sys
from math import ceil, floor, log, exp

import numpy
from scipy.stats import linregress
import h5py

from fragment import Fragment
from fivec_data import FiveCData
import libraries._fivec_binning as _binning
import libraries._fivec_optimize as _optimize
import fivec_binning


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
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`FiveC <hifive.fivec.FiveC>` class object.
    """

    def __init__(self, filename, mode='r', silent=False):
        """
        Create a FiveC object.
        """
        self.file = os.path.abspath(filename)
        self.silent = silent        
        self.gc_bins = None
        self.gc_corrections = None
        self.len_bins = None
        self.len_corrections = None
        self.corrections = None
        self.region_means = None
        self.gamma = None
        self.sigma = None
        self.trans_mean = None
        self.normalization = 'none'
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
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename.split('/')[-1]),
            return None
        self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                       os.path.dirname(self.file)), os.path.basename(filename))
        self.data = h5py.File(filename, 'r')
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
            if not self.silent:
                print >> sys.stderr, ("Could not find %s.\n") % (fragfilename),
            return None
        self.frags = h5py.File(fragfilename, 'r')
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.frags['chromosomes']):
            self.chr2int[chrom] = i
        # create arrays
        self.filter = numpy.ones(self.frags['fragments'].shape[0], dtype=numpy.int32)
        self.corrections = numpy.zeros(self.frags['fragments'].shape[0], dtype=numpy.float32)
        return None

    def save(self, out_fname=None):
        """
        Save analysis parameters to h5dict.

        :param filename: Specifies the file name of the :class:`FiveC <hifive.fivec.FiveC>` object to save this analysis to.
        :type filename: str.
        :returns: None
        """ 
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
            if 'fragfilename' in self.__dict__:
                fragfilename = self.fragfilename
                if fragfilename[:2] == './':
                    fragfilename = fragfilename[2:]
                parent_count = fragfilename.count('../')
                fragfilename = '/'.join(original_file.split('/')[:-(1 + parent_count)] +
                                        fragfilename.lstrip('/').split('/')[parent_count:])
                self.fragfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fragfilename)),
                                               os.path.dirname(self.file)), os.path.basename(fragfilename))
        datafile = h5py.File(self.file, 'w')
        for key in self.__dict__.keys():
            if key in ['data', 'frags', 'file', 'chr2int', 'silent']:
                continue
            elif self[key] is None:
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
                if not self.silent:
                    print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (datafilename),
            else:
                self.data = h5py.File(datafilename, 'r')
        # ensure fragment h5dict exists
        if 'fragfilename' in self.__dict__:
            fragfilename = self.fragfilename
            if fragfilename[:2] == './':
                fragfilename = fragfilename[2:]
            parent_count = fragfilename.count('../')
            fragfilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                fragfilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(fragfilename):
                if not self.silent:
                    print >> sys.stderr, ("Could not find %s. No fragments loaded.\n") % (fragfilename),
            else:
                self.frags = Fragment(fragfilename).fragments
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.frags['chromosomes']):
            self.chr2int[chrom] = i
        datafile.close()
        return None

    def filter_fragments(self, mininteractions=20, mindistance=0, maxdistance=0):
        """
        Iterate over the dataset and remove fragments that do not have 'minobservations' using only unfiltered fragments and interactions falling with the distance limits specified.

        In order to create a set of fragments that all have the necessary number of interactions, after each round of filtering, fragment interactions are retallied using only interactions that have unfiltered fragments at both ends.

        :param mininteractions: The required number of interactions for keeping a fragment in analysis.
        :type mininteractions: int.
        :param mindistance: The minimum inter-fragment distance to be included in filtering.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fragment distance to be included in filtering. A value of zero indicates no maximum cutoff.
        :type maxdistance: int.
        :returns: None
        """
        if not self.silent:
            print >> sys.stderr, ("Filtering fragments..."),
        original_count = numpy.sum(self.filter)
        previous_valid = original_count + 1
        current_valid = original_count
        coverage = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
         # copy needed arrays
        data = self.data['cis_data'][...]
        distances = self.frags['fragments']['mid'][data[:, 1]] - self.frags['fragments']['mid'][data[:, 0]]
        if maxdistance == 0:
            maxdistance = numpy.amax(distances) + 1
        valid = numpy.where((self.filter[data[:, 0]] * self.filter[data[:, 1]]) *
                            (distances >= mindistance) * (distances < maxdistance))[0]
        data = data[valid, :]
        # repeat until all remaining fragments have minobservation valid observations
        while current_valid < previous_valid:
            previous_valid = current_valid
            coverage = numpy.bincount(data[:, 0], minlength=self.filter.shape[0])
            coverage += numpy.bincount(data[:, 1], minlength=self.filter.shape[0])
            invalid = numpy.where(coverage < mininteractions)[0]
            self.filter[invalid] = 0
            valid = numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]])[0]
            data = data[valid, :]
            current_valid = numpy.sum(self.filter)
        if not self.silent:
            print >> sys.stderr, ("Removed %i of %i fragments\n") % (original_count - current_valid, original_count),
        return None

    def find_distance_parameters(self):
        """
        Regress log counts versus inter-fragment distances to find slope and intercept values and then find the standard deviation of corrected counts.

        :returns: None
        """
        if not self.silent:
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
        self.region_means = numpy.zeros(self.frags['regions'].shape[0], dtype=numpy.float32) + temp[1]
        self.sigma = float(numpy.std(log_counts - self.mu + self.gamma * log_distances))
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        return None

    def find_probability_fragment_corrections(self, mindistance=0, maxdistance=0, burnin_iterations=5000,
                                              annealing_iterations=10000, learningrate=0.1, precalculate=True,
                                              regions=[], precorrect=False, display=10):
        """
         Using gradient descent, learn correction values for each valid fragment based on a Log-Normal distribution of observations.

        :param mindistance: The minimum inter-fragment distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fragment distance to be included in modeling.
        :type maxdistance: int.
        :param burnin_iterations: The number of iterations to use with constant learning rate in gradient descent for learning fragment corrections.
        :type burnin_iterations: int.
        :param annealing_iterations: The number of iterations to use with a linearly-decreasing learning rate in gradient descent for learning fragment corrections.
        :type annealing_iterations: int.
        :param learningrate: The gradient scaling factor for parameter updates.
        :type learningrate: float
        :param precalculate: Specifies whether the correction values should be initialized at the fragment means.
        :type precalculate: bool.
        :param regions: A list of regions to calculate corrections for. If set as None, all region corrections are found.
        :type regions: list
        :param precorrect: Use regression-based corrections in expected value calculations, resulting in a chained normalization approach.
        :type precorrect: bool.
        :param display: Specifies how many iterations between when cost is calculated and displayed as model is learned. If 'display' is zero, the cost is not calculated of displayed.
        :type display: int.
        :returns: None
        """
        if precorrect and self.gc_corrections is None and self.len_corrections is None:
            if not self.silent:
                print >> sys.stderr, ("Precorrection can only be used in project has previously run 'find_regression_fragment_corrections'.\n"),
            return None
        if self.corrections is None:
            self.corrections = numpy.zeros(self.frags['fragments'].shape[0], dtype=numpy.float32)
        # if regions not given, set to all regions
        if regions == None or len(regions) == 0:
            regions = numpy.arange(self.frags['regions'].shape[0])
        # determine if distance parameters have been calculated
        if self.gamma is None:
            self.find_distance_parameters()
        if precorrect:
            if not self.gc_corrections is None:
                gc_bins = int((0.25 + 2 * self.gc_corrections.shape[0]) ** 0.5 - 0.5)
                gc_corrections = numpy.zeros((gc_bins, gc_bins), dtype=numpy.float32)
                indices = numpy.triu_indices(gc_bins, 0)
                gc_corrections[indices] = self.gc_corrections
                gc_corrections[indices[1], indices[0]] = gc_corrections[indices]
            else:
                gc_corrections = None
                gc_indices = None
            if not self.len_corrections is None:
                len_bins = int((0.25 + 2 * self.len_corrections.shape[0]) ** 0.5 - 0.5)
                len_corrections = numpy.zeros((len_bins, len_bins), dtype=numpy.float32)
                indices = numpy.triu_indices(len_bins, 0)
                len_corrections[indices] = self.len_corrections
                len_corrections[indices[1], indices[0]] = len_corrections[indices]
            else:
                len_corrections = None
                len_indices = None
        # limit corrections to only requested regions
        filt = numpy.copy(self.filter)
        for i in range(self.frags['regions'].shape[0]):
            if i not in regions:
                filt[self.frags['regions']['start_frag'][i]:self.frags['regions']['stop_frag'][i]] = 0
        # copy and calculate needed arrays
        if not self.silent:
            print >> sys.stderr, ("\r%s\rCopying needed data...") % (' ' * 80),
        data = self.data['cis_data'][...]
        distances = self.frags['fragments']['mid'][data[:, 1]] - self.frags['fragments']['mid'][data[:, 0]]
        if maxdistance == 0:
            maxdistance = numpy.amax(distances) + 1
        valid = numpy.where((filt[data[:, 0]] * filt[data[:, 1]]) *
                            (distances >= mindistance) * (distances < maxdistance))[0]
        data = data[valid, :]
        distances = numpy.log(distances[valid])
        log_counts_n = numpy.log(data[:, 2] - 0.5).astype(numpy.float32)
        log_counts = numpy.log(data[:, 2]).astype(numpy.float32)
        log_counts_p = numpy.log(data[:, 2] + 0.5).astype(numpy.float32)
        distance_signal = (-self.gamma * distances).astype(numpy.float32)
        for i in range(self.frags['regions'].shape[0]):
            distance_signal[self.frags['regions']['start_frag'][i]:self.frags['regions']['stop_frag'][i]] += (
                self.region_means[self.frags['regions']['index'][i]])
        # create empty arrays
        gradients = numpy.zeros(self.filter.shape[0], dtype=numpy.float32)
        sigma_gradient = numpy.zeros(1, dtype=numpy.float32)
        # find number of interactions for each fragment
        interactions = numpy.bincount(data[:, 0], minlength=self.filter.shape[0]).astype(numpy.int32)
        interactions += numpy.bincount(data[:, 1], minlength=self.filter.shape[0]).astype(numpy.int32)
        all_interactions = numpy.sum(interactions) / 2
        # if precalculation requested, find fragment means
        if precalculate:
            enrichments = log_counts / distance_signal
            count_sums = numpy.bincount(data[:, 0], weights=enrichments, minlength=gradients.shape[0])
            count_sums += numpy.bincount(data[:, 1], weights=enrichments, minlength=gradients.shape[0])
            self.corrections = ((count_sums / numpy.maximum(1, interactions)) ** 0.5).astype(numpy.float32)
        if precorrect:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding regression corrections...") % (' ' * 80),
            if not gc_corrections is None:
                gc_indices = numpy.searchsorted(self.gc_bins,
                             self.frags['fragments']['gc'][...]).astype(numpy.int32)
            if not len_corrections is None:
                len_indices = numpy.searchsorted(self.len_bins,
                              self.fends['fragments']['stop'][...] -
                              self.fends['fragments']['start'][...]).astype(numpy.int32)
            _optimize.find_regression_correction_adjustment(distance_signal,
                                                            data,
                                                            gc_corrections,
                                                            len_corrections,
                                                            gc_indices,
                                                            len_indices)
        # cycle through learning phases
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning corrections...") % (' ' * 80),
        find_variance = 0
        for phase in ['burnin', 'annealing']:
            learningstep = learningrate / max(1, annealing_iterations)
            if (phase == 'burnin' and burnin_iterations == 0) or (phase == 'annealing' and annealing_iterations == 0):
                cont = False
            else:
                cont = True
            iteration = 0
            while cont:
                iteration += 1
                # if requested and correct iteration, indicate cost is to be calculated
                if display > 0 and iteration%display == 0:
                    find_cost = 1
                else:
                    find_cost = 0
                # find gradients
                gradients.fill(0.0)
                cost = _optimize.calculate_gradients(data,
                                                     log_counts_n,
                                                     log_counts,
                                                     log_counts_p,
                                                     distance_signal,
                                                     self.corrections,
                                                     gradients,
                                                     self.sigma,
                                                     find_cost)
                # update gradients
                _optimize.update_corrections(filt,
                                             self.corrections,
                                             gradients,
                                             interactions,
                                             learningrate)
                # if appropriate iteration and requested, update display
                if find_cost and not self.silent:
                    print >> sys.stderr, ("\r%s\rLearning corrections... phase:%s iteration:%i cost:%06f") %\
                                         (' ' * 80, phase, iteration, cost),
                if phase == 'annealing':
                    learningrate -= learningstep
                    if iteration >= annealing_iterations:
                        cont = False
                elif iteration >= burnin_iterations:
                    cont = False
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning corrections... Final Cost: %f  Done\n") % (' ' * 80, cost),
        # Calculate region means
        data = self.data['cis_data'][...]
        data = data[numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]])[0], :]
        for i in regions:
            where = numpy.where((data[:, 0] >= self.frags['regions']['start_frag'][i]) *
                                (data[:, 0] < self.frags['regions']['stop_frag'][i]))[0]
            corrected_sum = numpy.sum(corrections[data[where, 0]] + corrections[data[where, 1]])
            self.region_means[i] += corrected_sum
            where = numpy.where(self.filter[self.frags['regions']['start_frag']:self.frags['regions']['stop_frag']] ==
                                1)[0] + self.frags['regions']['start_frag']
            self.corrections[where] -= corrected_sum / 2.0
        if precorrect:
            self.normalization = 'regression-probability'
        else:
            self.normalization = 'probability'
        return None

    def find_express_fragment_corrections(self, mindistance=0, maxdistance=0, iterations=1000, remove_distance=False,
                                          usereads='cis', regions=[], precorrect=False):
        """
        Using iterative approximation, learn correction values for each valid fragment.

        :param mindistance: The minimum inter-fragment distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fragment distance to be included in modeling.
        :type maxdistance: int.
        :param iterations: The number of iterations to use for learning fragment corrections.
        :type iterations: int.
        :param remove_distance: Specifies whether the estimated distance-dependent portion of the signal is removed prior to learning fragment corrections.
        :type remove_distance: bool.
        :param usereads: Specifies which set of interactions to use, 'cis', 'trans', or 'all'.
        :type usereads: str.
        :param regions: A list of regions to calculate corrections for. If set as None, all region corrections are found.
        :type regions: list
        :param precorrect: Use regression-based corrections in expected value calculations, resulting in a chained normalization approach.
        :type precorrect: bool.
        :returns: None
        """
        if precorrect and self.gc_corrections is None and self.len_corrections is None:
            if not self.silent:
                print >> sys.stderr, ("Precorrection can only be used in project has previously run 'find_regression_fragment_corrections'.\n"),
            return None
        # make sure usereads has a valid value
        if usereads not in ['cis', 'trans', 'all']:
            if not self.silent:
                print >> sys.stderr, ("'usereads' does not have a valid value.\n"),
            return None
        # if regions not given, set to all regions
        if regions == None or len(regions) == 0:
            regions = numpy.arange(self.frags['regions'].shape[0])
        if self.corrections is None:
            self.corrections = numpy.zeros(self.frags['fragments'].shape[0], dtype=numpy.float32)
        if precorrect:
            if not self.gc_corrections is None:
                gc_bins = int((0.25 + 2 * self.gc_corrections.shape[0]) ** 0.5 - 0.5)
                gc_corrections = numpy.zeros((gc_bins, gc_bins), dtype=numpy.float32)
                indices = numpy.triu_indices(gc_bins, 0)
                gc_corrections[indices] = self.gc_corrections
                gc_corrections[indices[1], indices[0]] = gc_corrections[indices]
            else:
                gc_corrections = None
                gc_indices = None
            if not self.len_corrections is None:
                len_bins = int((0.25 + 2 * self.len_corrections.shape[0]) ** 0.5 - 0.5)
                len_corrections = numpy.zeros((len_bins, len_bins), dtype=numpy.float32)
                indices = numpy.triu_indices(len_bins, 0)
                len_corrections[indices] = self.len_corrections
                len_corrections[indices[1], indices[0]] = len_corrections[indices]
            else:
                len_corrections = None
                len_indices = None
        # limit corrections to only requested regions
        filt = numpy.copy(self.filter)
        for i in range(self.frags['regions'].shape[0]):
            if i not in regions:
                filt[self.frags['regions']['start_frag'][i]:self.frags['regions']['stop_frag'][i]] = 0
        if not self.silent:
            print >> sys.stderr, ("\r%s\rCopying needed data...") % (' ' * 80),
        # copy and calculate needed arrays
        data = None
        trans_data = None
        log_counts = None
        trans_log_counts = None
        distance_signal = None
        trans_signal = None
        corrections = numpy.copy(self.corrections)
        if usereads in ['cis', 'all']:
            data = self.data['cis_data'][...]
            distances = (self.frags['fragments']['mid'][data[:, 1]] -
                         self.frags['fragments']['mid'][data[:, 0]]).astype(numpy.float32)
            if maxdistance == 0:
                maxdistance = numpy.amax(distances) + 1
            valid = numpy.where((filt[data[:, 0]] * filt[data[:, 1]]) *
                                (distances >= mindistance) * (distances < maxdistance))[0]
            data = data[valid, :]
            log_counts = numpy.log(data[:, 2]).astype(numpy.float32)
            if remove_distance:
                if self.gamma is None:
                    self.find_distance_parameters()
                distance_signal = (-self.gamma * numpy.log(distances)).astype(numpy.float32)
                for i in range(self.frags['regions'].shape[0]):
                    distance_signal[self.frags['regions']['start_frag'][i]:self.frags['regions']['stop_frag'][i]] += (
                        self.region_means[self.frags['regions']['index'][i]])
                distance_signal = distance_signal[valid]
        if usereads in ['trans', 'all']:
            trans_data =  self.data['trans_data'][...]
            valid = numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0]
            trans_data = trans_data[valid, :]
            trans_log_counts = numpy.log(trans_data[:, 2]).astype(numpy.float32)
            if remove_distance:
                if self.trans_mean is None:
                    self.find_trans_mean()
                trans_signal = numpy.zeros(trans_data.shape[0], dtype=numpy.float32) + self.trans_mean
        if precorrect:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding regression corrections...") % (' ' * 80),
            if not gc_corrections is None:
                gc_indices = numpy.searchsorted(self.gc_bins,
                             self.frags['fragments']['gc'][...]).astype(numpy.int32)
            if not len_corrections is None:
                len_indices = numpy.searchsorted(self.len_bins,
                              self.fends['fragments']['stop'][...] -
                              self.fends['fragments']['start'][...]).astype(numpy.int32)
            if usereads in ['cis', 'all']:
                if distance_signal is None:
                    distance_signal = numpy.zeros(data.shape[0], dtype=numpy.float32)
                _optimize.find_regression_correction_adjustment(distance_signal,
                                                                data,
                                                                gc_corrections,
                                                                len_corrections,
                                                                gc_indices,
                                                                len_indices)
            if usereads in ['trans', 'all']:
                if trans_signal is None:
                    trans_signal = numpy.zeros(trans_data.shape[0], dtype=numpy.float32)
                _optimize.find_regression_correction_adjustment(trans_signal,
                                                                trans_data,
                                                                gc_corrections,
                                                                len_corrections,
                                                                gc_indices,
                                                                len_indices)
        # create empty arrays
        fragment_means = numpy.zeros(self.filter.shape[0], dtype=numpy.float64)
        interactions = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
        # find number of interactions for each fragment
        for i in range(self.frags['regions'].shape[0]):
            if not data is None:
                interactions += (numpy.bincount(data[:, 0], minlength=interactions.shape[0]) +
                                 numpy.bincount(data[:, 1], minlength=interactions.shape[0])).astype(numpy.int32)
            if not trans_data is None:
                interactions += (numpy.bincount(trans_data[:, 0], minlength=interactions.shape[0]) +
                                 numpy.bincount(trans_data[:, 1], minlength=interactions.shape[0])).astype(numpy.int32)
        # learn corrections
        for iteration in range(iterations):
            # update corrections
            cost = _optimize.find_fragment_means(distance_signal,
                                                 trans_signal,
                                                 interactions,
                                                 fragment_means,
                                                 data,
                                                 trans_data,
                                                 log_counts,
                                                 trans_log_counts,
                                                 corrections)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLearning corrections... iteration:%i  cost:%f ") % (' ' * 80, iteration,
                                                                                                 cost),
        where = numpy.where(filt)[0]
        self.corrections[where] = corrections[where]
        # Calculate region means
        data = self.data['cis_data'][...]
        data = data[numpy.where(self.filter[data[:, 0]] * self.filter[data[:, 1]])[0], :]
        for i in regions:
            where = numpy.where((data[:, 0] >= self.frags['regions']['start_frag'][i]) *
                                (data[:, 0] < self.frags['regions']['stop_frag'][i]))[0]
            corrected_sum = numpy.sum(corrections[data[where, 0]] + corrections[data[where, 1]])
            self.region_means[i] += corrected_sum
            where = numpy.where(self.filter[self.frags['regions']['start_frag']:self.frags['regions']['stop_frag']]
                    )[0] + self.frags['regions']['start_frag']
            self.corrections[where] -= corrected_sum / 2.0
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning corrections... Final Cost: %f  Done\n") % (' ' * 80, cost),
        if precorrect:
            self.normalization = 'regression-express'
        else:
            self.normalization = 'express'
        return None

    def find_regression_fragment_corrections(self, mindistance=0, maxdistance=0, num_bins=[10, 10, 10],
                                             model=['gc', 'len', 'distance'], learning_threshold=1.0,
                                             max_iterations=100,  usereads='cis', regions=[]):
        """
        Using multinomial regression, learn correction values for combinations of model parameter bins.

        :param mindistance: The minimum inter-fend distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fend distance to be included in modeling.
        :type maxdistance: int.
        :param num_bins: The number of approximately equal-sized bins two divide model components into.
        :type num_bins: int.
        :param remove_distance: Use distance dependence curve in prior probability calculation for each observation.
        :type remove_distance: bool.
        :param model: A list of fend features to be used in model. Valid values are 'gc', 'len', and 'distance'. The 'distance' parameter is only good with 'cis' or 'all' reads. If used with 'all', distances will be partitioned into n - 1 bins and the final distance bin will contain all trans data.
        :type model: list
        :param learning_threshold: The minimum change in log-likelihood needed to continue iterative learning process.
        :type learning_threshold: float
        :param max_iterations: The maximum number of iterations to use for learning model parameters.
        :type max_iterations: int.
        :param usereads: Specifies which set of interactions to use, 'cis', 'trans', or 'all'.
        :type usereads: str.
        :param regions: A list of regions to calculate corrections for. If set as None, all region corrections are found.
        :type regions: list
        :returns: None
        """
        present = True
        for parameter in model:
            if not parameter in ['len', 'distance'] and parameter not in self.frags['fragments'].dtype.names:
                if not self.silent:
                    print >> sys.stderr, ("Fragment feature %s not found in fragment object. Try removing it from model or creating a new fragment object with feature data.\n"),
                return None
        if 'distance' in model and usereads == 'trans':
            if not self.silent:
                print >> sys.stderr, ("The 'distance' parameter can only be used in conjuction with 'cis' or 'all' reads.\n"),
            return None
        if len(model) != len(num_bins):
            if not self.silent:
                print >> sys.stderr, ("The number of items in the 'model' parameter must be the same as the number in the 'num_bins' parameter.\n"),
            return None
        # make sure usereads has a valid value
        if usereads not in ['cis', 'trans', 'all']:
            if not self.silent:
                print >> sys.stderr, ("'usereads' does not have a valid value.\n"),
            return None
        # if regions not given, set to all regions
        if regions == None or len(regions) == 0:
            regions = numpy.arange(self.frags['regions'].shape[0])
        # limit corrections to only requested regions
        filt = numpy.copy(self.filter)
        for i in range(self.frags['regions'].shape[0]):
            if i not in regions:
                filt[self.frags['regions']['start_frag'][i]:self.frags['regions']['stop_frag'][i]] = 0
        if maxdistance == 0:
            for i in range(self.frags['regions'].shape[0]):
                maxdistance = max(maxdistance, self.frags['regions']['stop'][i] -
                                               self.frags['regions']['start'][i]) + 1
        valid = numpy.where(filt == 1)[0]
        # Create requested bin cutoffs for each feature of the regression model
        if not self.silent:
            print >> sys.stderr, ("\r%s\rPartitioning features into bins...") % (' ' * 80),
        total_bins = 1
        if 'gc' in model:
            gc_bins = num_bins[model.index('gc')]
            gc_values = self.frags['fragments']['gc'][...][valid]
            gc_values.sort()
            self.gc_bins = gc_values[numpy.round(numpy.linspace(0, gc_values.shape[0],
                                     gc_bins + 1)).astype(numpy.int32)[1:] - 1]
            self.gc_bins[-1] = numpy.inf
            self.gc_corrections = numpy.zeros(gc_bins * (gc_bins + 1) / 2, dtype=numpy.float64)
            total_bins *= self.gc_corrections.shape[0]
            gc_div = 1
            gc_indices = numpy.searchsorted(self.gc_bins, self.frags['fragments']['gc'][...]).astype(numpy.int32)
        else:
            gc_div = 0
            gc_indices = None
            gc_bins = 0
        if 'len' in model:
            len_bins = num_bins[model.index('len')]
            len_values = (self.frags['fragments']['stop'][...][valid] -
                          self.frags['fragments']['start'][...][valid]).astype(numpy.float32)
            len_values.sort()
            self.len_bins = len_values[numpy.round(numpy.linspace(0, len_values.shape[0],
                                       len_bins + 1)).astype(numpy.int32)[1:] - 1]
            self.len_bins[-1] = numpy.inf
            self.len_corrections = numpy.zeros(len_bins * (len_bins + 1) / 2, dtype=numpy.float64)
            len_div = total_bins
            total_bins *= self.len_corrections.shape[0]
            len_indices = numpy.searchsorted(self.len_bins, self.frags['fragments']['stop'][...] -
                                             self.frags['fragments']['start'][...]).astype(numpy.int32)
        else:
            len_div = 0
            len_indices = None
            len_bins = 0
        mids = self.frags['fragments']['mid'][...]
        if 'distance' in model:
            distance_bins = num_bins[model.index('distance')]
            if maxdistance == 0:
                max_dist = 0
                dists = mids[chr_indices[1:] - 1] - mids[chr_indices[:-1]]
                for chrom in chroms:
                    max_dist = max(max_dist, dists[self.chr2int[chrom]])
            else:
                max_dist = maxdistance
            if usereads == 'cis':
                distance_cutoffs = numpy.ceil(numpy.exp(numpy.linspace(numpy.log(max(1, mindistance)),
                                              numpy.log(max_dist), distance_bins + 1))[1:]).astype(numpy.int32)
            else:
                distance_cutoffs = numpy.ceil(numpy.exp(numpy.linspace(numpy.log(max(1, mindistance)),
                                              numpy.log(max_dist), distance_bins))[1:]).astype(numpy.int32)
            distance_cutoffs[-1] += 1
            distance_corrections = numpy.zeros(distance_bins, dtype=numpy.float64)
            distance_div = total_bins
            total_bins *= distance_bins
        else:
            distance_cutoffs = None
            distance_corrections = None
            distance_div = 0
            distance_bins = 0
        bin_counts = numpy.zeros((total_bins, 2), dtype=numpy.int64)
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinding bin counts...") % (' ' * 80),
        # Find number of observations in each bin
        strands = self.frags['fragments']['strand'][...]
        data = None
        trans_data = None
        if usereads in ['cis', 'all']:
            data = self.data['cis_data'][...]
            data = data[numpy.where(filt[data[:, 0]] * filt[data[:, 1]])[0], :]
        if usereads in ['trans', 'all']:
            trans_data = self.data['trans_data'][...]
            trans_data = trans_data[numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0], :]
        _binning.regression_bin_observed(data,
                                         trans_data,
                                         mids,
                                         bin_counts,
                                         gc_indices,
                                         len_indices,
                                         distance_cutoffs,
                                         gc_div,
                                         len_div,
                                         distance_div,
                                         gc_bins,
                                         len_bins,
                                         distance_bins,
                                         mindistance,
                                         maxdistance)
        if usereads in ['cis', 'all']:
            for i in regions:
                # Find number of possible interactions in each bin
                startfrag = self.frags['regions']['start_frag'][i]
                stopfrag = self.frags['regions']['stop_frag'][i]
                _binning.regression_bin_cis_expected(filt,
                                                     mids,
                                                     strands,
                                                     bin_counts,
                                                     gc_indices,
                                                     len_indices,
                                                     distance_cutoffs,
                                                     gc_div,
                                                     len_div,
                                                     distance_div,
                                                     gc_bins,
                                                     len_bins,
                                                     distance_bins,
                                                     mindistance,
                                                     maxdistance,
                                                     startfrag,
                                                     stopfrag)
        if usereads in ['trans', 'all']:
            for i in range(len(regions) - 1):
                # Find number of possible interactions in each bin
                startfrag1 = self.frags['regions']['start_frag'][regions[i]]
                stopfrag1 = self.frags['regions']['stop_frag'][regions[i]]
                for j in range(i + 1, len(regions)):
                    startfrag1 = self.frags['regions']['start_frag'][regions[j]]
                    stopfrag1 = self.frags['regions']['stop_frag'][regions[j]]
                    _binning.regression_bin_trans_expected(filt,
                                                           mids,
                                                           strands,
                                                           bin_counts,
                                                           gc_indices,
                                                           len_indices,
                                                           gc_div,
                                                           len_div,
                                                           distance_div,
                                                           gc_bins,
                                                           len_bins,
                                                           distance_bins,
                                                           mindistance,
                                                           maxdistance,
                                                           startfrag1,
                                                           stopfrag1,
                                                           startfrag2,
                                                           stopfrag2)
        # Find seed values
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinding seed values...") % (' ' * 80),
        prior = numpy.sum(bin_counts[:, 0]) / numpy.sum(bin_counts[:, 1]).astype(numpy.float64)
        log_prior = numpy.log2(prior)
        log_2 = numpy.log(2.0)
        if gc_bins > 0:
            gc_indices = numpy.arange(bin_counts.shape[0], dtype=numpy.int32) % self.gc_corrections.shape[0]
            self.gc_corrections = (
                numpy.bincount(gc_indices, weights=bin_counts[:, 0],
                               minlength=self.gc_corrections.shape[0]).astype(numpy.float64) /
                (numpy.maximum(1, numpy.bincount(gc_indices, weights=bin_counts[:, 1],
                               minlength=self.gc_corrections.shape[0])) *
                prior)).astype(numpy.float64)
        else:
            gc_indices = None
        if len_bins > 0:
            len_indices = ((numpy.arange(bin_counts.shape[0], dtype=numpy.int32) / len_div) %
                           self.len_corrections.shape[0])
            self.len_corrections = (
                numpy.bincount(len_indices, weights=bin_counts[:, 0],
                minlength=self.len_corrections.shape[0]).astype(numpy.float64) /
                (numpy.maximum(1, numpy.bincount(len_indices, weights=bin_counts[:, 1],
                minlength=self.len_corrections.shape[0])) * prior)).astype(numpy.float64)
        else:
            len_indices = None
        if distance_bins > 0:
            distance_indices = ((numpy.arange(bin_counts.shape[0], dtype=numpy.int32) / distance_div) %
                                distance_corrections.shape[0])
            distance_corrections = (
                numpy.bincount(distance_indices, weights=bin_counts[:, 0],
                minlength=distance_corrections.shape[0]).astype(numpy.float64) /
                (numpy.maximum(1, numpy.bincount(distance_indices, weights=bin_counts[:, 1],
                minlength=distance_corrections.shape[0])) * prior)).astype(numpy.float64)
        else:
            distance_indices = None

        def find_ll(gc_c, gc_i, len_c, len_i, distance_c, distance_i, counts, prior):
            sum_log = numpy.zeros(counts.shape[0], dtype=numpy.float64)
            prod = numpy.ones(counts.shape[0], dtype=numpy.float64)
            if not gc_c is None:
                temp = gc_c[gc_i]
                sum_log += numpy.log2(temp)
                prod *= temp
            if not len_c is None:
                temp = len_c[len_i]
                sum_log += numpy.log2(temp)
                prod *= temp
            if not distance_c is None:
                temp = distance_c[distance_i]
                sum_log += numpy.log2(temp)
                prod *= temp
            return (-numpy.sum(counts[:, 0] * (numpy.log2(prior) + sum_log) +
                    (counts[:, 1] - counts[:, 0]) * numpy.log2(1.0 - prior * prod)))

        ll = find_ll(self.gc_corrections, gc_indices, self.len_corrections, len_indices, distance_corrections,
                     distance_indices, bin_counts, prior)
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning regression corrections... iteration:00  ll:%f") % (' ' * 80, ll),
        iteration = 0
        delta = numpy.inf
        pgtol = 1e-8

        def temp_ll(x, *args):
            counts, sum_log, prod, prior, log_prior, log_2 = args[:6]
            return -numpy.sum(counts[:, 0] * (log_prior + sum_log + numpy.log2(x[0])) +
                             (counts[:, 1] - counts[:, 0]) * numpy.log2(1.0 - prior * prod * x[0]))

        def temp_ll_grad(x, *args):
            counts, sum_log, prod, prior, log_prior, log_2 = args[:6]
            grad = numpy.array(-numpy.sum(counts[:, 0] / (log_2 * x[0]) - (counts[:, 1] - counts[:, 0]) *
                               prior * prod / (log_2 * (1.0 - prior * prod * x[0]))), dtype=numpy.float64)
            return grad

        numpy.seterr(invalid='ignore')
        while iteration < max_iterations and delta >= learning_threshold:
            if not self.gc_corrections is None:
                new_gc_corrections = numpy.zeros(self.gc_corrections.shape[0], dtype=numpy.float64)
                for i in range(selc.gc_corrections.shape[0]):
                    where = numpy.where(gc_indices == i)[0]
                    temp_bin_counts = bin_counts[where, :]
                    temp_sum_log = numpy.zeros(where.shape[0], dtype=numpy.float64)
                    temp_prod = numpy.ones(where.shape[0], dtype=numpy.float64)
                    if not self.len_corrections is None:
                        temp = self.len_corrections[len_indices[where]]
                        temp_sum_log += numpy.log2(temp)
                        temp_prod *= temp
                    if not distance_corrections is None:
                        temp = distance_corrections[distance_indices[where]]
                        temp_sum_log += numpy.log2(temp)
                        temp_prod *= temp
                    x0 = self.gc_corrections[i:(i+1)]
                    x, f, d = bfgs(func=temp_ll, x0=x0, fprime=temp_ll_grad, pgtol=pgtol,
                                   args=(temp_bin_counts, temp_sum_log, temp_prod, prior, log_prior, log_2))
                    new_gc_corrections[i] = x[0]
            if not self.len_corrections is None:
                for i in range(self.len_corrections.shape[0]):
                    where = numpy.where(len_indices == i)[0]
                    temp_bin_counts = bin_counts[where, :]
                    temp_sum_log = numpy.zeros(where.shape[0], dtype=numpy.float64)
                    temp_prod = numpy.ones(where.shape[0], dtype=numpy.float64)
                    if not self.map_corrections is None:
                        temp = self.map_corrections[map_indices[where]]
                        temp_sum_log += numpy.log2(temp)
                        temp_prod *= temp
                    if not self.gc_corrections is None:
                        temp = self.gc_corrections[gc_indices[where]]
                        temp_sum_log += numpy.log2(temp)
                        temp_prod *= temp
                    if not distance_corrections is None:
                        temp = distance_corrections[distance_indices[where]]
                        temp_sum_log += numpy.log2(temp)
                        temp_prod *= temp
                    x0 = self.len_corrections[i:(i+1)]
                    x, f, d = bfgs(func=temp_ll, x0=x0, fprime=temp_ll_grad, pgtol=pgtol,
                                   args=(temp_bin_counts, temp_sum_log, temp_prod, prior, log_prior, log_2))
                    new_len_corrections[i] = x[0]
            if not self.gc_corrections is None:
                self.gc_corrections = new_gc_corrections
            if not self.len_corrections is None:
                self.len_corrections = new_len_corrections
            iteration += 1
            new_ll = find_ll(self.gc_corrections, gc_indices, self.len_corrections, len_indices,
                             distance_corrections, distance_indices, bin_counts, prior)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLearning regression corrections... iteration:%02i  ll:%f") % (' ' * 80,
                                      iteration, new_ll),
            delta = ll - new_ll
            if delta < 0.0:
                delta = numpy.inf
            ll = new_ll
        if not self.gc_corrections is None:
            self.gc_corrections = numpy.log(self.gc_corrections)
        if not self.len_corrections is None:
            self.len_corrections = numpy.log(self.len_corrections)
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning regression corrections... Final ll:%f\n") % (' ' * 80, ll),
        self.normalization = 'regression'
        return None

    def find_trans_mean(self):
        """
        Calculate the mean signal across all valid fragment-pair trans (inter-region) interactions.
        
        :returns: None
        """
        if not self.silent:
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
        if not self.silent:
            print >> sys.stderr, ('Done\n'),
        return None

    def cis_heatmap(self, region, binsize=0, binbounds=None, start=None, stop=None, startfrag=None, stopfrag=None,
                    datatype='enrichment', arraytype='full', skipfiltered=False, returnmapping=False,
                    dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                    removefailed=False, image_file=None, **kwargs):
        """
        Return a heatmap of cis data of the type and shape specified by the passed arguments.

        This function returns a heatmap for a single region, bounded by either 'start' and 'stop' or 'startfend' and 'stopfend' ('start' and 'stop' take precedence). If neither is given, the complete region is included. The data in the array is determined by the 'datatype', being raw, fragment-corrected, distance-corrected, enrichment, or expected data. The array shape is given by 'arraytype' and can be compact (if unbinned), upper, or full. See :mod:`fivec_binning <hifive.fivec_binning>` for further explanation of 'datatype' and 'arraytype'. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param region: The index of the region to obtain data from.
        :type region: int.
        :param binsize: This is the coordinate width of each bin. If 'binsize' is zero, unbinned data is returned.
        :type binsize: int.
        :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fragment not falling in a bin is ignored.
        :type binbounds: numpy array
        :param start: The smallest coordinate to include in the array, measured from fragment midpoints. If both 'start' and 'startfrag' are given, 'start' will override 'startfrag'. If unspecified, this will be set to the midpoint of the first fragment for 'region'. Optional.
        :type start: int.
        :param stop: The largest coordinate to include in the array, measured from fragment midpoints. If both 'stop' and 'stopfrag' are given, 'stop' will override 'stopfrag'. If unspecified, this will be set to the midpoint of the last fragment plus one for 'region'. Optional.
        :type stop: int.
        :param startfrag: The first fragment to include in the array. If unspecified and 'start' is not given, this is set to the first fragment in 'region'. In cases where 'start' is specified and conflicts with 'startfrag', 'start' is given preference. Optional
        :type startfrag: int.
        :param stopfrag: The first fragment not to include in the array. If unspecified and 'stop' is not given, this is set to the last fragment in 'region' plus one. In cases where 'stop' is specified and conflicts with 'stopfrag', 'stop' is given preference. Optional.
        :type stopfrag: str.
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact' (if unbinned), 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N and M are the number of forward and reverse probe fragments, respectively. 'full' returns a square, symmetric array of size N x N x 2 where N is the total number of fragments or bins. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2, where N is the total number of fragments or bins.
        :type arraytype: str.
        :param skipfiltered: If True, all interaction bins for filtered out fragments are removed and a reduced-size array is returned.
        :type skipfiltered: bool.
        :param returnmapping: If True, a list containing the data array and either a 1d array containing fragment numbers included in the data array if the array is not compact or two 1d arrays containin fragment numbers for forward and reverse fragments if the array is compact is return. Otherwise only the data array is returned.
        :type returnmapping: bool.
        :param dynamically_binned: If True, return dynamically binned data.
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
        if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
            if not self.silent:
                print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        if ((binsize != 0 and arraytype not in ['full', 'upper']) or
            (arraytype not in ['full', 'compact', 'upper'])):
            if not self.silent:
                print >> sys.stderr, ("Unrecognized or inappropriate array type. No data returned.\n"),
            return None
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            data = fivec_binning.bin_cis_signal(self, region, binsize=binsize, binbounds=binbounds,
                                                start=start, stop=stop, startfrag=startfrag,
                                                stopfrag=stopfrag, datatype=datatype, arraytype=arraytype,
                                                skipfiltered=skipfiltered, returnmapping=returnmapping,
                                                silent=self.silent)
        else:
            expansion, exp_mapping = fivec_binning.find_cis_signal(self, region, binsize=expansion_binsize,
                                                                   start=start, stop=stop, startfrag=startfrag,
                                                                   stopfrag=stopfrag, datatype=datatype,
                                                                   arraytype=arraytype, skipfiltered=True,
                                                                   returnmapping=True, silent=self.silent)
            binned, mapping = fivec_binning.find_cis_signal(self, region, binsize=binsize, binbounds=binbounds,
                                                            start=start, stop=stop, startfrag=startfrag,
                                                            stopfrag=stopfrag, datatype=datatype, arraytype=arraytype,
                                                            returnmapping=True, silent=self.silent)
            fivec_binning.dynamically_bin_cis_array(expansion, exp_mapping, binned, mapping,
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
            if arraytype == 'upper':
                img = plotting.plot_upper_array(binned, silent=self.silent, **kwargs)
            else:
                img = plotting.plot_full_array(binned, silent=self.silent, **kwargs)
            img.save(image_file, format='png')
        return data

    def trans_heatmap(self, region1, region2, binsize=1000000, binbounds1=None, start1=None, stop1=None,
                      startfrag1=None, stopfrag1=None, binbounds2=None, start2=None, stop2=None, startfrag2=None,
                      stopfrag2=None, datatype='enrichment', arraytype='full', returnmapping=False,
                      dynamically_binned=False, minobservations=0,  searchdistance=0, expansion_binsize=0,
                      removefailed=False, skipfiltered=False, image_file=None, **kwargs):
        """
        Return a heatmap of trans data of the type and shape specified by the passed arguments.

        This function returns a heatmap for trans interactions between two regions, bounded by either 'start1', 'stop1', 'start2' and 'stop2' or 'startfrag1', 'stopfrag1', 'startfrag2', and 'stopfrag2' ('start' and 'stop' take precedence). The data in the array is determined by the 'datatype', being raw, fragment-corrected, distance-corrected, enrichment, or expected data. The array shape is always rectangular but can be either compact (which returns two arrays) or full. See :mod:`fivec_binning <hifive.fivec_binning>` for further explanation of 'datatype' and 'arraytype'. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param region1: The index of the first region to obtain data from.
        :type region1: int.
        :param region2: The index of the second region to obtain data from.
        :type region2: int.
        :param binsize: This is the coordinate width of each bin.
        :type binsize: int.
        :param binbounds1: An array containing start and stop coordinates for a set of user-defined bins for 'region1'. Any fragment not falling in a bin is ignored.
        :type binbounds1: numpy array
        :param start1: The coordinate at the beginning of the smallest bin from 'region1'. If unspecified, 'start1' will be the first multiple of 'binsize' below the 'startfrag1' mid. If there is a conflict between 'start1' and 'startfrag1', 'start1' is given preference. Optional.
        :type start1: int.
        :param stop1: The largest coordinate to include in the array from 'region1', measured from fragment midpoints. If both 'stop1' and 'stopfrag1' are given, 'stop1' will override 'stopfrag1'. 'stop1' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop1: int.
        :param startfrag1: The first fragment from 'region1' to include in the array. If unspecified and 'start1' is not given, this is set to the first valid fend in 'region1'. In cases where 'start1' is specified and conflicts with 'startfrag1', 'start1' is given preference. Optional.
        :type startfrag1: int.
        :param stopfrag1: The first fragment not to include in the array from 'region1'. If unspecified and 'stop1' is not given, this is set to the last valid fragment in 'region1' + 1. In cases where 'stop1' is specified and conflicts with 'stopfrag1', 'stop1' is given preference. Optional.
        :type stopfrag1: int.
        :param start1: The coordinate at the beginning of the smallest bin from 'region1'. If unspecified, 'start1' will be the first multiple of 'binsize' below the 'startfrag1' mid. If there is a conflict between 'start1' and 'startfrag1', 'start1' is given preference. Optional.
        :param binbounds2: An array containing start and stop coordinates for a set of user-defined bins for 'region2'. Any fragment not falling in a bin is ignored.
        :type binbounds2: numpy array
        :type start2: int.
        :param stop2: The largest coordinate to include in the array from 'region2', measured from fragment midpoints. If both 'stop2' and 'stopfrag2' are given, 'stop2' will override 'stopfrag2'. 'stop2' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop2: int.
        :param startfrag2: The first fragment from 'region2' to include in the array. If unspecified and 'start2' is not given, this is set to the first valid fend in 'region2'. In cases where 'start2' is specified and conflicts with 'startfrag2', 'start2' is given preference. Optional.
        :type startfrag2: int.
        :param stopfrag2: The first fragment not to include in the array from 'region2'. If unspecified and 'stop2' is not given, this is set to the last valid fragment in 'region2' + 2. In cases where 'stop2' is specified and conflicts with 'stopfrag2', 'stop2' is given preference. Optional.
        :type stopfrag2: int.
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment', 'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-filtered bins return value of 1. Expected values are returned for 'distance', 'fragment', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact' (if unbinned) and 'full'. 'compact' means data are arranged in a N x M x 2 array where N and M are the number of forward and reverse probe fragments, respectively. If compact is selected, only data for the forward primers of 'region1' and reverse primers of 'region2' are returned. 'full' returns a square, symmetric array of size N x N x 2 where N is the total number of fragments or bins.
        :type arraytype: str.
        :param returnmapping: If 'True', a list containing the data array and mapping information is returned. Otherwise only a data array(s) is returned.
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
        :param skipfiltered: If 'True', all interaction bins for filtered out fragments are removed and a reduced-size array is returned.
        :type skipfiltered: bool.
        :param image_file: If a filename is specified, a PNG image file is written containing the heatmap data. Arguments for the appearance of the image can be passed as additional keyword arguments.
        :type image_file: str.
        :returns: Array in format requested with 'arraytype' containing inter-region data requested with 'datatype'. If 'returnmapping' is True, a list is returned with mapping information. If 'arraytype' is 'full', a single data array and two 1d arrays of fragments corresponding to rows and columns, respectively is returned. If 'arraytype' is 'compact', two data arrays are returned (forward1 by reverse2 and forward2 by reverse1) along with forward and reverse fragment positions for each array for a total of 5 arrays.
        """
        # check that all values are acceptable
        if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
            if not self.silent:
                print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        if arraytype != 'full' or (not dynamically_binned and arraytype not in ['full', 'compact']):
            if not self.silent:
                print >> sys.stderr, ("Unrecognized or inappropriate array type. No data returned.\n"),
            return None
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            data = fivec_binning.find_trans_signal(self, region1, region2, binsize=binsize, binbounds1=binbounds1,
                                                   start1=start1, stop1=stop1, startfrag1=startfrag1,
                                                   stopfrag1=stopfrag1, binbounds2=binbounds2, start2=start2,
                                                   stop2=stop2, startfrag2=startfrag2, stopfrag2=stopfrag2,
                                                   datatype=datatype, arraytype=arraytype,
                                                   skipfiltered=skipfiltered, returnmapping=returnmapping,
                                                   silent=self.silent)
        else:
            expansion, exp_map1, exp_map2 = fivec_binning.find_trans_signal(self, region1, region2, start1=start1,
                                                                            stop1=stop1, startfrag1=startfrag1,
                                                                            stopfrag1=stopfrag1, start2=start2,
                                                                            stop2=stop2, startfrag2=startfrag2,
                                                                            stopfrag2=stopfrag2, datatype=datatype,
                                                                            arraytype='full', skipfiltered=True,
                                                                            expansion_binsize=expansion_binsize,
                                                                            returnmapping=True, silent=self.silent)
            binned, mapping1, mapping2 = fivec_binning.find_trans_signal(self, region1, region2, binsize=binsize,
                                                                         binbounds1=binbounds1, start1=start1,
                                                                         stop1=stop1, startfrag1=startfrag1,
                                                                         stopfrag1=stopfrag1, binbounds2=binbounds2,
                                                                         start2=start2, stop2=stop2,
                                                                         startfrag2=startfrag2, stopfrag2=stopfrag2,
                                                                         datatype=datatype, arraytype='full',
                                                                         skipfiltered=True, returnmapping=True,
                                                                         silent=self.silent)
            fivec_binning.dynamically_bin_trans_array(expansion, exp_map1, exp_map2, binned, mapping1, mapping2,
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

    def write_heatmap_dict(self, filename, binsize, includetrans=True, datatype='enrichment', arraytype='full',
                           regions=[]):
        """
        Create an h5dict file containing binned interaction arrays, bin positions, and an index of included regions.

        :param filename: Location to write h5dict object to.
        :type filename: str.
        :param binsize: Size of bins for interaction arrays. If "binsize" is zero, fragment interactions are returned without binning.
        :type binsize: int.
        :param includetrans: Indicates whether trans interaction arrays should be calculated and saved.
        :type includetrans: bool.
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment', 'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-filtered bins return value of 1. Expected values are returned for 'distance', 'fragment', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact' and 'full'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2.
        :type arraytype: str.
        :param regions: If given, indicates which regions should be included. If left empty, all regions are included.
        :type regions: list.
        :returns: None
        """
        if (regions is None or
                (isinstance(regions, list) and
                (len(regions) == 0 or
                (len(regions) == 1 and regions[0] == ''))) or
                regions == ''):
            regions = list(numpy.arange(self.frags['regions'].shape[0]))
        else:
            for i in range(len(regions)):
                regions[i] = int(regions[i])
        fivec_binning.write_heatmap_dict(self, filename, binsize, includetrans=includetrans,
                                         datatype=datatype, arraytype=arraytype,
                                         regions=regions, silent=self.silent)
        return None
