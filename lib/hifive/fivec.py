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

    def cis_heatmap(self, region, start=None, stop=None, startfrag=None, stopfrag=None, binsize=0,
                    datatype='enrichment', arraytype='full', skipfiltered=False, returnmapping=False,
                    dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                    removefailed=False):
        """
        Return a heatmap of cis data of the type and shape specified by the passed arguments.

        This function returns a heatmap for a single region, bounded by either 'start' and 'stop' or 'startfend' and 'stopfend' ('start' and 'stop' take precedence). If neither is given, the complete region is included. The data in the array is determined by the 'datatype', being raw, fragment-corrected, distance-corrected, enrichment, or expected data. The array shape is given by 'arraytype' and can be compact (if unbinned), upper, or full. See :mod:`fivec_binning <hifive.fivec_binning>` for further explanation of 'datatype' and 'arraytype'. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param region: The index of the region to obtain data from.
        :type region: int.
        :param start: The smallest coordinate to include in the array, measured from fragment midpoints. If both 'start' and 'startfrag' are given, 'start' will override 'startfrag'. If unspecified, this will be set to the midpoint of the first fragment for 'region'. Optional.
        :type start: int.
        :param stop: The largest coordinate to include in the array, measured from fragment midpoints. If both 'stop' and 'stopfrag' are given, 'stop' will override 'stopfrag'. If unspecified, this will be set to the midpoint of the last fragment plus one for 'region'. Optional.
        :type stop: int.
        :param startfrag: The first fragment to include in the array. If unspecified and 'start' is not given, this is set to the first fragment in 'region'. In cases where 'start' is specified and conflicts with 'startfrag', 'start' is given preference. Optional
        :type startfrag: int.
        :param stopfrag: The first fragment not to include in the array. If unspecified and 'stop' is not given, this is set to the last fragment in 'region' plus one. In cases where 'stop' is specified and conflicts with 'stopfrag', 'stop' is given preference. Optional.
        :type stopfrag: str.
        :param binsize: This is the coordinate width of each bin. If 'binsize' is zero, unbinned data is returned.
        :type binsize: int.
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
        :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'. If returnmapping is True, a list is returned containined the requested data array and an array of associated positions (dependent on the binning options selected).
        """
        # check that all values are acceptable
        datatypes = {'raw': 0, 'fragment': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
        if datatype not in datatypes:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        else:
            datatype_int = datatypes[datatype]
        if ((dynamically_binned and arraytype not in ['full', 'upper']) or
            (binsize != 0 and arraytype not in ['full', 'upper']) or
            (arraytype not in ['full', 'compact', 'upper'])):
            print >> sys.stderr, ("Unrecognized or inappropriate array type. No data returned.\n"),
            return None
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            # determine if data is to be unbinned or binned
            if binsize == 0:
                # data should be unbinned
                data = fivec_binning.unbinned_cis_signal(self, region, start=start, stop=stop, startfrag=startfrag,
                                                         stopfrag=stopfrag, datatype=datatype, arraytype=arraytype,
                                                         skipfiltered=skipfiltered, returnmapping=returnmapping)
            else:
                # data should be binned
                data = fivec_binning.bin_cis_signal(self, region, start=start, stop=stop, startfrag=startfrag,
                                                    stopfrag=stopfrag, binsize=binsize, datatype=datatype,
                                                    arraytype=arraytype, returnmapping=returnmapping)
        else:
            if expansion_binsize == 0:
                # data should be dynamically binned with unbinned expansion data
                expansion, frags = fivec_binning.unbinned_cis_signal(self, region, start=start, stop=stop,
                                                                     startfrag=startfrag, stopfrag=stopfrag,
                                                                     datatype=datatype, arraytype=arraytype,
                                                                     skipfiltered=True, returnmapping=True)
                mids = self.frags['fragments']['mid'][frags]
            else:
                expansion, mapping = fivec_binning.bin_cis_signal(self, region, start=start, stop=stop,
                                                                  startfrag=startfrag, stopfrag=stopfrag,
                                                                  binsize=expansion_binsize, datatype=datatype,
                                                                  arraytype=arraytype, returnmapping=True)
                mids = (mapping[:, 2] + mapping[:, 3]) / 2
            binned, mapping = fivec_binning.bin_cis_signal(self, region, start=start, stop=stop, startfrag=startfrag,
                                                           stopfrag=stopfrag, binsize=binsize, datatype=datatype,
                                                           arraytype=arraytype, returnmapping=True)
            fivec_binning.dynamically_bin_cis_array(expansion, mids, binned, mapping[:, 2:],
                                                    minobservations=minobservations, searchdistance=searchdistance,
                                                    removefailed=removefailed)
            if returnmapping:
                data = [binned, mapping]
            else:
                data = binned
        return data

    def trans_heatmap(self, region1, region2, start1=0, stop1=None, startfrag1=None, stopfrag1=None, start2=0,
                      stop2=None, startfrag2=None, stopfrag2=None, binsize=1000000, datatype='enrichment',
                      arraytype='full', returnmapping=False, dynamically_binned=False, minobservations=0,
                      searchdistance=0, expansion_binsize=0, removefailed=False):
        """
        Return a heatmap of trans data of the type and shape specified by the passed arguments.

        This function returns a heatmap for trans interactions between two regions, bounded by either 'start1', 'stop1', 'start2' and 'stop2' or 'startfrag1', 'stopfrag1', 'startfrag2', and 'stopfrag2' ('start' and 'stop' take precedence). The data in the array is determined by the 'datatype', being raw, fragment-corrected, distance-corrected, enrichment, or expected data. The array shape is always rectangular but can be either compact (which returns two arrays) or full. See :mod:`fivec_binning <hifive.fivec_binning>` for further explanation of 'datatype' and 'arraytype'. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

        :param region1: The index of the first region to obtain data from.
        :type region1: int.
        :param region2: The index of the second region to obtain data from.
        :type region2: int.
        :param start1: The coordinate at the beginning of the smallest bin from 'region1'. If unspecified, 'start1' will be the first multiple of 'binsize' below the 'startfrag1' mid. If there is a conflict between 'start1' and 'startfrag1', 'start1' is given preference. Optional.
        :type start1: int.
        :param stop1: The largest coordinate to include in the array from 'region1', measured from fragment midpoints. If both 'stop1' and 'stopfrag1' are given, 'stop1' will override 'stopfrag1'. 'stop1' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop1: int.
        :param startfrag1: The first fragment from 'region1' to include in the array. If unspecified and 'start1' is not given, this is set to the first valid fend in 'region1'. In cases where 'start1' is specified and conflicts with 'startfrag1', 'start1' is given preference. Optional.
        :type startfrag1: int.
        :param stopfrag1: The first fragment not to include in the array from 'region1'. If unspecified and 'stop1' is not given, this is set to the last valid fragment in 'region1' + 1. In cases where 'stop1' is specified and conflicts with 'stopfrag1', 'stop1' is given preference. Optional.
        :type stopfrag1: int.
        :param start1: The coordinate at the beginning of the smallest bin from 'region1'. If unspecified, 'start1' will be the first multiple of 'binsize' below the 'startfrag1' mid. If there is a conflict between 'start1' and 'startfrag1', 'start1' is given preference. Optional.
        :type start2: int.
        :param stop2: The largest coordinate to include in the array from 'region2', measured from fragment midpoints. If both 'stop2' and 'stopfrag2' are given, 'stop2' will override 'stopfrag2'. 'stop2' will be shifted higher as needed to make the last bin of size 'binsize'. Optional.
        :type stop2: int.
        :param startfrag2: The first fragment from 'region2' to include in the array. If unspecified and 'start2' is not given, this is set to the first valid fend in 'region2'. In cases where 'start2' is specified and conflicts with 'startfrag2', 'start2' is given preference. Optional.
        :type startfrag2: int.
        :param stopfrag2: The first fragment not to include in the array from 'region2'. If unspecified and 'stop2' is not given, this is set to the last valid fragment in 'region2' + 2. In cases where 'stop2' is specified and conflicts with 'stopfrag2', 'stop2' is given preference. Optional.
        :type stopfrag2: int.
        :param binsize: This is the coordinate width of each bin.
        :type binsize: int.
        :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment', 'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-filtered bins return value of 1. Expected values are returned for 'distance', 'fragment', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
        :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact' (if unbinned) and 'full'. 'compact' means data are arranged in a N x M x 2 array where N and M are the number of forward and reverse probe fragments, respectively. Two arrays will be returned for this format, the first with forward probe fragments from region1 and reverse probe fragments from region2. The second is the compliment of the first. 'full' returns a square, symmetric array of size N x N x 2 where N is the total number of fragments or bins.
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
        :returns: Array in format requested with 'arraytype' containing inter-region data requested with 'datatype'. If 'returnmapping' is True, a list is returned with mapping information. If 'arraytype' is 'full', a single data array and two 1d arrays of fragments corresponding to rows and columns, respectively is returned. If 'arraytype' is 'compact', two data arrays are returned (forward1 by reverse2 and forward2 by reverse1) along with forward and reverse fragment positions for each array for a total of 5 arrays.
        """
        # check that all values are acceptable
        datatypes = {'raw': 0, 'fragment': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
        if datatype not in datatypes:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        else:
            datatype_int = datatypes[datatype]
        if ((dynamically_binned and arraytype != 'full') or
            (binsize != 0 and arraytype != 'full') or
            (arraytype not in ['full', 'compact'])):
            print >> sys.stderr, ("Unrecognized or inappropriate array type. No data returned.\n"),
            return None
        # determine if data is to be dynamically binned
        if not dynamically_binned:
            if binsize == 0:
                data = fivec_binning.unbinned_trans_signal(self, region1, region2, start1=start1, stop1=stop1,
                                                           startfrag1=startfrag1, stopfrag1=stopfrag1, start2=start2,
                                                           stop2=stop2, startfrag2=startfrag2, stopfrag2=stopfrag2,
                                                           datatype=datatype, arraytype=arraytype,
                                                           skipfiltered=skipfiltered, returnmapping=returnmapping)
            else:
                # data should be binned
                data = fivec_binning.bin_trans_signal(self, region1, region2, start1=start1, stop1=stop1,
                                                      startfrag1=startfrag1, stopfrag1=stopfrag1, start2=start2,
                                                      stop2=stop2, startfrag2=startfrag2, stopfrag2=stopfrag2,
                                                      binsize=binsize, datatype=datatype, returnmapping=returnmapping)
        else:
            if expansion_binsize == 0:
                expansion, frags1, frags2 = fivec_binning.unbinned_trans_signal(self, region1, region2, start1=start1,
                                                                                stop1=stop1, startfrag1=startfrag1,
                                                                                stopfrag1=stopfrag1, start2=start2,
                                                                                stop2=stop2, startfrag2=startfrag2,
                                                                                stopfrag2=stopfrag2, datatype=datatype,
                                                                                arraytype='full', skipfiltered=True,
                                                                                returnmapping=True)
                mids1 = self.frags['fragments']['mid'][frags1]
                mids2 = self.frags['fragments']['mid'][frags2]
            else:
                expansion, map1, map2 = fivec_binning.bin_trans_signal(self, region1, region2, start1=start1,
                                                                       stop1=stop1, startfrag1=startfrag1,
                                                                       stopfrag1=stopfrag1, start2=start2, stop2=stop2,
                                                                       startfrag2=startfrag2, stopfrag2=stopfrag2,
                                                                       binsize=expansion_binsize, datatype=datatype,
                                                                       returnmapping=True)
                mids1 = (map1[:, 2] + map1[:, 3]) / 2
                mids2 = (map2[:, 2] + map2[:, 3]) / 2
            if binsize == 0:
                binned, mapping1, mapping2 = fivec_binning.unbinned_trans_signal(self, region1, region2, start1=start1,
                                                                                 stop1=stop1, startfrag1=startfrag1,
                                                                                 stopfrag1=stopfrag1, start2=start2,
                                                                                 stop2=stop2, startfrag2=startfrag2,
                                                                                 stopfrag2=stopfrag2,
                                                                                 datatype=datatype,
                                                                                 arraytype='full', skipfiltered=True,
                                                                                 returnmapping=True)
                bounds1 = numpy.hstack((self.frags['fragments']['start'][mapping1].reshape(-1, 1),
                                        self.frags['fragments']['stop'][mapping1].reshape(-1, 1)))
                bounds2 = numpy.hstack((self.frags['fragments']['start'][mapping2].reshape(-1, 1),
                                        self.frags['fragments']['stop'][mapping2].reshape(-1, 1)))
            else:
                binned, mapping1, mapping2 = fivec_binning.bin_trans_signal(self, region1, region2, start1=start1,
                                                                            stop1=stop1, startfrag1=startfrag1,
                                                                            stopfrag1=stopfrag1, start2=start2,
                                                                            stop2=stop2, startfrag2=startfrag2,
                                                                            stopfrag2=stopfrag2,
                                                                            binsize=binsize, datatype=datatype,
                                                                            returnmapping=True)
                bounds1 = mapping1[:, 2:4]
                bounds2 = mapping2[:, 2:4]
            fivec_binning.dynamically_bin_trans_array(expansion, mids1, mids2, binned, bounds1, bounds2,
                                                      minobservations=minobservations, searchdistance=searchdistance,
                                                      removefailed=removefailed)
            if returnmapping:
                data = [binned, mapping1, mapping2]
            else:
                data = binned
        return data

    def write_heatmap_dict(self, filename, binsize, includetrans=True, remove_distance=False, arraytype='full',
                           regions=[]):
        """
        Create an h5dict file containing binned interaction arrays, bin positions, and an index of included regions.

        :param filename: Location to write h5dict object to.
        :type filename: str.
        :param binsize: Size of bins for interaction arrays. If "binsize" is zero, fragment interactions are returned without binning.
        :type binsize: int.
        :param includetrans: Indicates whether trans interaction arrays should be calculated and saved.
        :type includetrans: bool.
        :param remove_distance: If 'True', the expected value is calculated including the expected distance mean. Otherwise, only fragment corrections are used.
        :type remove_distance: bool.
        :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact' and 'full'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2.
        :type arraytype: str.
        :param regions: If given, indicates which regions should be included. If left empty, all regions are included.
        :type regions: list.
        :returns: None
        """
        fivec_binning.write_heatmap_dict(self, filename, binsize, includetrans=includetrans,
                                         remove_distance=remove_distance, arraytype=arraytype,
                                         regions=regions)
        return None
