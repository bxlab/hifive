#!/usr/bin/env python

"""A class for handling 5C analysis."""

import os
import sys
from math import log

import numpy
from scipy.stats import linregress
import h5py
from scipy.optimize import fmin_l_bfgs_b as bfgs

import libraries._fivec_binning as _binning
import libraries._fivec_optimize as _optimize
import fivec_binning
import plotting


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

    :attributes: * **file** (*str.*) - A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcome.
                 * **normalization** (*str.*) - A string stating which type of normalization has been performed on this object. This starts with the value 'none'.

    In addition, many other attributes are initialized to the 'None' state.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a FiveC object."""
        self.file = os.path.abspath(filename)
        self.filetype = 'fivec_project'
        self.silent = silent        
        self.binning_corrections = None
        self.binning_correction_indices = None
        self.binning_frag_indices = None
        self.binning_num_bins = None
        self.model_parameters = None
        self.corrections = None
        self.region_means = None
        self.gamma = None
        self.sigma = None
        self.trans_mean = None
        self.normalization = 'none'
        self.history = ''
        if mode != 'w':
            self.load()
        return None

    def __getitem__(self, key):
        """Dictionary-like lookup."""
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        """Dictionary-like value setting."""
        self.__dict__[key] = value
        return None

    def load_data(self, filename):
        """
        Load fragment-pair counts and fragment object from :class:`FiveCData <hifive.fivec_data.FiveCData>` object.

        :param filename: Specifies the file name of the :class:`FiveCData <hifive.fivec_data.FiveCData>` object to associate with this analysis.
        :type filename: str.
        :returns: None

        :Attributes: * **datafilename** (*str.*) - A string containing the relative path of the FiveCData file.
                     * **fragfilename** (*str.*) - A string containing the relative path of the Fragment file associated with the FiveCData file.
                     * **frags** (*filestream*) - A filestream to the hdf5 Fragment file such that all saved Fragment attributes can be accessed through this class attribute.
                     * **data** (*filestream*) - A filestream to the hdf5 FiveCData file such that all saved FiveCData attributes can be accessed through this class attribute.
                     * **chr2int** (*dict.*) - A dictionary that converts chromosome names to chromosome indices.
                     * **filter** (*ndarray*) - A numpy array of type int32 and size N where N is the number of fragments. This contains the inclusion status of each fragment with a one indicating included and zero indicating excluded and is initialized with all fragments included.

        When a FiveCData object is associated with the project file, the 'history' attribute is updated with the history of the FiveCData object.
        """
        self.history += "FiveC.load_data(filename='%s') - " % filename
        # ensure data h5dict exists
        if not os.path.exists(filename):
            if not self.silent:
                print >> sys.stderr, ("Could not find %s. No data loaded.\n") % (filename.split('/')[-1]),
            self.history += "Error: '%s' not found\n" % filename
            return None
        self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                       os.path.dirname(self.file)), os.path.basename(filename))
        self.data = h5py.File(filename, 'r')
        self.history = self.data['/'].attrs['history'] + self.history
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
            self.history += "Error: '%s' not found\n" % fragfilename
            return None
        self.frags = h5py.File(fragfilename, 'r')
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.frags['chromosomes']):
            self.chr2int[chrom] = i
        # create arrays
        self.filter = numpy.ones(self.frags['fragments'].shape[0], dtype=numpy.int32)
        self.history += 'Success\n'
        return None

    def save(self, out_fname=None):
        """
        Save analysis parameters to h5dict.

        :param filename: Specifies the file name of the :class:`FiveC <hifive.fivec.FiveC>` object to save this analysis to.
        :type filename: str.
        :returns: None
        """ 
        self.history.replace("'None'", "None")
        if not out_fname is None:
            original_file = os.path.abspath(self.file)
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
        else:
            out_fname = self.file
        datafile = h5py.File(out_fname, 'w')
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
        # return attributes to init state
        self.binning_corrections = None
        self.binning_correction_indices = None
        self.binning_frag_indices = None
        self.binning_num_bins = None
        self.model_parameters = None
        self.corrections = None
        self.region_means = None
        self.gamma = None
        self.sigma = None
        self.trans_mean = None
        self.normalization = 'none'
        self.history = ''
        # load data hdf5 dict
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
                self.frags = h5py.File(fragfilename, 'r')
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
        self.history += "FiveC.filter_fragments(mininteractions=%i, mindistance=%s, maxdistance=%s) - " % (mininteractions, str(mindistance), str(maxdistance))
        if not self.silent:
            print >> sys.stderr, ("Filtering fragments..."),
        original_count = numpy.sum(self.filter)
        previous_valid = original_count + 1
        current_valid = original_count
        coverage = numpy.zeros(self.filter.shape[0], dtype=numpy.int32)
         # copy needed arrays
        data = self.data['cis_data'][...]
        distances = self.frags['fragments']['mid'][data[:, 1]] - self.frags['fragments']['mid'][data[:, 0]]
        if maxdistance == 0 or maxdistance is None:
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
        self.history += "Success\n"
        return None

    def find_distance_parameters(self):
        """
        Regress log counts versus inter-fragment distances to find slope and intercept values and then find the standard deviation of corrected counts.

        :returns: None

        :Attributes: * **gamma** (*float*) - A float denoting the negative slope of the distance-dependence regression line.
                     * **sigma** (*float*) - A float denoting the standard deviation of nonzero data about the distance-dependence regression line.
                     * **region_means** (*ndarray*) - A numpy array of type float32 and length equal to the number of regions. This is initialized to zeros until fragment correction values are found.
        """
        self.history += "FiveC.find_distance_parameters() - "
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
        counts = numpy.log(data[valid, 2])
        if not self.corrections is None:
            counts -= self.corrections[data[valid, 0]] - self.corrections[data[valid, 1]]
        temp = linregress(log_distances, counts)[:2]
        self.gamma = -float(temp[0])
        if self.region_means is None:
            self.region_means = numpy.zeros(self.frags['regions'].shape[0], dtype=numpy.float32) + temp[1]
        self.sigma = float(numpy.std(counts - temp[1] + self.gamma * log_distances))
        if not self.silent:
            print >> sys.stderr, ("Done\n"),
        self.history += "Success\n"
        return None

    def find_probability_fragment_corrections(self, mindistance=0, maxdistance=0, max_iterations=1000,
                                              minchange=0.0005, learningstep=0.1, precalculate=True, regions=[],
                                              precorrect=False):
        """
         Using gradient descent, learn correction values for each valid fragment based on a Log-Normal distribution of observations.

        :param mindistance: The minimum inter-fragment distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fragment distance to be included in modeling.
        :type maxdistance: int.
        :param max_iterations: The maximum number of iterations to carry on gradient descent for.
        :type max_iterations: int.
        :type annealing_iterations: int.
        :param minchange: The cutoff threshold for early learning termination for the maximum absolute gradient value.
        :type minchange: float
        :param learningstep: The scaling factor for decreasing learning rate by if step doesn't meet armijo criterion.
        :type learningstep: float
        :param precalculate: Specifies whether the correction values should be initialized at the fragment means.
        :type precalculate: bool.
        :param regions: A list of regions to calculate corrections for. If set as None, all region corrections are found.
        :type regions: list
        :param precorrect: Use binning-based corrections in expected value calculations, resulting in a chained normalization approach.
        :type precorrect: bool.
        :returns: None

        :Attributes: * **corrections** (*ndarray*) - A numpy array of type float32 and length equal to the number of fragments. All invalid fragments have an associated correction value of zero.

        The 'normalization' attribute is updated to 'probability' or 'binning-probability', depending on if the 'precorrect' option is selected. In addition, the 'region_means' attribute is updated such that the mean correction (sum of all valid regional correction value pairs) is adjusted to zero and the corresponding region mean is adjusted the same amount but the opposite sign. 
        """
        self.history += "FiveC.find_probability_fragment_corrections(mindistance=%s, maxdistance=%s, max_iterations=%i, minchange=%f, learningstep=%f, precalculate=%s, regions=%s, precorrect=%s) - " % (str(mindistance), str(maxdistance), max_iterations, minchange, learningstep, precalculate, str(regions), precorrect)
        if precorrect and self.binning_corrections is None:
            if not self.silent:
                print >> sys.stderr, ("Precorrection can only be used in project has previously run 'find_binning_fragment_corrections'.\n"),
            self.history += "Error: 'find_binning_fragment_corrections()' not run yet\n"
            return None
        if self.corrections is None:
            self.corrections = numpy.zeros(self.frags['fragments'].shape[0], dtype=numpy.float32)
        # if regions not given, set to all regions
        if regions == None or len(regions) == 0:
            regions = numpy.arange(self.frags['regions'].shape[0])
        # determine if distance parameters have been calculated
        if self.gamma is None:
            self.find_distance_parameters()
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
        if maxdistance == 0 or maxdistance is None:
            maxdistance = numpy.amax(distances) + 1
        valid = numpy.where((filt[data[:, 0]] * filt[data[:, 1]]) *
                            (distances >= mindistance) * (distances < maxdistance))[0]
        data = data[valid, :]
        distances = numpy.log(distances[valid])
        counts_n = numpy.log(data[:, 2] - 0.5).astype(numpy.float32)
        counts = numpy.log(data[:, 2]).astype(numpy.float32)
        counts_p = numpy.log(data[:, 2] + 0.5).astype(numpy.float32)
        distance_signal = (-self.gamma * distances).astype(numpy.float32)
        distance_signal += self.region_means[self.frags['fragments']['region'][data[:, 0]]]
        # create empty arrays
        gradients = numpy.zeros(self.filter.shape[0], dtype=numpy.float32)
        valid = numpy.where(filt)[0]
        # find number of interactions for each fragment
        interactions = numpy.bincount(data[:, 0], minlength=self.filter.shape[0]).astype(numpy.int32)
        interactions += numpy.bincount(data[:, 1], minlength=self.filter.shape[0]).astype(numpy.int32)
        interactions = numpy.maximum(1, interactions)
        # if precalculation requested, find fragment means
        if precalculate:
            enrichments = counts - distance_signal
            count_sums = numpy.bincount(data[:, 0], weights=enrichments, minlength=gradients.shape[0])
            count_sums += numpy.bincount(data[:, 1], weights=enrichments, minlength=gradients.shape[0])
            self.corrections = ((count_sums / numpy.maximum(1, interactions)) * 0.5).astype(numpy.float32)
        if precorrect:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding binning corrections...") % (' ' * 80),
            _optimize.find_binning_correction_adjustment(distance_signal,
                                                            data,
                                                            self.binning_corrections,
                                                            self.binning_correction_indices,
                                                            self.binning_num_bins,
                                                            self.binning_frag_indices)
        # cycle through learning phases
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning corrections...") % (' ' * 80),
        iteration = 0
        cont = True
        change = numpy.inf
        new_corrections = numpy.copy(self.corrections)
        start_cost = _optimize.calculate_prob_cost(data,
                                                   counts_n,
                                                   counts,
                                                   counts_p,
                                                   distance_signal,
                                                   self.corrections,
                                                   self.sigma)
        previous_cost = start_cost
        while cont:
            iteration += 1
            # find gradients
            gradients.fill(0.0)
            _optimize.calculate_gradients(data,
                                          counts_n,
                                          counts,
                                          counts_p,
                                          distance_signal,
                                          self.corrections,
                                          gradients,
                                          self.sigma)
            # find best step size
            armijo = numpy.inf
            t = 0.1
            gradients /= interactions
            gradient_norm = numpy.sum(gradients[valid] ** 2.0)
            j = 0
            best_score = numpy.inf
            best_t = 0.1
            while armijo > 0.0:
                # update gradients
                _optimize.update_corrections(filt,
                                             self.corrections,
                                             new_corrections,
                                             gradients,
                                             t)
                cost = _optimize.calculate_prob_cost(data,
                                                     counts_n,
                                                     counts,
                                                     counts_p,
                                                     distance_signal,
                                                     new_corrections,
                                                     self.sigma)
                if numpy.isnan(cost):
                    cost = numpy.inf
                    armijo = numpy.inf
                else:
                    armijo = cost - previous_cost + t * gradient_norm
                    if cost < best_score:
                        best_score = cost
                        best_t = t
                if not self.silent:
                    print >> sys.stderr, ("\r%s iteration:%i cost:%f change:%f armijo: %f %s") %\
                                         ('Learning corrections...', iteration, previous_cost,
                                          change, armijo, ' ' * 20),
                t *= learningstep
                j += 1
                if j == 20:
                    armijo = -numpy.inf
                    t = best_t
                    _optimize.update_corrections(filt,
                                                 self.corrections,
                                                 new_corrections,
                                                 gradients,
                                                 t)
                    cost = _optimize.calculate_prob_cost(data,
                                                         counts_n,
                                                         counts,
                                                         counts_p,
                                                         distance_signal,
                                                         new_corrections,
                                                         self.sigma)
            previous_cost = cost
            self.corrections = new_corrections
            change = numpy.amax(numpy.abs(gradients[valid] / new_corrections[valid]))
            if not self.silent:
                print >> sys.stderr, ("\r%s iteration:%i cost:%f change:%f %s") %\
                                     ('Learning corrections...', iteration, cost, change, ' ' * 40),
            iteration += 1
            if iteration >= max_iterations or change <= minchange:
                cont = False
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning corrections... Initial Cost: %f  Final Cost: %f  Done\n") %\
                                 (' ' * 80, start_cost, cost),
        # Calculate region means
        if self.region_means is None:
            self.region_means = numpy.zeros(self.frags['regions'].shape[0], dtype=numpy.float32)
        for i in regions:
            start = self.frags['regions']['start_frag'][i]
            stop = self.frags['regions']['stop_frag'][i]
            forward = (numpy.where(self.filter[start:stop] * 
                                   (self.frags['fragments']['strand'][start:stop] == 0))[0] + start)
            reverse = (numpy.where(self.filter[start:stop] * 
                                   (self.frags['fragments']['strand'][start:stop] == 1))[0] + start)
            if forward.shape[0] == 0 or reverse.shape[0] == 0:
                continue
            region_mean = (numpy.sum(self.corrections[forward]) * reverse.shape[0] +
                           numpy.sum(self.corrections[reverse]) * forward.shape[0])
            region_mean /= forward.shape[0] * reverse.shape[0]
            self.corrections[forward] -= region_mean / 2.0
            self.corrections[reverse] -= region_mean / 2.0
            self.region_means[i] += region_mean
        if precorrect:
            self.normalization = 'binning-probability'
        else:
            self.normalization = 'probability'
        self.history += 'Succcess\n'
        return None

    def find_express_fragment_corrections(self, mindistance=0, maxdistance=0, iterations=1000, remove_distance=False,
                                          usereads='cis', regions=[], precorrect=False, logged=True, kr=False):
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
        :param precorrect: Use binning-based corrections in expected value calculations, resulting in a chained normalization approach.
        :type precorrect: bool.
        :param logged: Use log-counts instead of counts for learning.
        :type logged: bool.
        :param kr: Use the Knight Ruiz matrix balancing algorithm instead of weighted matrix balancing. This option ignores 'iterations' and 'logged'.
        :type kr: bool.
        :returns: None

        Calling this function creates the following attributes:

        :Attributes: * **corrections** (*ndarray*) - A numpy array of type float32 and length equal to the number of fragments. All invalid fragments have an associated correction value of zero.

        The 'normalization' attribute is updated to 'express' or 'binning-express', depending on if the 'precorrect' option is selected. In addition, if the 'remove_distance' option is selected, the 'region_means' attribute is updated such that the mean correction (sum of all valid regional correction value pairs) is adjusted to zero and the corresponding region mean is adjusted the same amount but the opposite sign. 
        """
        self.history += "FiveC.find_express_fragment_corrections(mindistance=%s, maxdistance=%s, iterations=%i, remove_distance=%s, usereads='%s', regions=%s, precorrect=%s, logged=%s, kr=%s) - " % (str(mindistance), str(maxdistance), iterations, remove_distance, usereads, str(regions), precorrect, logged, kr)
        if precorrect and self.binning_corrections is None:
            if not self.silent:
                print >> sys.stderr, ("Precorrection can only be used in project has previously run 'find_binning_fragment_corrections'.\n"),
            self.history += "Error: 'find_binning_fragment_corrections()' not run yet\n"
            return None
        # make sure usereads has a valid value
        if usereads not in ['cis', 'trans', 'all']:
            if not self.silent:
                print >> sys.stderr, ("'usereads' does not have a valid value.\n"),
            self.history += "Error: '%s' not a valid value for 'usereads'\n" % usereads
            return None
        # if regions not given, set to all regions
        if regions == None or len(regions) == 0:
            regions = numpy.arange(self.frags['regions'].shape[0])
        if self.corrections is None:
            self.corrections = numpy.zeros(self.frags['fragments'].shape[0], dtype=numpy.float32)
        if kr:
            self._find_kr_corrections(mindistance, maxdistance, remove_distance, 
                             usereads, regions, precorrect, logged)
            return None
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
        counts = None
        trans_counts = None
        distance_signal = None
        trans_signal = None
        corrections = numpy.copy(self.corrections)
        if usereads in ['cis', 'all']:
            data = self.data['cis_data'][...]
            distances = (self.frags['fragments']['mid'][data[:, 1]] -
                         self.frags['fragments']['mid'][data[:, 0]]).astype(numpy.float32)
            if maxdistance == 0 or maxdistance is None:
                maxdistance = numpy.amax(distances) + 1
            valid = numpy.where((filt[data[:, 0]] * filt[data[:, 1]]) *
                                (distances >= mindistance) * (distances < maxdistance))[0]
            data = data[valid, :]
            counts = numpy.log(data[:, 2]).astype(numpy.float64)
            distances = distances[valid]
            if remove_distance:
                if self.gamma is None:
                    self.find_distance_parameters()
                distance_signal = (-self.gamma * numpy.log(distances)).astype(numpy.float32)
                distance_signal += self.region_means[self.frags['fragments']['region'][data[:, 0]]]
        if usereads in ['trans', 'all']:
            trans_data =  self.data['trans_data'][...]
            valid = numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0]
            trans_data = trans_data[valid, :]
            trans_counts = numpy.log(trans_data[:, 2]).astype(numpy.float64)
            if remove_distance:
                if self.trans_mean is None:
                    self.find_trans_mean()
                trans_signal = numpy.zeros(trans_data.shape[0], dtype=numpy.float32) + self.trans_mean
        if precorrect:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding binning corrections...") % (' ' * 80),
            if usereads in ['cis', 'all']:
                if distance_signal is None:
                    distance_signal = numpy.zeros(data.shape[0], dtype=numpy.float32)
                _optimize.find_binning_correction_adjustment(distance_signal,
                                                             data,
                                                             self.binning_corrections,
                                                             self.binning_correction_indices,
                                                             self.binning_num_bins,
                                                             self.binning_frag_indices)
            if usereads in ['trans', 'all']:
                if trans_signal is None:
                    trans_signal = numpy.zeros(trans_data.shape[0], dtype=numpy.float32)
                _optimize.find_binning_correction_adjustment(trans_signal,
                                                             trans_data,
                                                             self.binning_corrections,
                                                             self.binning_correction_indices,
                                                             self.binning_num_bins,
                                                             self.binning_frag_indices)
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
            if logged:
                cost = _optimize.find_log_fragment_means(distance_signal,
                                                     trans_signal,
                                                     interactions,
                                                     fragment_means,
                                                     data,
                                                     trans_data,
                                                     counts,
                                                     trans_counts,
                                                     corrections)
            else:
                cost = _optimize.find_fragment_means(distance_signal,
                                                     trans_signal,
                                                     interactions,
                                                     fragment_means,
                                                     data,
                                                     trans_data,
                                                     counts,
                                                     trans_counts,
                                                     corrections)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLearning corrections... iteration:%i  cost:%f ") % (' ' * 80, iteration,
                                                                                                 cost),
        where = numpy.where(filt)[0]
        self.corrections[where] = corrections[where]
        # Calculate region means
        if self.region_means is None:
            self.region_means = numpy.zeros(self.frags['regions'].shape[0], dtype=numpy.float32)
        for i in regions:
            start = self.frags['regions']['start_frag'][i]
            stop = self.frags['regions']['stop_frag'][i]
            forward = (numpy.where(self.filter[start:stop] * 
                                   (self.frags['fragments']['strand'][start:stop] == 0))[0] + start)
            reverse = (numpy.where(self.filter[start:stop] * 
                                   (self.frags['fragments']['strand'][start:stop] == 1))[0] + start)
            if forward.shape[0] == 0 or reverse.shape[0] == 0:
                continue
            region_mean = (numpy.sum(self.corrections[forward]) * reverse.shape[0] +
                           numpy.sum(self.corrections[reverse]) * forward.shape[0])
            region_mean /= forward.shape[0] * reverse.shape[0]
            self.corrections[forward] -= region_mean / 2.0
            self.corrections[reverse] -= region_mean / 2.0
            if remove_distance:
                self.region_means[i] += region_mean
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning corrections... Final Cost: %f  Done\n") % (' ' * 80, cost),
        if precorrect:
            self.normalization = 'binning-express'
        else:
            self.normalization = 'express'
        self.history += 'Success\n'
        return None

    def _find_kr_corrections(self, mindistance=0, maxdistance=0, remove_distance=True, 
                             usereads='cis', regions=[], precorrect=False, logged=False):
        if self.gamma is None:
            self.find_distance_parameters()
        all_regions = numpy.copy(regions)
        filt = numpy.copy(self.filter)
        if maxdistance == 0 or maxdistance is None:
            maxdistance = 99999999999
        if usereads != 'cis':
            for i in range(self.frags['regions'].shape[0]):
                if i not in regions:
                    filt[self.frags['regions']['start_frag'][i]:self.frags['regions']['stop_frag'][i]] = 0
            regions = ['all']
        for region in regions:
            if region == 'all':
                startfrag = 0
                stopfrag = self.frags['chr_indices'][-1]
                regfilt = filt
            else:
                startfrag = self.frags['regions']['start_frag'][region]
                stopfrag = self.frags['regions']['stop_frag'][region]
                regfilt = filt[startfrag:stopfrag]
            # create needed arrays
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLoading needed data...") % (' ' * 80),
            mids = self.frags['fragments']['mid'][startfrag:stopfrag]
            strands = self.frags['fragments']['strand'][startfrag:stopfrag]
            if usereads in ['cis', 'all']:
                start_index = self.data['cis_indices'][startfrag]
                stop_index = self.data['cis_indices'][stopfrag]
                data = self.data['cis_data'][start_index:stop_index, :]
                distances = mids[data[:, 1] - startfrag] - mids[data[:, 0] - startfrag]
                valid = numpy.where(regfilt[data[:, 0] - startfrag] * regfilt[data[:, 1] - startfrag] *
                                    (distances >= mindistance) * (distances < maxdistance))[0]
                data = data[valid, :]
            else:
                data = None
            if usereads in ['trans', 'all']:
                trans_data = self.data['trans_data'][...]
                valid = numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0]
                trans_data = trans_data[valid, :]
            else:
                trans_data = None
            trans_means = None
            # remapped data
            rev_mapping = numpy.where(regfilt)[0]
            mapping = numpy.zeros(regfilt.shape[0], dtype=numpy.int32) - 1
            mapping[rev_mapping] = numpy.arange(rev_mapping.shape[0])
            if not data is None:
                data[:, 0] = mapping[data[:, 0] - startfrag]
                data[:, 1] = mapping[data[:, 1] - startfrag]
            if not trans_data is None:
                trans_data[:, 0] = mapping[trans_data[:, 0]]
                trans_data[:, 1] = mapping[trans_data[:, 1]]
            mids = mids[rev_mapping]
            strands = strands[rev_mapping]
            if not self.silent:
                print >> sys.stderr, ("\r%s\rChecking for fragment interaction count...") % (' ' * 80),
            # precalculate interaction distance means for all included interactions
            if not data is None:
                counts = data[:, 2].astype(numpy.float64)
            else:
                counts = None
            if not trans_data is None:
                trans_counts = trans_data[:, 2].astype(numpy.float64)
            else:
                trans_counts = None
            trans_means = None
            distance_means = None
            if remove_distance:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rPrecalculating distances...") % (' ' * 80),
                if usereads != 'cis':
                    trans_mean = numpy.sum(trans_counts).astype(numpy.float64)
                    ffrags = numpy.where(strands == 0)[0]
                    rfrags = numpy.where(strands == 1)[0]
                    interactions = ffrags.shape[0] * rfrags.shape[0]
                    all_ints = self.frags['fragments']['region'][rev_mapping]
                    fints = all_ints[ffrags]
                    rints = all_ints[rfrags]
                    interactions -= numpy.sum(
                            numpy.bincount(fints, minlength=self.frags['regions'].shape[0]) *
                            numpy.bincount(rints, minlength=self.frags['regions'].shape[0]))
                    trans_mean /= interactions
                    trans_means = numpy.empty(trans_data.shape[0], dtype=numpy.float32).fill(trans_mean)
                    if not data is None:
                        distance_means = numpy.zeros(data.shape[0], dtype=numpy.float32)
                        findices = numpy.r_[0, numpy.bincount(fints)]
                        rindices = numpy.r_[0, numpy.bincount(rints)]
                        for i in range(1, findices.shape[0]):
                            findices[i] += findices[i - 1]
                            rindices[i] += rindices[i - 1]
                        for i in range(findices.shape[0] - 1):
                            if findices[i] < findices[i + 1] and rindices[i] < rindices[i + 1]:
                                distance_means[:] = numpy.exp(-self.gamma *
                                        numpy.log(mids[data[:, 1]] - mids[data[:, 0]]) +
                                        self.region_means[all_ints[data[:, 0]]])
                else:
                    distance_means = numpy.zeros(data.shape[0], dtype=numpy.float32)
                    distance_means[:] = (-self.gamma * numpy.log(mids[data[:, 1]] - mids[data[:, 0]]) +
                                         self.region_means[region])
            if precorrect:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rFinding binning corrections...") % (' ' * 80),
                if not data is None:
                    if distance_means is None:
                        distance_means = numpy.ones(data.shape[0], dtype=numpy.float32)
                    _optimize.find_binning_correction_adjustment(distance_means,
                                                                 data,
                                                                 self.binning_corrections,
                                                                 self.binning_correction_indices,
                                                                 self.binning_num_bins,
                                                                 self.binning_frag_indices)
                if not trans_data is None:
                    if trans_means is None:
                        trans_means = numpy.ones(trans_data.shape[0], dtype=numpy.float32)
                    _optimize.find_binning_correction_adjustment(trans_means,
                                                                 trans_data,
                                                                 self.binning_corrections,
                                                                 self.binning_correction_indices,
                                                                 self.binning_num_bins,
                                                                 self.binning_frag_indices)
            if not distance_means is None:
                counts /= numpy.exp(distance_means)
            if not trans_means is None:
                trans_counts /= numpy.exp(trans_means)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding fend corrections...") % (' ' * 80),
            # add psuedo-count diagonal
            if data is None:
                data = numpy.zeros((rev_mapping.shape[0], 2), dtype=numpy.int32)
                data[:, 0] = numpy.arange(rev_mapping.shape[0])
                data[:, 1] = numpy.arange(rev_mapping.shape[0])
                counts = numpy.ones(data.shape[0], dtype=numpy.float64) * 0.5
            else:
                temp = numpy.zeros((rev_mapping.shape[0], 3), dtype=numpy.int32)
                temp[:, 0] = numpy.arange(rev_mapping.shape[0])
                temp[:, 1] = numpy.arange(rev_mapping.shape[0])
                data = numpy.vstack((data, temp))
                counts = numpy.hstack((counts, numpy.ones(data.shape[0], dtype=numpy.float64) * 0.5))
            # calculate corrections
            corrections = numpy.ones((rev_mapping.shape[0], 1), dtype=numpy.float64)
            g = 0.9
            eta = etamax = 0.1
            tol = 1e-12
            stop_tol = tol * 0.5
            rt = tol ** 2.0
            delta = 0.1
            Delta = 3
            v = numpy.zeros((corrections.shape[0], 1), dtype=numpy.float64)
            w = numpy.zeros((corrections.shape[0], 1), dtype=numpy.float64)
            _optimize.calculate_v(data, trans_data, counts, trans_counts, corrections, v)
            rk = 1.0 - v
            rho_km1 = numpy.dot(rk.T, rk)[0, 0]
            rho_km2 = rho_km1
            rold = rout = rho_km1
            i = MVP = 0
            while rout > rt:
                i += 1
                k = 0
                y = numpy.ones((rev_mapping.shape[0], 1), dtype=numpy.float64)
                innertol = max(eta ** 2.0 * rout, rt)
                while rho_km1 > innertol:
                    k += 1
                    if k == 1:
                        Z = rk / v
                        p = numpy.copy(Z)
                        rho_km1 = numpy.dot(rk.T, Z)
                    else:
                        beta = rho_km1 / rho_km2
                        p = Z + beta * p
                    # Update search direction efficiently
                    w.fill(0.0)
                    _optimize.calculate_w(data, trans_data, counts, trans_counts, corrections, p, w)
                    w += v * p
                    alpha = rho_km1 / numpy.dot(p.T, w)[0, 0]
                    ap = alpha * p
                    # Test distance to boundary of cone
                    ynew = y + ap
                    if numpy.amin(ynew) <= delta:
                        if delta == 0:
                            break
                        ind = numpy.where(ap < 0.0)[0]
                        gamma = numpy.amin((delta - y[ind]) / ap[ind])
                        y += gamma * ap
                        break
                    if numpy.amax(ynew) >= Delta:
                        ind = numpy.where(ynew > Delta)[0]
                        gamma = numpy.amin((Delta - y[ind]) / ap[ind])
                        y += gamma * ap
                        break
                    y = numpy.copy(ynew)
                    rk -= alpha * w
                    rho_km2 = rho_km1
                    Z = rk / v
                    rho_km1 = numpy.dot(rk.T, Z)[0, 0]
                corrections *= y
                v.fill(0.0)
                _optimize.calculate_v(data, trans_data, counts, trans_counts, corrections, v)
                rk = 1.0 - v
                rho_km1 = numpy.dot(rk.T, rk)[0, 0]
                rout = rho_km1
                MVP += k + 1
                # Update inner iteration stopping criterion
                rat = rout / rold
                rold = rout
                res_norm = rout ** 0.5
                eta_o = eta
                eta = g * rat
                if g * eta_o ** 2.0 > 0.1:
                    eta = max(eta, g * eta_o ** 2.0)
                eta = max(min(eta, etamax), stop_tol / res_norm)
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rIteration %i Residual: %e") % (" " * 80, i, rout),
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding fragment corrections... Region: %s Done\n") % (' ' * 80, str(region)),
            self.corrections[rev_mapping + startfrag] = numpy.log(1.0 / corrections)
        # calculate chromosome mean
        if self.region_means is None:
            self.region_means = numpy.zeros(self.frags['regions'].shape[0], dtype=numpy.float32)
        for i in all_regions:
            start = self.frags['regions']['start_frag'][i]
            stop = self.frags['regions']['stop_frag'][i]
            forward = (numpy.where(self.filter[start:stop] * 
                                   (self.frags['fragments']['strand'][start:stop] == 0))[0] + start)
            reverse = (numpy.where(self.filter[start:stop] * 
                                   (self.frags['fragments']['strand'][start:stop] == 1))[0] + start)
            if forward.shape[0] == 0 or reverse.shape[0] == 0:
                continue
            region_mean = (numpy.sum(self.corrections[forward]) * reverse.shape[0] +
                           numpy.sum(self.corrections[reverse]) * forward.shape[0])
            region_mean /= forward.shape[0] * reverse.shape[0]
            self.corrections[forward] -= region_mean / 2.0
            self.corrections[reverse] -= region_mean / 2.0
            if remove_distance:
                self.region_means[i] += region_mean
        if not self.silent:
            print >> sys.stderr, ("\r%s\rCompleted learning express corrections.\n") % (' ' * 80),
        if precorrect:
            self.normalization = 'binning-express'
        else:
            self.normalization = 'express'
        self.history += "Succcess\n"
        return None

    def find_binning_fragment_corrections(self, mindistance=0, maxdistance=0, model=['gc', 'len'], num_bins=[10, 10],
                                          parameters=['even', 'even'], learning_threshold=1.0, max_iterations=100,
                                          usereads='cis', regions=[], precorrect=False):
        """
        Using multivariate binning model, learn correction values for combinations of model parameter bins.

        :param mindistance: The minimum inter-fend distance to be included in modeling.
        :type mindistance: int.
        :param maxdistance: The maximum inter-fend distance to be included in modeling.
        :type maxdistance: int.
        :param model: A list of fragment features to be used in model. Valid values are 'len' and any features included in the creation of the associated Fragment object.
        :type model: list
        :param num_bins: The number of approximately equal-sized bins two divide model components into.
        :type num_bins: int.
        :param parameters: A list of types, one for each model parameter. Types can be either 'even' or 'fixed', indicating whether each parameter bin should contain approximately even numbers of interactions or be of fixed width spanning 1 / Nth of the range of the parameter's values, respectively. Parameter types can also have the suffix '-const' to indicate that the parameter should not be optimized.
        :type parameters: list
        :param remove_distance: Use distance dependence curve in prior probability calculation for each observation.
        :type remove_distance: bool.
        :param learning_threshold: The minimum change in log-likelihood needed to continue iterative learning process.
        :type learning_threshold: float
        :param max_iterations: The maximum number of iterations to use for learning model parameters.
        :type max_iterations: int.
        :param usereads: Specifies which set of interactions to use, 'cis', 'trans', and 'all'.
        :type usereads: str.
        :param regions: A list of regions to calculate corrections for. If set as None, all region corrections are found.
        :type regions: list
        :param precorrect: Use fragment-based corrections in expected value calculations, resulting in a chained normalization approach.
        :type precorrect: bool.
        :returns: None

        :Attributes: * **model_parameters** (*ndarray*) - A numpy array of strings containing model parameter names.
                     * **binning_num_bins** (*ndarray*) - A numpy array of type int32 containing the number of bins for each model parameter.
                     * **binning corrections** (*ndarray*) - A numpy array of type float32 and length equal to the sum of binning_num_bins * (binning_num_bins - 1) / 2. This array contains a 1D stack of correction values, ordered according to the parameter order in the 'model_parameters' attribute.
                     * **binning_correction_indices** (*ndarray*) - A numpy array of type int32 and length equal to the number of model parameters plus one. This array contains the first position in 'binning_corrections' for the first bin of the model parameter in the corresponding position in the 'model_parameters' array. The last position in the array contains the total number of binning correction values.
                     * **binning_frag_indices** (*ndarray*) - A numpy array of type int32 and size N x M where M is the number of model parameters and N is the number of fragments. This array contains the binning index for each parameter for each fragment.
        
        The 'normalization' attribute is updated to 'binning', 'probability-binning', or 'express-binning', depending on if the 'precorrect' option is selected and which normalization has been previously run.
        """
        self.history += "FiveC.find_binning_fragment_corrections(mindistance=%s, maxdistance=%s, num_bins=%s, model=%s, learning_threshold=%f, max_iterations=%i, usereads='%s', regions=%s) - " % (str(mindistance), str(maxdistance), str(num_bins), str(model), learning_threshold, max_iterations, usereads, str(regions))
        for parameter in model:
            if not parameter in ['len'] and parameter not in self.frags['fragments'].dtype.names:
                if not self.silent:
                    print >> sys.stderr, ("Fragment feature %s not found in fragment object. Try removing it from model or creating a new fragment object with feature data.\n") % (parameter),
                self.history += "Error: model parameter '%s' not found in fragment data\n" % parameter
                return None
        for parameter in parameters:
            if parameter not in ['even', 'fixed', 'even-const', 'fixed-const']:
                if not self.silent:
                    print >> sys.stderr, ("Fragment feature type %s is not valid.") % (parameter),
                self.history += "Error: model feature type '%s' not valid\n" % parameter
                return None
        if len(model) != len(num_bins):
            if not self.silent:
                print >> sys.stderr, ("The number of items in the 'model' parameter must be the same as the number in the 'num_bins' parameter.\n"),
            self.history += "Error: mismatch between lengths of 'num_bins' and 'model'\n"
            return None
        # make sure usereads has a valid value
        if usereads not in ['cis', 'trans', 'all']:
            if not self.silent:
                print >> sys.stderr, ("'usereads' does not have a valid value.\n"),
            self.history += "Error: '%s' not a valid value for 'usereads'\n" % usereads
            return None
        if precorrect and self.corrections is None:
            if not self.silent:
                print >> sys.stderr, ("Precorrection can only be used in project has previously run 'find_probability_fragment_corrections' or 'find_express_fragment_corrections'.\n"),
            self.history += "Error: 'find_binning_fragment_corrections()' or 'find_binning_fragment_corrections()' not run yet\n"
            return None
        # if regions not given, set to all regions
        if regions == None or len(regions) == 0:
            regions = numpy.arange(self.frags['regions'].shape[0])
        # limit corrections to only requested regions
        filt = numpy.copy(self.filter)
        for i in range(self.frags['regions'].shape[0]):
            if i not in regions:
                filt[self.frags['regions']['start_frag'][i]:self.frags['regions']['stop_frag'][i]] = 0
        if maxdistance == 0 or maxdistance is None:
            for i in range(self.frags['regions'].shape[0]):
                maxdistance = max(maxdistance, self.frags['regions']['stop'][i] -
                                               self.frags['regions']['start'][i]) + 1
        if not self.silent:
            print >> sys.stderr, ("\r%s\rPartitioning features into bins...") % (' ' * 80),
        num_bins = numpy.array(num_bins, dtype=numpy.int32)
        total_bins = 1
        all_bins = numpy.zeros(0, dtype=numpy.float32)
        all_corrections = numpy.ones(0, dtype=numpy.float64)
        all_indices = numpy.zeros((filt.shape[0], len(model)), dtype=numpy.int32)
        bin_indices = numpy.zeros(len(model) + 1, dtype=numpy.int32)
        correction_indices = numpy.zeros(len(model) + 1, dtype=numpy.int32)
        bin_divs = numpy.zeros(len(model), dtype=numpy.int32)
        for i in range(len(model)):
            if model[i] == 'len':
                values = (self.frags['fragments']['stop'][...] -
                          self.frags['fragments']['start'][...]).astype(numpy.float32)
            else:
                values = self.frags['fragments'][model[i]][...]
            if parameters[i].count('even') > 0:
                temp = numpy.copy(values)
                temp.sort()
                all_bins = numpy.hstack((all_bins, temp[numpy.round(numpy.linspace(0, values.shape[0],
                                        num_bins[i] + 1)).astype(numpy.int32)[1:] - 1])).astype(numpy.float32)
            else:
                all_bins = numpy.hstack((all_bins, numpy.linspace(numpy.amin(values),
                                        numpy.amax(values), num_bins[i] + 1)[1:])).astype(numpy.float32)
            all_bins[-1] = numpy.inf
            bin_indices[i + 1] = all_bins.shape[0]
            all_corrections = numpy.hstack((all_corrections, numpy.zeros(num_bins[i] * (num_bins[i] + 1) / 2,
                                           dtype=numpy.float64)))
            correction_indices[i + 1] = all_corrections.shape[0]
            bin_divs[i] = total_bins
            total_bins *= num_bins[i] * (num_bins[i] + 1) / 2
            all_indices[:, i] = numpy.searchsorted(all_bins[bin_indices[i]:bin_indices[i + 1]],
                                                   values).astype(numpy.int32)
        self.binning_frag_indices = all_indices
        self.binning_num_bins = num_bins
        mids = self.frags['fragments']['mid'][...]
        region_ints = self.frags['fragments']['region'][...]

        bin_counts = numpy.zeros(total_bins, dtype=numpy.int64)
        bin_sums = numpy.zeros((total_bins, 2), dtype=numpy.float64)
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinding bin counts...") % (' ' * 80),
        # Find number of observations in each bin
        data = None
        trans_data = None
        trans_mean = 0.0
        gamma = 0.0
        trans_mean = 0.0
        frag_corrections = None
        if precorrect:
            frag_corrections = self.corrections
        if usereads in ['cis', 'all']:
            data = self.data['cis_data'][...]
            data = data[numpy.where(filt[data[:, 0]] * filt[data[:, 1]])[0], :]
            if self.gamma is None:
                self.find_distance_parameters()
            gamma = self.gamma
        if usereads in ['trans', 'all']:
            trans_data = self.data['trans_data'][...]
            trans_data = trans_data[numpy.where(filt[trans_data[:, 0]] * filt[trans_data[:, 1]])[0], :] 
            if self.trans_mean is None:
                self.find_trans_mean()
            trans_mean = self.trans_mean
        _binning.binning_bin_observed(data,
                                      trans_data,
                                      region_ints,
                                      mids,
                                      bin_counts,
                                      bin_sums,
                                      all_indices,
                                      num_bins,
                                      bin_divs,
                                      self.region_means,
                                      frag_corrections,
                                      mindistance,
                                      maxdistance,
                                      gamma,
                                      trans_mean)
        # Find seed values
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinding seed values...") % (' ' * 80),
        all_indices = numpy.zeros((bin_counts.shape[0], len(model)), dtype=numpy.int32)
        n = numpy.sum(bin_counts)
        for i in range(correction_indices.shape[0] - 1):
            all_indices[:, i] = ((numpy.arange(bin_counts.shape[0], dtype=numpy.int32) / bin_divs[i]) %
                                 (correction_indices[i + 1] - correction_indices[i]))
            temp0 = numpy.bincount(all_indices[:, i], weights=bin_sums[:, 0], minlength=num_bins[i])
            temp1 = numpy.bincount(all_indices[:, i], weights=bin_counts, minlength=num_bins[i])
            where = numpy.where(temp1 > 0)[0]
            all_corrections[where + correction_indices[i]] = temp0[where] / temp1[where].astype(numpy.float64)

        def find_cor_sums(index, indices, corrections, correction_indices, cor_sums):
            cor_sums.fill(0.0)
            for i in range(indices.shape[1]):
                if i == index:
                    continue
                cor_sums += corrections[correction_indices[i] + indices[:, i]]
            return

        def find_sigma2(counts, sums, cor_sums, n):
            return 1.0 / (n - 1.0) * numpy.sum(counts * cor_sums ** 2.0 + sums[:, 1] - 2 * cor_sums * sums[:, 0])

        def find_ll(counts, sums, cor_sums, n, sigma, sigma2):
            return (n * (log(sigma) + 0.5 * log(2.0 * numpy.pi)) + 1.0 / (2.0 * sigma2) *
                    numpy.sum(counts * cor_sums ** 2.0 + sums[:, 1] - 2.0 * cor_sums * sums[:, 0]))
        
        def temp_ll(x, *args):
            counts, sums, cor_sums, n, ni, sigma, sigma2 = args[:7]
            return (ni * (log(sigma) + 0.5 * log(2.0 * numpy.pi)) + 1.0 / (2.0 * sigma2) *
                    numpy.sum(counts * (cor_sums + x[0]) ** 2.0 + sums[:, 1] - 2.0 * (cor_sums + x[0]) * sums[:, 0]))

        def temp_ll_grad(x, *args):
            counts, sums, cor_sums, n, ni, sigma, sigma2 = args[:7]
            cor_sums1 = cor_sums + x[0]
            cor_sums2 = cor_sums1 ** 2.0
            M = numpy.sum(counts * cor_sums2 + sums[:, 1] - 2.0 * cor_sums1 * sums[:, 0])
            N = 1.0 / (n - 1.0) * numpy.sum(2.0 * counts * cor_sums1 - 2.0 * sums[:, 0])
            O = N / (2.0 * sigma)
            grad = (ni / sigma * O - M * N / (2.0 * sigma2 ** 2.0) + 1.0 / sigma2 *
                    (ni * x[0] + numpy.sum(counts * cor_sums) - numpy.sum(sums[:, 0])))
            return numpy.array(grad, dtype=numpy.float64)

        cor_sums = numpy.zeros(bin_counts.shape[0], dtype=numpy.float64)
        find_cor_sums(-1, all_indices, all_corrections, correction_indices, cor_sums)
        sigma2 = find_sigma2(bin_counts, bin_sums, cor_sums, n)
        sigma = sigma2 ** 0.5
        ll = find_ll(bin_counts, bin_sums, cor_sums, n, sigma, sigma2)
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning binning corrections... iteration:00  ll:%f") % (' ' * 80, ll),
        iteration = 0
        delta = numpy.inf
        pgtol = 1e-8
        old_settings = numpy.seterr(invalid='ignore', divide='ignore')
        while iteration < max_iterations and delta >= learning_threshold:
            new_corrections = numpy.copy(all_corrections)
            for h in range(len(model)):
                # don't learn constant parameters
                if parameters[h].count('const') > 0:
                    continue
                num_cor = correction_indices[h + 1] - correction_indices[h]
                bins = bin_counts.shape[0] / num_cor
                temp_cor_sums = numpy.zeros(bins, dtype=numpy.float64)
                for i in range(num_cor):
                    where = numpy.where(all_indices[:, h] == i)[0]
                    temp_bin_counts = bin_counts[where]
                    temp_bin_sums = bin_sums[where]
                    find_cor_sums(h, all_indices[where, :], all_corrections, correction_indices, temp_cor_sums)
                    ni = numpy.sum(temp_bin_counts)
                    x0 = all_corrections[(correction_indices[h] + i):(correction_indices[h] + i + 1)]
                    x, f, d = bfgs(func=temp_ll, x0=x0, fprime=temp_ll_grad, pgtol=pgtol,
                                   args=(temp_bin_counts, temp_bin_sums, temp_cor_sums, n, ni, sigma, sigma2))
                    new_corrections[correction_indices[h] + i] = x[0]
            all_corrections = new_corrections
            iteration += 1
            find_cor_sums(-1, all_indices, all_corrections, correction_indices, cor_sums)
            sigma2 = find_sigma2(bin_counts, bin_sums, cor_sums, n)
            sigma = sigma2 ** 0.5
            new_ll = find_ll(bin_counts, bin_sums, cor_sums, n, sigma, sigma2)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rLearning binning corrections... iteration:%02i  ll:%f\n") % (' ' * 80,
                                      iteration, new_ll),
            delta = ll - new_ll
            if delta < 0.0:
                delta = numpy.inf
            ll = new_ll
        numpy.seterr(**old_settings)
        self.binning_corrections = all_corrections.astype(numpy.float32)
        self.model_parameters = numpy.array(model)
        self.binning_correction_indices = correction_indices
        if not self.silent:
            print >> sys.stderr, ("\r%s\rLearning binning corrections... Final ll:%f\n") % (' ' * 80, ll),
        if precorrect:
            if self.normalization == 'probability':
                self.normalization = 'probability-binning'
            else:
                self.normalization = 'express-binning'
        else:
            self.normalization = 'binning'
        self.history += "Success\n"
        return None

    def find_trans_mean(self):
        """
        Calculate the mean signal across all valid fragment-pair trans (inter-region) interactions.
        
        :returns: None

        :Attributes: * **trans_mean** (*float*) - A float corresponding to the mean signal of inter-region interactions.
        """
        self.history += "FiveC.find_trans_mean() - "
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
        self.history += "Success\n"
        return None

    def cis_heatmap(self, region, binsize=0, binbounds=None, start=None, stop=None, startfrag=None, stopfrag=None,
                    datatype='enrichment', arraytype='full', skipfiltered=False, returnmapping=False,
                    dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                    removefailed=False, image_file=None, **kwargs):
        """
        Return a heatmap of cis data of the type and shape specified by the passed arguments.

        This function returns a heatmap for a single region, bounded by either 'start' and 'stop' or 'startfrag' and 'stopfrag' ('start' and 'stop' take precedence). If neither is given, the complete region is included. The data in the array is determined by the 'datatype', being raw, fragment-corrected, distance-corrected, enrichment, or expected data. The array shape is given by 'arraytype' and can be compact (if unbinned), upper, or full. See :mod:`fivec_binning <hifive.fivec_binning>` for further explanation of 'datatype' and 'arraytype'. If using dynamic binning ('dynamically_binned' is set to True), 'minobservations', 'searchdistance', 'expansion_binsize', and 'removefailed' are used to control the dynamic binning process. Otherwise these arguments are ignored.

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
        if datatype not in ['raw', 'fragment', 'distance', 'enrichment', 'expected']:
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
            data = fivec_binning.find_cis_signal(self, region, binsize=binsize, binbounds=binbounds,
                                                 start=start, stop=stop, startfrag=startfrag,
                                                 stopfrag=stopfrag, datatype=datatype, arraytype=arraytype,
                                                 skipfiltered=skipfiltered, returnmapping=returnmapping,
                                                 silent=self.silent)
        else:
            expansion, exp_mapping = fivec_binning.find_cis_signal(self, region, binsize=expansion_binsize,
                                                                   start=start, stop=stop, startfrag=startfrag,
                                                                   stopfrag=stopfrag, datatype=datatype,
                                                                   arraytype='full', skipfiltered=True,
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
        if datatype not in ['raw', 'fragment', 'distance', 'enrichment', 'expected']:
            if not self.silent:
                print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
            return None
        if arraytype not in ['compact', 'full'] or (dynamically_binned and arraytype == 'compact'):
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
                                                                         datatype=datatype, arraytype=arraytype,
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

    def write_heatmap(self, filename, binsize, includetrans=True, datatype='enrichment', arraytype='full',
                      regions=[], dynamically_binned=False, minobservations=0,  searchdistance=0, expansion_binsize=0,
                      removefailed=False):
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
        :param dynamically_binned: If 'True', return dynamically binned data.
        :type dynamically_binned: bool.
        :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
        :type minobservations: int.
        :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
        :type searchdistance: int.
        :param expansion_binsize: The size of bins to use for data to pull from when expanding dynamic bins. If set to zero, unbinned data is used.
        :type expansion_binsize: int.
        :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
        :returns: None

        The following attributes are created within the hdf5 dictionary file. Arrays are accessible as datasets while the resolution is held as an attribute.

        :Attributes: * **resolution** (*int.*) - The bin size that data are accumulated in.
                     * **regions** (*ndarray*) - A numpy array containing region data for each region included in the heatmaps.
                     * **N.positions** (*ndarray*) - A series of numpy arrays of type int32, one for each region where N is the region index, containing one row for each bin and four columns denoting the start and stop coordinates and first fragment and last fragment plus one for each bin. This is included if data is in the 'full' format.
                     * **N.forward_positions** (*ndarray*) - A series of numpy arrays of type int32, one for each region where N is the region index, containing one row for each bin along the first axis and four columns denoting the start and stop coordinates and first fragment and last fragment plus one for each bin. This is included if data is in the 'compact' format and corresponds to only forward strand fragments.
                     * **N.reverse_positions** (*ndarray*) - A series of numpy arrays of type int32, one for each region where N is the region index, containing one row for each bin along the second axis and four columns denoting the start and stop coordinates and first fragment and last fragment plus one for each bin. This is included if data is in the 'compact' format and corresponds to only reverse strand fragments.
                     * **N.counts** (*ndarray*) - A series of numpy arrays of type int32, one for each region where N is the region index, containing the observed counts for valid fragment combinations. If arrays are in the 'compact' format, the first axis corresponds to forward fragments and the second axis corresponds to reverse fragments. If the array is in the 'upper' format,data are in an upper-triangle format such that they have N * (N - 1) / 2 entries where N is the number of fragments or bins in the region.
                     * **N.expected** (*ndarray*) - A series of numpy arrays of type float32, one for each region where N is the region index, containing the expected counts for valid fragment combinations. If the array is in the 'upper' format,data are in an upper-triangle format such that they have N * (N - 1) / 2 entries where N is the number of fragments or bins in the region.
                     * **N_by_M.counts** (*ndarray*) - A series of numpy arrays of type int32, one for each region pair N and M if trans data are included, containing the observed counts for valid fragment combinations. The region index order specifies which axis corresponds to which region. If data are in the 'compact' format, both region index orders will be present.
                     * **N_by_M.expected** (*ndarray*) - A series of numpy arrays of type float32, one for each region pair N and M if trans data are included, containing the expected counts for valid fend combinations. The chromosome name order specifies which axis corresponds to which region. If data are in the 'compact' format, both region index orders will be present.
        """
        history = self.history
        history += "FiveC.write_heatmap(filename='%s', binsize=%i, includetrans=%s, datatype='%s', arraytype='%s', regions=%s, dynamically_binned=%s, minobservations=%i, searchdistance=%i, expansion_binsize=%i, removefailed=%s)" % (filename, binsize, includetrans, datatype, arraytype, str(regions), dynamically_binned, minobservations, searchdistance, expansion_binsize, removefailed)
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
                                         regions=regions, silent=self.silent, history=history,
                                         dynamically_binned=dynamically_binned, minobservations=minobservations,
                                         searchdistance=searchdistance, expansion_binsize=expansion_binsize,
                                         removefailed=removefailed)
        return None
