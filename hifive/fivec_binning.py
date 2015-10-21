#!/usr/bin/env python

"""
This is a module contains scripts for generating compact, upper-triangle and full matrices of 5C interaction data.

Concepts
--------

Data can either be arranged in compact, complete, or flattened (row-major) upper-triangle arrays. Compact arrays are N x M, where N is the number of forward probe fragments and M is the number of reverse probe fragments. Data can be raw, fragment-corrected, distance-dependence removed, or enrichment values. Arrays are 3-dimensional with observed values in the first layer of d3, expected values in the second layer of d3. The exception to this is upper-triangle arrays, which are 2d, dividing observed and expected along the second axis.

API documentation
-----------------
"""

import os
import sys
import subprocess

import numpy
import h5py

import libraries._fivec_binning as _fivec_binning


def find_cis_signal(fivec, region, binsize=0, binbounds=None, start=None, stop=None, startfrag=None, stopfrag=None,
                    datatype='enrichment', arraytype='full', skipfiltered=False, returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    :param fivec:  A :class:`FiveC <hifive.fivec.FiveC>` class object containing fragment and count data.
    :type fivec: :class:`FiveC <hifive.fivec.FiveC>`
    :param region: The index of the region to pull data from.
    :type region: int.
    :param binsize: This is the coordinate width of each bin. A value of zero indicates unbinned. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fragment not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The smallest coordinate to include in the array, measured from fragment midpoints. If 'binbounds' is given, this value is ignored. If both 'start' and 'startfrag' are given, 'start' will override 'startfrag'. If unspecified, this will be set to the midpoint of the first fragment for 'region', adjusted to the first multiple of 'binsize' if not zero. Optional.
    :type start: int.
    :param stop: The largest coordinate to include in the array, measured from fragment midpoints. If 'binbounds' is given, this value is ignored. If both 'stop' and 'stopfrag' are given, 'stop' will override 'stopfrag'. If unspecified, this will be set to the midpoint of the last fragment plus one for 'region', adjusted to the last multiple of 'start' + 'binsize' if not zero. Optional.
    :type stop: int.
    :param startfrag: The first fragment to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'start' is not given, this is set to the first fragment in 'region'. In cases where 'start' is specified and conflicts with 'startfrag', 'start' is given preference. Optional.
    :type startfrag: int.
    :param stopfrag: The first fragment not to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'stop' is not given, this is set to the last fragment in 'region' plus one. In cases where 'stop' is specified and conflicts with 'stopfrag', 'stop' is given preference. Optional.
    :type stopfrag: int.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment', 'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fragments return value of one. Expected values are returned for 'distance', 'fragment', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both correction and distance mean values. 'enrichment' also scales both observed and expected by the standard deviation, giving a completely normalized set of values.
    :type datatype: str.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact' (though only when 'binned' is zero), 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N and M are the number of forward and reverse probe fragments, respectively. 'full' returns a square, symmetric array of size N x N x 2 where N is the total number of fragments. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2, where N is the total number of fragments.
    :type arraytype: str.
    :param skipfiltered: If 'True', all interaction bins for filtered out fragments are removed and a reduced-size array is returned.
    :type skipfiltered: bool.
    :param returnmapping: If 'True', a list containing the data array and either one or two 2d arrays containing first coordinate included and excluded from each bin, and the first fragment included and excluded from each bin  corresponding to both axes or the first and second axis for an upper or compact array, respectively, is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that all values are acceptable
    if datatype not in ['raw', 'fragment', 'distance', 'enrichment', 'expected']:
        if not silent:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    elif datatype in ['fragment', 'enrichment'] and fivec.normalization == 'none':
        if not silent:
            print >> sys.stderr, ("Normalization has not been performed yet on this project. Select either 'raw' or 'distance' for datatype. No data returned\n"),
        return None
    elif datatype in ['distance', 'enrichment'] and fivec.gamma is None:
        fivec.find_distance_parameters()
    if arraytype not in ['full', 'compact', 'upper']:
        if not silent:
            print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    if arraytype == 'compact' and (binsize > 0 or not binbounds is None):
        if not silent:
            print >> sys.stderr, ("'Compact' array can only be used with unbinned data. No data returned.\n"),
        return None
    # Determine start, stop, startfrag, and stopfrag
    chrom = fivec.frags['regions']['chromosome'][region]
    chrint = fivec.chr2int[chrom]
    if not binbounds is None:
        start = binbounds[0, 0]
        stop = binbounds[-1, 1]
        startfrag = _find_frag_from_coord(fivec, chrint, start)
        stopfrag = _find_frag_from_coord(fivec, chrint, stop)
    else:
        if start is None and startfrag is None:
            startfrag = fivec.frags['regions']['start_frag'][region]
            start = fivec.frags['fragments']['mid'][startfrag]
            if binsize > 0:
                start = (start / binsize) * binsize
        elif start is None:
            start = fivec.frags['fragments']['mid'][startfrag]
            if binsize > 0:
                start = (start / binsize) * binsize
        else:
            startfrag = _find_frag_from_coord(fivec, chrint, start)
        if (stop is None or stop == 0) and stopfrag is None:
            stopfrag = fivec.frags['regions']['stop_frag'][region]
            stop = fivec.frags['fragments']['mid'][stopfrag - 1] + 1
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        elif stop is None or stop == 0:
            stop = fivec.frags['fragments']['mid'][stopfrag - 1] + 1
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        else:
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
            stopfrag = _find_frag_from_coord(fivec, chrint, stop)
    if stopfrag - startfrag == 0:
        if not silent:
            print >> sys.stderr, ("Insufficient data, no data returned.\n"),
        return None
    if not silent:
        print >> sys.stderr, ("Finding %s %s array for %s:%i-%i...") % (datatype, arraytype, chrom, start, stop),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        start_index = fivec.data['cis_indices'][startfrag]
        stop_index = fivec.data['cis_indices'][stopfrag]
        if stop_index - start_index == 0:
            if not silent:
                print >> sys.stderr, ("Insufficient data, no data returned.\n"),
            return None
        data_indices = fivec.data['cis_indices'][startfrag:(stopfrag + 1)]
        data_indices -= data_indices[0]
        data = fivec.data['cis_data'][start_index:stop_index, :]
        data[:, :2] -= startfrag
    else:
        data_indices = None
        data = None
    # Determine mapping of valid fends to bins
    mids = fivec.frags['fragments']['mid'][startfrag:stopfrag]
    strands = fivec.frags['fragments']['strand'][startfrag:stopfrag]
    mapping = numpy.zeros(stopfrag - startfrag, dtype=numpy.int32) - 1
    valid = numpy.where(fivec.filter[startfrag:stopfrag] > 0)[0]
    if valid.shape[0] == 0:
        if not silent:
            print >> sys.stderr, ("Insufficient data\n"),
        return None
    if binsize == 0 and binbounds is None:
        if arraytype == 'compact':
            for_valid = numpy.where((strands == 0) * (fivec.filter[startfrag:stopfrag] > 0))[0]
            rev_valid = numpy.where((strands == 1) * (fivec.filter[startfrag:stopfrag] > 0))[0]
            mapping.fill(0)
            if skipfiltered:
                mapping[for_valid] = numpy.arange(for_valid.shape[0]) + 1
                mapping[rev_valid] = -1 - numpy.arange(rev_valid.shape[0])
                num_for_bins = for_valid.shape[0]
                num_rev_bins = rev_valid.shape[0]
            else:
                num_for_bins = numpy.sum(strands == 0)
                num_rev_bins = numpy.sum(strands == 1)
                mapping[for_valid] = numpy.arange(num_for_bins) + 1
                mapping[rev_valid] = -1 - numpy.arange(num_rev_bins)
        else:
            if skipfiltered:
                mapping[valid] = numpy.arange(valid.shape[0])
                num_bins = valid.shape[0]
            else:
                mapping[valid] = valid
                num_bins = mapping.shape[0]
    elif not binbounds is None:
        start_indices = numpy.searchsorted(binbounds[:, 0], mids[valid], side='right') - 1
        stop_indices = numpy.searchsorted(binbounds[:, 1], mids[valid], side='right')
        where = numpy.where(start_indices == stop_indices)[0]
        mapping[valid[where]] = start_indices[where]
        num_bins = binbounds.shape[0]
    else:
        mapping[valid] = (mids[valid] - start) / binsize
        num_bins = (stop - start) / binsize
    if arraytype != 'compact' and num_bins < 2:
        if not silent:
            print >> sys.stderr, ("Insufficient data\n"),
        return None
    # If correction is required, determine what type and get appropriate data
    corrections = None
    binning_corrections = None
    correction_indices = None
    binning_num_bins = None
    frag_indices = None
    if datatype in ['fragment', 'enrichment', 'expected']:
        if fivec.normalization.count('probability') + fivec.normalization.count('express') > 0:
            corrections = fivec.corrections[startfrag:stopfrag]
        if fivec.normalization.count('binning') > 0:
            binning_corrections = fivec.binning_corrections
            correction_indices = fivec.binning_correction_indices
            binning_num_bins = fivec.binning_num_bins
            frag_indices = fivec.binning_frag_indices
    if datatype in ['distance', 'enrichment', 'expected']:
        gamma = fivec.gamma
        region_mean = fivec.region_means[region]
    else:
        gamma = 0.0
        region_mean = 0.0
    # Create requested array and corresponding mapping
    if binsize == 0 and arraytype == 'compact':
        data_array = numpy.zeros((num_for_bins, num_rev_bins, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in data values
    if arraytype == 'compact':
        if datatype != 'expected':
            _fivec_binning.find_cis_compact_observed(data, data_indices, mapping, data_array)
        if datatype != 'raw':
            _fivec_binning.find_cis_compact_expected(mapping, corrections, binning_corrections, correction_indices,
                                                     binning_num_bins, frag_indices, mids, strands, data_array, gamma,
                                                     region_mean, startfrag)
        else:
            where = numpy.where(data_array[:, :, 0] > 0.0)
            data_array[where[0], where[1], 1] = 1.0
        if datatype == 'expected':
            where = numpy.where(data_array[:, :, 1] > 0.0)
            data_array[where[0], where[1], 0] = 1.0
    else:
        if datatype != 'expected':
            _fivec_binning.find_cis_upper_observed(data, data_indices, mapping, data_array)
        if datatype != 'raw':
            _fivec_binning.find_cis_upper_expected(mapping, corrections, binning_corrections, correction_indices,
                                                   binning_num_bins, frag_indices, mids, strands, data_array, gamma,
                                                   region_mean, startfrag)
        else:
            where = numpy.where(data_array[:, 0] > 0.0)
            data_array[where[0], 1] = 1.0
        if datatype == 'expected':
            where = numpy.where(data_array[:, 1] > 0.0)
            data_array[where[0], 0] = 1.0
    # If requesting 'full' array, convert 'upper' array type to 'full'
    if arraytype == 'full':
        indices = numpy.triu_indices(num_bins, 1)
        full_data_array = numpy.zeros((num_bins, num_bins, 2), dtype=numpy.float32)
        full_data_array[indices[1], indices[0], :] = data_array
        full_data_array[indices[0], indices[1], :] = data_array
        del data_array
        data_array = full_data_array
    if returnmapping:
        if arraytype == 'compact':
            bin_mapping1 = numpy.zeros((data_array.shape[0], 4), dtype=numpy.int32)
            bin_mapping2 = numpy.zeros((data_array.shape[1], 4), dtype=numpy.int32)
            if skipfiltered:
                bin_mapping1[:, 2] = for_valid + startfrag
                bin_mapping2[:, 2] = rev_valid + startfrag
            else:
                bin_mapping1[:, 2] = numpy.where(strands == 0)[0] + startfrag
                bin_mapping2[:, 2] = numpy.where(strands == 1)[0] + startfrag
            bin_mapping1[:, 3] = bin_mapping1[:, 2] + 1
            bin_mapping1[:, 0] = fivec.frags['fragments']['start'][bin_mapping1[:, 2]]
            bin_mapping1[:, 1] = fivec.frags['fragments']['stop'][bin_mapping1[:, 2]]
            bin_mapping2[:, 3] = bin_mapping2[:, 2] + 1
            bin_mapping2[:, 0] = fivec.frags['fragments']['start'][bin_mapping2[:, 2]]
            bin_mapping2[:, 1] = fivec.frags['fragments']['stop'][bin_mapping2[:, 2]]
            if not silent:
                print >> sys.stderr, ("Done\n"),
            return [data_array, bin_mapping1, bin_mapping2]
        else:
            bin_mapping = numpy.zeros((num_bins, 4), dtype=numpy.int32)
            if binsize == 0 and binbounds is None:
                if skipfiltered:
                    bin_mapping[:, 2] = valid + startfrag
                else:
                    bin_mapping[:, 2] = numpy.arange(startfrag, stopfrag)
                bin_mapping[:, 3] = bin_mapping[:, 2] + 1
                bin_mapping[:, 0] = fivec.frags['fragments']['start'][bin_mapping[:, 2]]
                bin_mapping[:, 1] = fivec.frags['fragments']['stop'][bin_mapping[:, 2]]
            else:
                if binbounds is None:
                    bin_mapping[:, 0] = start + binsize * numpy.arange(num_bins)
                    bin_mapping[:, 1] = bin_mapping[:, 0] + binsize
                else:
                    bin_mapping[:, :2] = binbounds
                bin_mapping[:, 2] = numpy.searchsorted(mids, bin_mapping[:, 0]) + startfrag
                bin_mapping[:, 3] = numpy.searchsorted(mids, bin_mapping[:, 1]) + startfrag
            if not silent:
                print >> sys.stderr, ("Done\n"),
            return [data_array, bin_mapping]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return data_array


def _find_frag_from_coord(fivec, chrint, coord):
    """Find the next fragment after the coordinate on chromosome 'chrint'."""
    first_frag = fivec.frags['chr_indices'][chrint]
    last_frag = fivec.frags['chr_indices'][chrint + 1]
    return numpy.searchsorted(fivec.frags['fragments']['mid'][first_frag:last_frag], coord) + first_frag


def bin_cis_array(data_array, data_mapping, binsize=10000, binbounds=None, start=None, stop=None, arraytype='full',
                  returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill 'binsize' bins or bins defined by 'binbounds' with data provided in the array passed by 'data_array'.

    :param data_array: A 2d (upper) or 3d (full) array containing data to be binned. Array format will be determined from the number of dimensions.
    :type data_array: numpy array
    :param data_mapping: An N x 4 2d integer array containing the start and stop coordinates, and start and stop fragments for each of the N bin ranges in 'data_array'.
    :type data_mapping: numpy array
    :param binsize: This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any bin from 'data_array' not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The coordinate at the beginning of the first bin of the binned data. If unspecified, 'start' will be the first multiple of 'binsize' below the first coordinate from 'data_mapping'. If 'binbounds' is given, 'start' is ignored. Optional.
    :type start: int.
    :param stop: The coordinate at the end of the last bin of the binned data. If unspecified, 'stop' will be the first multiple of 'binsize' after the last coordinate from 'data_mapping'. If needed, 'stop' is adjusted upward to create a complete last bin. If 'binbounds' is given, 'stop' is ignored. Optional.
    :type stop: int.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'full' and 'upper'. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    :type arraytype: str.
    :param returnmapping: If 'True', a list containing the data array and a 2d array containing first coordinate included and excluded from each bin, and the first fragment included and excluded from each bin is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing binned data requested with 'datatype' pulled from 'data_array' or list of binned data array and mapping array.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that arraytype value is acceptable
    if arraytype not in ['full', 'upper']:
        if not silent:
            print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine input array type
    if len(data_array.shape) == 2 and data_mapping.shape[0] * (data_mapping.shape[0] - 1) / 2 == data_array.shape[0]:
        input_type = 'upper'
    elif len(data_array.shape) == 3 and data_array.shape[0] == data_mapping.shape[0]:
        input_type = 'full'
    else:
        if not silent:
            print >> sys.stderr, ("Unrecognized input array type. No data returned.\n"),
        return None
    # Determine start and stop, if necessary
    if binbounds is None:
        if start is None:
            start = (data_mapping[0, 0] / binsize) * binsize
        if stop is None:
            stop = ((data_mapping[-1, 1] - 1) / binsize + 1) * binsize
        else:
            stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        num_bins = (stop - start) / binsize
        binbounds = numpy.zeros((num_bins, 2), dtype=numpy.int32)
        binbounds[:, 0] = numpy.arange(num_bins) * binsize + start
        binbounds[:, 1] = binbounds[:, 0] + binsize
    else:
        num_bins = binbounds.shape[0]
        start = binbounds[0, 0]
        stop = binbounds[0, 1]
    mids = (data_mapping[:, 0] + data_mapping[:, 1]) / 2
    if not silent:
        print >> sys.stderr, ("Finding binned %s array...") % (arraytype),
    # Find bin mapping for each fend
    mapping = numpy.zeros(mids.shape[0], dtype=numpy.int32) - 1
    frag_ranges = numpy.zeros((binbounds.shape[0], 2), dtype=numpy.int32)
    for i in range(binbounds.shape[0]):
        firstbin = numpy.searchsorted(mids, binbounds[i, 0])
        lastbin = numpy.searchsorted(mids, binbounds[i, 1])
        mapping[firstbin:lastbin] = i
        frag_ranges[i, 0] = data_mapping[firstbin, 2]
        frag_ranges[i, 1] = data_mapping[lastbin, 3]
    # Create requested array
    binned_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in binned data values
    if input_type == 'upper':
        indices = numpy.triu_indices(data_array.shape[0], 1)
        _fivec_binning.bin_upper_to_upper(binned_array, data_array[indices[0], indices[1], :], mapping, num_bins)
    else:
        _fivec_binning.bin_upper_to_upper(binned_array, data_array, mapping, num_bins)
    # If requesting 'full' array, convert 'upper' array type to 'full'
    if arraytype == 'full':
        indices = numpy.triu_indices(num_bins, 1)
        full_binned_array = numpy.zeros((num_bins, num_bins, 2), dtype=numpy.float32)
        full_binned_array[indices[1], indices[0], :] = binned_array
        full_binned_array[indices[0], indices[1], :] = binned_array
        del binned_array
        binned_array = full_binned_array
    # If mapping requested, calculate bin bounds
    if returnmapping:
        mapping = numpy.zeros((num_bins, 4), dtype=numpy.int32)
        mapping[:, 0] = binbounds[:, 0]
        mapping[:, 1] = binbounds[:, 1]
        mapping[:, 2:4] = frag_ranges
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [binned_array, mapping]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return binned_array


def dynamically_bin_cis_array(unbinned, unbinnedpositions, binned, binbounds, minobservations=50,
                              searchdistance=0, removefailed=True, **kwargs):
    """
    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations', or 'searchdistance' criteria.

    :param unbinned: A full or upper array containing data to be binned. Array format will be determined from the number of dimensions.
    :type unbinned: numpy array
    :param unbinnedpositions: A 2d integer array indicating the first and last coordinate of each bin in 'unbinned' array.
    :type unbinnedpositions: numpy array
    :param binned: A full or upper array containing binned data to be dynamically binned. Array format will be determined from the number of dimensions. Data in this array will be altered by this function.
    :type binned: numpy array
    :param binbounds: A N x 2 integer array indicating the start and end position of each of N bins in 'binned' array.
    :type binbounds: numpy array
    :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
    :type minobservations: int.
    :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
    :type searchdistance: int.
    :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
    :type removefailed: bool.
    :returns: None
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # Determine unbinned array type
    if len(unbinned.shape) == 2 and (unbinnedpositions.shape[0] * (unbinnedpositions.shape[0] - 1) / 2 ==
                                     unbinned.shape[0]):
        ub_signal = unbinned
    elif len(unbinned.shape) == 3 and unbinned.shape[0] == unbinnedpositions.shape[0]:
        ub_signal = numpy.zeros((unbinned.shape[0] * (unbinned.shape[0] - 1) / 2, 2), dtype=numpy.float32)
        indices = numpy.triu_indices(unbinned.shape[0], 1)
        ub_signal[:, 0] = unbinned[indices[0], indices[1], 0]
        ub_signal[:, 1] = unbinned[indices[0], indices[1], 1]
    else:
        if not silent:
            print >> sys.stderr, ("Unrecognized unbinned array type. No data returned.\n"),
        return None
    # Determine binned array type
    if len(binned.shape) == 2 and binbounds.shape[0] * (binbounds.shape[0] - 1) / 2 == binned.shape[0]:
        binned_type = 'upper'
        b_signal = binned
    elif len(binned.shape) == 3 and binned.shape[0] == binbounds.shape[0]:
        binned_type = 'full'
        b_signal = numpy.zeros((binned.shape[0] * (binned.shape[0] - 1) / 2, 2), dtype=numpy.float32)
        indices = numpy.triu_indices(binned.shape[0], 1)
        b_signal[:, 0] = binned[indices[0], indices[1], 0]
        b_signal[:, 1] = binned[indices[0], indices[1], 1]
    else:
        if not silent:
            print >> sys.stderr, ("Unrecognized binned array type. No data returned.\n"),
        return None
    if not silent:
        print >> sys.stderr, ("Dynamically binning data..."),
    # Determine bin edges relative to unbinned positions
    unbinnedmids = (unbinnedpositions[:, 0] + unbinnedpositions[:, 1]) / 2
    binedges = numpy.zeros(binbounds.shape, dtype=numpy.int32)
    binedges[:, 0] = numpy.searchsorted(unbinnedmids, binbounds[:, 0])
    binedges[:, 1] = numpy.searchsorted(unbinnedmids, binbounds[:, 1])
    # Determine bin midpoints
    mids = (binbounds[:, 0] + binbounds[:, 1]) / 2
    # Dynamically bin using appropriate array type combination
    _fivec_binning.dynamically_bin_upper_from_upper(ub_signal, unbinnedmids, b_signal, binedges,
                                              mids, minobservations, searchdistance, int(removefailed))
    if binned_type == 'full':
        binned[indices[0], indices[1], 0] = b_signal[:, 0]
        binned[indices[0], indices[1], 1] = b_signal[:, 1]
        binned[indices[1], indices[0], 0] = b_signal[:, 0]
        binned[indices[1], indices[0], 1] = b_signal[:, 1]
    if not silent:
        print >> sys.stderr, ("Done\n"),
    return None


def find_trans_signal(fivec, region1, region2, binsize=0, binbounds1=None, start1=None, stop1=None, startfrag1=None,
                    stopfrag1=None, binbounds2=None, start2=None, stop2=None, startfrag2=None, stopfrag2=None,
                    datatype='enrichment', arraytype='compact', skipfiltered=False, returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    :param fivec:  A :class:`FiveC <hifive.fivec.FiveC>` class object containing fragment and count data.
    :type fivec: :class:`FiveC <hifive.fivec.FiveC>`
    :param region1: The index of the first region to pull data from.
    :type region1: int.
    :param region2: The index of the second region to pull data from.
    :type region2: int.
    :param binsize: This is the coordinate width of each bin. A value of zero indicates unbinned. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds1: An array containing start and stop coordinates for a set of user-defined bins for region1. Any fragment not falling in a bin is ignored.
    :type binbounds1: numpy array
    :param start1: The smallest coordinate to include in the array from 'region1', measured from fragment midpoints. If 'binbounds1' is given, this value is ignored. If both 'start1' and 'startfrag1' are given, 'start1' will override 'startfrag1'. If unspecified, this will be set to the midpoint of the first fragment for 'region1', adjusted to the first multiple of 'binsize' if not zero. Optional.
    :type start1: int.
    :param stop1: The largest coordinate to include in the array from 'region1', measured from fragment midpoints. If 'binbounds1' is given, this value is ignored. If both 'stop1' and 'stopfrag1' are given, 'stop1' will override 'stopfrag1'. If unspecified, this will be set to the midpoint of the last fragment plus one for 'region1', adjusted to the last multiple of 'start1' + 'binsize' if not zero. Optional.
    :type stop1: int.
    :param startfrag1: The first fragment to include in the array from 'region1'. If 'binbounds1' is given, this value is ignored. If unspecified and 'start1' is not given, this is set to the first fragment in 'region1'. In cases where 'start1' is specified and conflicts with 'startfrag1', 'start1' is given preference. Optional.
    :type startfrag1: int.
    :param stopfrag1: The first fragment not to include in the array from 'region1'. If 'binbounds1' is given, this value is ignored. If unspecified and 'stop1' is not given, this is set to the last fragment in 'region1' plus one. In cases where 'stop1' is specified and conflicts with 'stopfrag1', 'stop1' is given preference. Optional.
    :type stopfrag1: int.
    :param binbounds2: An array containing start and stop coordinates for a set of user-defined bins for region2. Any fragment not falling in a bin is ignored.
    :type binbounds2: numpy array
    :param start2: The smallest coordinate to include in the array from 'region2', measured from fragment midpoints. If 'binbounds2' is given, this value is ignored. If both 'start2' and 'startfrag2' are given, 'start2' will override 'startfrag2'. If unspecified, this will be set to the midpoint of the first fragment for 'region2', adjusted to the first multiple of 'binsize' if not zero. Optional.
    :type start2: int.
    :param stop2: The largest coordinate to include in the array from 'region2', measured from fragment midpoints. If 'binbounds2' is given, this value is ignored. If both 'stop2' and 'stopfrag2' are given, 'stop2' will override 'stopfrag2'. If unspecified, this will be set to the midpoint of the last fragment plus one for 'region2', adjusted to the last multiple of 'start2' + 'binsize' if not zero. Optional.
    :type stop2: int.
    :param startfrag2: The first fragment to include in the array from 'region2'. If 'binbounds2' is given, this value is ignored. If unspecified and 'start2' is not given, this is set to the first fragment in 'region2'. In cases where 'start2' is specified and conflicts with 'startfrag2', 'start2' is given preference. Optional.
    :type startfrag2: int.
    :param stopfrag2: The first fragment not to include in the array from 'region2'. If 'binbounds2' is given, this value is ignored. If unspecified and 'stop2' is not given, this is set to the last fragment in 'region2' plus one. In cases where 'stop2' is specified and conflicts with 'stopfrag2', 'stop2' is given preference. Optional.
    :type stopfrag2: int.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment', 'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fragments return value of one. Expected values are returned for 'distance', 'fragment', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both correction and distance mean values. 'enrichment' also scales both observed and expected by the standard deviation, giving a completely normalized set of values.
    :type datatype: str.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact' (though only when 'binned' is zero), and 'full'. 'compact' means data are arranged in a N x M x 2 array where N and M are the number of forward and reverse probe fragments, respectively. This will only return the array of forward primers from 'region1' and reverse primers from 'region2'. 'full' returns a square, symmetric array of size N x N x 2 where N is the total number of fragments.
    :type arraytype: str.
    :param skipfiltered: If 'True', all interaction bins for filtered out fragments are removed and a reduced-size array is returned.
    :type skipfiltered: bool.
    :param returnmapping: If 'True', a list containing the data array and either one or four 2d arrays containing first coordinate included and excluded from each bin, and the first fragment included and excluded from each bin  corresponding to both axes or the first and second axis for 'region1' forward fragments by 'region2' reverse fragments and 'region1' reverse fragments by 'region2' forward fragments for a full or compact array, respectively, is returned. Otherwise only the data array (or data arrays is compact) is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that all values are acceptable
    if datatype not in ['raw', 'fragment', 'distance', 'enrichment', 'expected']:
        if not silent:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    elif datatype in ['fragment', 'enrichment'] and fivec.normalization == 'none':
        if not silent:
            print >> sys.stderr, ("Normalization has not been performed yet on this project. Select either 'raw' or 'distance' for datatype. No data returned\n"),
        return None
    elif datatype in ['distance', 'enrichment'] and fivec.trans_mean is None:
        fivec.find_trans_mean()
    if arraytype not in ['full', 'compact']:
        if not silent:
            print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    if arraytype == 'compact' and (binsize > 0 or not binbounds1 is None or not binbounds2 is None):
        if not silent:
            print >> sys.stderr, ("'Compact' array can only be used with unbinned data. No data returned.\n"),
        return None
    # Determine start, stop, startfrag, and stopfrag
    chrom1 = fivec.frags['regions']['chromosome'][region1]
    chrint1 = fivec.chr2int[chrom1]
    if not binbounds1 is None:
        start1 = binbounds1[0, 0]
        stop1 = binbounds1[-1, 1]
        startfrag1 = _find_frag_from_coord(fivec, chrint1, start1)
        stopfrag1 = _find_frag_from_coord(fivec, chrint1, stop1)
    else:
        if start1 is None and startfrag1 is None:
            startfrag1 = fivec.frags['regions']['start_frag'][region1]
            start1 = fivec.frags['fragments']['mid'][startfrag1]
            if binsize > 0:
                start1 = (start1 / binsize) * binsize
        elif start1 is None:
            start1 = fivec.frags['fragments']['mid'][startfrag1]
            if binsize > 0:
                start1 = (start1 / binsize) * binsize
        else:
            startfrag1 = _find_frag_from_coord(fivec, chrint1, start1)
        if (stop1 is None or stop1 == 0) and stopfrag1 is None:
            stopfrag1 = fivec.frags['regions']['stop_frag'][region1]
            stop1 = fivec.frags['fragments']['mid'][stopfrag1 - 1] + 1
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        elif stop1 is None or stop1 == 0:
            stop1 = fivec.frags['fragments']['mid'][stopfrag1 - 1] + 1
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        else:
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
            stopfrag1 = _find_frag_from_coord(fivec, chrint1, stop1)
    chrom2 = fivec.frags['regions']['chromosome'][region2]
    chrint2 = fivec.chr2int[chrom2]
    if not binbounds2 is None:
        start2 = binbounds2[0, 0]
        stop2 = binbounds2[-1, 1]
        startfrag2 = _find_frag_from_coord(fivec, chrint2, start2)
        stopfrag2 = _find_frag_from_coord(fivec, chrint2, stop2)
    else:
        if start2 is None and startfrag2 is None:
            startfrag2 = fivec.frags['regions']['start_frag'][region2]
            start2 = fivec.frags['fragments']['mid'][startfrag2]
            if binsize > 0:
                start2 = (start2 / binsize) * binsize
        elif start2 is None:
            start2 = fivec.frags['fragments']['mid'][startfrag2]
            if binsize > 0:
                start2 = (start2 / binsize) * binsize
        else:
            startfrag2 = _find_frag_from_coord(fivec, chrint2, start2)
        if (stop2 is None or stop2 == 0) and stopfrag2 is None:
            stopfrag2 = fivec.frags['regions']['stop_frag'][region2]
            stop2 = fivec.frags['fragments']['mid'][stopfrag2 - 1] + 1
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        elif stop2 is None or stop2 == 0:
            stop2 = fivec.frags['fragments']['mid'][stopfrag2 - 1] + 1
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        else:
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
            stopfrag2 = _find_frag_from_coord(fivec, chrint2, stop2)
    if stopfrag1 - startfrag1 == 0 or stopfrag2 - startfrag2 == 0:
        if not silent:
            print >> sys.stderr, ("Insufficient data, no data returned.\n"),
        return None
    if not silent:
        print >> sys.stderr, ("Finding %s %s array for %s:%i-%i by %s:%i-%i...") % (datatype, arraytype, chrom1,
                                                                                    start1, stop1, chrom2, start2,
                                                                                    stop2),
    # Determine mapping of valid fends to bins
    mids1 = fivec.frags['fragments']['mid'][startfrag1:stopfrag1]
    mids2 = fivec.frags['fragments']['mid'][startfrag2:stopfrag2]
    strands1 = fivec.frags['fragments']['strand'][startfrag1:stopfrag1]
    strands2 = fivec.frags['fragments']['strand'][startfrag2:stopfrag2]
    mapping1 = numpy.zeros(stopfrag1 - startfrag1, dtype=numpy.int32) - 1
    mapping2 = numpy.zeros(stopfrag2 - startfrag2, dtype=numpy.int32) - 1
    valid1 = numpy.where(fivec.filter[startfrag1:stopfrag1] > 0)[0]
    valid2 = numpy.where(fivec.filter[startfrag2:stopfrag2] > 0)[0]
    if binsize == 0 and binbounds1 is None and binbounds2 is None:
        if arraytype == 'compact':
            forward = numpy.where(strands1 == 0)[0]
            valid1 = numpy.where((fivec.filter[startfrag1:stopfrag1] > 0) * (strands1 == 0))[0]
            for_valid = numpy.where(fivec.filter[startfrag1 + forward] > 0)[0]
            reverse = numpy.where(strands2 == 1)[0]
            valid2 = numpy.where((fivec.filter[startfrag2:stopfrag2] > 0) * (strands2 == 1))[0]
            rev_valid = numpy.where(fivec.filter[startfrag2 + reverse] > 0)[0]
            if skipfiltered:
                mapping1[valid1] = numpy.arange(valid1.shape[0])
                num_bins1 = valid1.shape[0]
                mapping2[valid2] = numpy.arange(valid2.shape[0])
                num_bins2 = valid2.shape[0]
            else:
                num_bins1 = forward.shape[0]
                mapping1[valid1] = for_valid
                num_bins2 = reverse.shape[0]
                mapping2[valid2] = rev_valid
        else:
            if skipfiltered:
                mapping1[valid1] = numpy.arange(valid1.shape[0])
                num_bins1 = valid1.shape[0]
                mapping2[valid2] = numpy.arange(valid2.shape[0])
                num_bins2 = valid2.shape[0]
            else:
                mapping1[valid1] = valid1
                num_bins1 = mapping1.shape[0]
                mapping2[valid2] = valid2
                num_bins2 = mapping2.shape[0]
    else:
        if not binbounds1 is None:
            start_indices = numpy.searchsorted(binbounds1[:, 0], mids1[valid1], side='right') - 1
            stop_indices = numpy.searchsorted(binbounds1[:, 1], mids1[valid1], side='right')
            where = numpy.where(start_indices == stop_indices)[0]
            mapping1[valid1[where]] = start_indices[where]
            num_bins1 = binbounds1.shape[0]
        else:
            mapping1[valid1] = (mids1[valid1] - start1) / binsize
            num_bins1 = (stop1 - start1) / binsize
        if not binbounds2 is None:
            start_indices = numpy.searchsorted(binbounds2[:, 0], mids2[valid2], side='right') - 1
            stop_indices = numpy.searchsorted(binbounds2[:, 1], mids2[valid2], side='right')
            where = numpy.where(start_indices == stop_indices)[0]
            mapping2[valid2[where]] = start_indices[where]
            num_bins2 = binbounds2.shape[0]
        else:
            mapping2[valid2] = (mids2[valid2] - start2) / binsize
            num_bins2 = (stop2 - start2) / binsize
    if num_bins1 < 1 or num_bins2 < 1:
        if not silent:
            print >> sys.stderr, ("Insufficient data\n"),
        return None
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        if startfrag1 < startfrag2:
            start_index = fivec.data['trans_indices'][startfrag1]
            stop_index = fivec.data['trans_indices'][stopfrag1]
        else:
            start_index = fivec.data['trans_indices'][startfrag2]
            stop_index = fivec.data['trans_indices'][stopfrag2]
        if stop_index - start_index == 0:
            if not silent:
                print >> sys.stderr, ("Insufficient data, no data returned.\n"),
            return None
        if startfrag1 < startfrag2:
            data_indices = fivec.data['trans_indices'][startfrag1:(stopfrag1 + 1)]
        else:
            data_indices = fivec.data['trans_indices'][startfrag2:(stopfrag2 + 1)]
        data_indices -= data_indices[0]
        data = fivec.data['trans_data'][start_index:stop_index, :]
        if startfrag1 < startfrag2:
            data[:, 0] -= startfrag1
            data[:, 1] -= startfrag2
        else:
            data[:, 0] -= startfrag2
            data[:, 1] -= startfrag1
    else:
        data_indices = None
        data = None
    # If correction is required, determine what type and get appropriate data
    corrections1 = None
    corrections2 = None
    binning_corrections = None
    correction_indices = None
    binning_num_bins = None
    frag_indices = None
    if datatype in ['fragment', 'enrichment', 'expected']:
        if fivec.normalization.count('probability') + fivec.normalization.count('express') > 0:
            corrections1 = fivec.corrections[startfrag1:stopfrag1]
            corrections2 = fivec.corrections[startfrag2:stopfrag2]
        if fivec.normalization.count('binning') > 0:
            binning_corrections = fivec.binning_corrections
            correction_indices = fivec.binning_correction_indices
            binning_num_bins = fivec.binning_num_bins
            frag_indices = fivec.binning_frag_indices
    if datatype in ['distance', 'enrichment', 'expected']:
        if fivec.trans_mean is None:
            fivec.find_trans_mean()
        trans_mean = fivec.trans_mean
    else:
        trans_mean = 0.0
    # Create requested array and corresponding mapping
    if startfrag1 < startfrag2:
        data_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins2, num_bins1, 2), dtype=numpy.float32)
    # Fill in data values
    if datatype != 'expected':
        if startfrag1 < startfrag2:
            _fivec_binning.find_trans_observed(data, data_indices, mapping1, mapping2, data_array)
        else:
            _fivec_binning.find_trans_observed(data, data_indices, mapping2, mapping1, data_array)
    if datatype != 'raw':
        if startfrag1 < startfrag2:
            _fivec_binning.find_trans_expected(mapping1, mapping2, corrections1, corrections2, binning_corrections,
                                               correction_indices, binning_num_bins, frag_indices, strands1, strands2,
                                               data_array, trans_mean, startfrag1, startfrag2)
        else:
            _fivec_binning.find_trans_expected(mapping2, mapping1, corrections2, corrections1, binning_corrections,
                                               correction_indices, binning_num_bins, frag_indices, strands2, strands1,
                                               data_array, trans_mean, startfrag2, startfrag1)
    else:
        where = numpy.where(data_array[:, :, 0] > 0.0)
        data_array[where[0], where[1], 1] = 1.0
    if datatype == 'expected':
        where = numpy.where(data_array[:, :, 1] > 0.0)
        data_array[where[0], where[1], 0] = 1.0
    # if startfrag2 < startfrag1, transpose data_array
    if startfrag1 > startfrag2:
        data_array = numpy.transpose(data_array, (1, 0, 2))
    if returnmapping:
        if arraytype == 'compact':
            bin_mapping1 = numpy.zeros((data_array.shape[0], 4), dtype=numpy.int32)
            bin_mapping2 = numpy.zeros((data_array.shape[1], 4), dtype=numpy.int32)
            if skipfiltered:
                bin_mapping1[:, 2] = valid1 + startfrag1
                bin_mapping2[:, 2] = valid2 + startfrag2
            else:
                bin_mapping1[:, 2] = numpy.where(strands1 == 0)[0] + startfrag1
                bin_mapping2[:, 2] = numpy.where(strands2 == 1)[0] + startfrag2
            bin_mapping1[:, 3] = bin_mapping1[:, 2] + 1
            bin_mapping1[:, 0] = fivec.frags['fragments']['start'][bin_mapping1[:, 2]]
            bin_mapping1[:, 1] = fivec.frags['fragments']['stop'][bin_mapping1[:, 2]]
            bin_mapping2[:, 3] = bin_mapping2[:, 2] + 1
            bin_mapping2[:, 0] = fivec.frags['fragments']['start'][bin_mapping2[:, 2]]
            bin_mapping2[:, 1] = fivec.frags['fragments']['stop'][bin_mapping2[:, 2]]
            if not silent:
                print >> sys.stderr, ("Done\n"),
            return [data_array, bin_mapping1, bin_mapping2]
        else:
            bin_mapping1 = numpy.zeros((num_bins1, 4), dtype=numpy.int32)
            bin_mapping2 = numpy.zeros((num_bins2, 4), dtype=numpy.int32)
            if binsize == 0 and binbounds1 is None:
                if skipfiltered:
                    bin_mapping1[:, 2] = valid1 + startfrag1
                else:
                    bin_mapping1[:, 2] = numpy.arange(startfrag1, stopfrag1)
                bin_mapping1[:, 3] = bin_mapping1[:, 2] + 1
                bin_mapping1[:, 0] = fivec.frags['fragments']['start'][bin_mapping1[:, 2]]
                bin_mapping1[:, 1] = fivec.frags['fragments']['stop'][bin_mapping1[:, 2]]
            else:
                if binbounds1 is None:
                    bin_mapping1[:, 0] = start1 + binsize * numpy.arange(num_bins1)
                    bin_mapping1[:, 1] = bin_mapping1[:, 0] + binsize
                else:
                    bin_mapping1[:, :2] = binbounds1
                bin_mapping1[:, 2] = numpy.searchsorted(mids1, bin_mapping1[:, 0]) + startfrag1
                bin_mapping1[:, 3] = numpy.searchsorted(mids1, bin_mapping1[:, 1]) + startfrag1
            if binsize == 0 and binbounds2 is None:
                if skipfiltered:
                    bin_mapping2[:, 2] = valid2 + startfrag2
                else:
                    bin_mapping2[:, 2] = numpy.arange(startfrag2, stopfrag2)
                bin_mapping2[:, 3] = bin_mapping2[:, 2] + 1
                bin_mapping2[:, 0] = fivec.frags['fragments']['start'][bin_mapping2[:, 2]]
                bin_mapping2[:, 1] = fivec.frags['fragments']['stop'][bin_mapping2[:, 2]]
            else:
                if binbounds2 is None:
                    bin_mapping2[:, 0] = start2 + binsize * numpy.arange(num_bins2)
                    bin_mapping2[:, 1] = bin_mapping2[:, 0] + binsize
                else:
                    bin_mapping2[:, :2] = binbounds2
                bin_mapping2[:, 2] = numpy.searchsorted(mids2, bin_mapping2[:, 0]) + startfrag2
                bin_mapping2[:, 3] = numpy.searchsorted(mids2, bin_mapping2[:, 1]) + startfrag2
            if not silent:
                print >> sys.stderr, ("Done\n"),
            return [data_array, bin_mapping1, bin_mapping2]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return data_array


def dynamically_bin_trans_array(unbinned, unbinnedpositions1, unbinnedpositions2, binned, binbounds1, binbounds2,
                                minobservations=50, searchdistance=0, removefailed=True, **kwargs):
    """
    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations', or 'searchdistance' criteria.

    :param unbinned: A full array containing data to be binned.
    :type unbinned: numpy array
    :param unbinnedpositions1: A 2d integer array indicating the first and last coordinate of each bin in 'unbinned' array along the first axis.
    :type unbinnedpositions1: numpy array
    :param unbinnedpositions2: A 2d integer array indicating the first and last coordinate of each bin in 'unbinned' array along the second axis.
    :type unbinnedpositions2: numpy array
    :param binned: A full array containing binned data to be dynamically binned. Data in this array will be altered by this function.
    :type binned: numpy array
    :param binbounds1: A N x 2 integer array indicating the start and end position of each of N bins in 'binned' array along the first axis.
    :type binbounds1: numpy array
    :param binbounds2: A N x 2 integer array indicating the start and end position of each of N bins in 'binned' array along the second axis.
    :type binbounds2: numpy array
    :param minobservations: The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
    :type minobservations: int.
    :param searchdistance: The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on expansion distance.
    :type searchdistance: int.
    :param removefailed: If a non-zero 'searchdistance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'removefailed' is True, the observed and expected values for that bin are zero.
    :type removefailed: bool.
    :returns: None
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if not silent:
        print >> sys.stderr, ("Dynamically binning data..."),
    # Determine bin edges relative to unbinned positions
    unbinnedmids1 = (unbinnedpositions1[:, 0] + unbinnedpositions1[:, 1]) / 2
    unbinnedmids2 = (unbinnedpositions2[:, 0] + unbinnedpositions2[:, 1]) / 2
    binedges1 = numpy.zeros(binbounds1.shape, dtype=numpy.int32)
    binedges1[:, 0] = numpy.searchsorted(unbinnedmids1, binbounds1[:, 0])
    binedges1[:, 1] = numpy.searchsorted(unbinnedmids1, binbounds1[:, 1])
    binedges2 = numpy.zeros(binbounds2.shape, dtype=numpy.int32)
    binedges2[:, 0] = numpy.searchsorted(unbinnedmids2, binbounds2[:, 0])
    binedges2[:, 1] = numpy.searchsorted(unbinnedmids2, binbounds2[:, 1])
    # Determine bin midpoints
    mids1 = (binbounds1[:, 0] + binbounds1[:, 1]) / 2
    mids2 = (binbounds2[:, 0] + binbounds2[:, 1]) / 2
    # Dynamically bin using appropriate array type combination
    _fivec_binning.dynamically_bin_trans(unbinned, unbinnedmids1, unbinnedmids2, binned, binedges1,
                                         binedges2, mids1, mids2, minobservations, searchdistance, int(removefailed))
    if not silent:
        print >> sys.stderr, ("Done\n"),
    return None


def write_heatmap_dict(fivec, filename, binsize, includetrans=True, datatype='enrichment',
                       regions=[], arraytype='full', dynamically_binned=False, minobservations=0,
                       searchdistance=0, expansion_binsize=0, removefailed=False, **kwargs):
    """
    Create an h5dict file containing binned interaction arrays, bin positions, and an index of included regions.

    :param fivec: A :class:`FiveC <hifive.fivec.FiveC>` class object containing fragment and count data.
    :type fivec: :class:`FiveC <hifive.fivec.FiveC>`
    :param filename: Location to write h5dict object to.
    :type filename: str.
    :param binsize: Size of bins for interaction arrays. If "binsize" is zero, fragment interactions are returned without binning.
    :type binsize: int.
    :param includetrans: Indicates whether trans interaction arrays should be calculated and saved.
    :type includetrans: bool.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment', 'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-filtered bins return value of 1. Expected values are returned for 'distance', 'fragment', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
        :type datatype: str.
    :param arraytype: This determines what shape of array data are returned in if unbinned heatmaps are requested. Acceptable values are 'compact' and 'full'. 'compact' means data are arranged in a N x M array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N.
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
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if binsize > 0:
        arraytype='full'
    elif dynamically_binned:
        arraytype='full'
    # Check if trans mean is needed and calculate if not already done
    if includetrans and datatype in ['distance', 'enrichment'] and 'trans_mean' not in fivec.__dict__.keys():
        fivec.find_trans_mean()
    # Check if filename already exists, and remove if it does
    if os.path.exists(filename):
        if not silent:
            print >> sys.stderr, ("%s already exists, overwriting.") % filename
        subprocess.call('rm %s' % filename, shell=True)
    if not silent:
        print >> sys.stderr, ("Creating binned heatmap...\n"),
    output = h5py.File(filename, 'w')
    # If regions not specified, fill list
    if len(regions) == 0:
        regions = list(numpy.arange(fivec.frags['regions'].shape[0]))
    if binsize > 0:
        output.attrs['resolution'] = binsize
    else:
        output.attrs['resolution'] = 'fragment'
    # Find cis heatmaps
    remove = []
    for region in regions:
        if dynamically_binned:
            temp = find_cis_signal(fivec, region, binsize=expansion_binsize, datatype=datatype, arraytype='full',
                                   returnmapping=True, silent=silent, skipfiltered=True)
            if temp is None:
                remove.append(region)
                continue
            expansion, exp_mapping = temp
            binned, mapping = find_cis_signal(fivec, region, binsize=binsize, datatype=datatype, arraytype=arraytype,
                                              returnmapping=True, silent=silent, skipfiltered=True)

            dynamically_bin_cis_array(expansion, exp_mapping, binned, mapping, minobservations=minobservations,
                                      searchdistance=searchdistance, removefailed=removefailed, silent=silent)
            results = [binned, mapping]
        else:
            results = find_cis_signal(fivec, region, binsize=binsize, datatype=datatype, arraytype=arraytype,
                                      returnmapping=True, silent=silent, skipfiltered=True)
        # Check if array contains data
        if results is None or results[0].shape[0] == 0:
            remove.append(region)
            continue
        output.create_dataset('%i.counts' % region, data=results[0][:, :, 0])
        output.create_dataset('%i.expected' % region, data=results[0][:, :, 1])
        if binsize > 0 or arraytype == 'full':
            output.create_dataset('%i.positions' % region, data=results[1][:, :2])
        else:
            output.create_dataset('%i.forward_positions' % region, data=results[1][:, :2])
            output.create_dataset('%i.reverse_positions' % region, data=results[2][:, :2])
    for region in remove:
        del regions[regions.index(region)]
    all_regions = fivec.frags['regions'][...]
    output.create_dataset('regions', data=all_regions[regions][...])
    # If requested, find trans heatmaps
    if includetrans:
        for i in range(len(regions)-1):
            for j in range(i + 1, len(regions)):
                if dynamically_binned:
                    expansion, exp_map1, exp_map2 = find_trans_signal(fivec, regions[i], regions[j],
                                                                      binsize=expansion_binsize, datatype=datatype,
                                                                      arraytype='full', skipfiltered=True,
                                                                      silent=silent)
                    if arraytype == 'compact':
                        binned, mapping1, mapping2 = find_trans_signal(fivec, regions[i], regions[j], binsize=binsize,
                                                                       datatype=datatype, arraytype=arraytype,
                                                                       skipfiltered=True, silent=silent)
                        dynamically_bin_trans_array(expansion, exp_map1, exp_map2, binned, mapping1, mapping2,
                                                    minobservations=minobservations, searchdistance=searchdistance,
                                                    removefailed=removefailed, silent=silent)
                        output.create_dataset('%s_by_%s.counts' % (regions[i], regions[j]), data=binned[:, :, 0])
                        output.create_dataset('%s_by_%s.expected' % (regions[i], regions[j]), data=binned[:, :, 1])
                        binned, mapping1, mapping2 = find_trans_signal(fivec, regions[j], regions[i], binsize=binsize,
                                                                       datatype=datatype, arraytype=arraytype,
                                                                       skipfiltered=True, silent=silent)
                        dynamically_bin_trans_array(expansion, exp_map1, exp_map2, binned, mapping1, mapping2,
                                                    minobservations=minobservations, searchdistance=searchdistance,
                                                    removefailed=removefailed, silent=silent)
                        output.create_dataset('%s_by_%s.counts' % (regions[j], regions[i]), data=binned[:, :, 0])
                        output.create_dataset('%s_by_%s.expected' % (regions[j], regions[i]), data=binned[:, :, 1])
                    else:
                        binned, mapping1, mapping2 = find_trans_signal(fivec, regions[i], regions[j], binsize=binsize,
                                                                       datatype=datatype, arraytype=arraytype,
                                                                       skipfiltered=True, silent=silent)
                        dynamically_bin_trans_array(expansion, exp_map1, exp_map2, binned, mapping1, mapping2,
                                                    minobservations=minobservations, searchdistance=searchdistance,
                                                    removefailed=removefailed, silent=silent)
                        output.create_dataset('%s_by_%s.counts' % (regions[i], regions[j]), data=results[:, :, 0])
                        output.create_dataset('%s_by_%s.expected' % (regions[i], regions[j]), data=results[:, :, 1])
                else:   
                    if arraytype == 'compact':
                        results = find_trans_signal(fivec, regions[i], regions[j], binsize=binsize, datatype=datatype,
                                                    arraytype=arraytype, skipfiltered=True, silent=silent)
                        output.create_dataset('%s_by_%s.counts' % (regions[i], regions[j]), data=results[:, :, 0])
                        output.create_dataset('%s_by_%s.expected' % (regions[i], regions[j]), data=results[:, :, 1])
                        results = find_trans_signal(fivec, regions[j], regions[i], binsize=binsize, datatype=datatype,
                                                    arraytype=arraytype, skipfiltered=True, silent=silent)
                        output.create_dataset('%s_by_%s.counts' % (regions[j], regions[i]), data=results[:, :, 0])
                        output.create_dataset('%s_by_%s.expected' % (regions[j], regions[i]), data=results[:, :, 1])
                    else:
                        results = find_trans_signal(fivec, regions[i], regions[j], binsize=binsize, datatype=datatype,
                                                    arraytype=arraytype, skipfiltered=True, silent=silent)
                        output.create_dataset('%s_by_%s.counts' % (regions[i], regions[j]), data=results[:, :, 0])
                        output.create_dataset('%s_by_%s.expected' % (regions[i], regions[j]), data=results[:, :, 1])
    if 'history' in kwargs:
        output.attrs['history'] = kwargs['history']
    output.attrs['filetype'] = 'fivec_heatmap'
    output.close()
    if not silent:
        print >> sys.stderr, ("Creating binned heatmap...Done\n"),
    return None
