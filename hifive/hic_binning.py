#!/usr/bin/env python

"""
This is a module contains scripts for generating compact, upper-triangle and full matrices of HiC interaction data.

Concepts
--------

These functions rely on the :class:`HiC` class in conjunction with the :class:`Fend` and :class:`HiCData` classes.

Data can either be arranged in compact, complete, or flattened (row-major) upper-triangle arrays. Compact arrays are N x M, where N is the number of fends or bins, and M is the maximum distance between fends or bins. This is useful for working with sets of short interactions. Data can be raw, fend-corrected, distance-dependence removed, or enrichment values. Arrays are 3-dimensional with observed values in the first layer of d3, expected values in the second layer of d3. The exception to this is upper-triangle arrays, which are 2d, divinding observed and expected along the second axis.

API Documentation
-----------------
"""

import os
import sys
import subprocess

import numpy
import h5py
try:
    from mpi4py import MPI
except:
    pass

from libraries._hic_interactions import find_max_fend
import libraries._hic_binning as _hic_binning


def find_cis_signal(hic, chrom, binsize=10000, binbounds=None, start=None, stop=None, startfend=None, stopfend=None,
                    datatype='enrichment', arraytype='compact', maxdistance=0, skipfiltered=False, returnmapping=False,
                    **kwargs):
    """
    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param chrom: The name of a chromosome contained in 'hic'.
    :type chrom: str.
    :param binsize: This is the coordinate width of each bin. A value of zero indicates unbinned. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fend not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The smallest coordinate to include in the array, measured from fend midpoints or the start of the first bin. If 'binbounds' is given, this value is ignored. If both 'start' and 'startfend' are given, 'start' will override 'startfend'. If unspecified, this will be set to the midpoint of the first fend for 'chrom', adjusted to the first multiple of 'binsize' if not zero. Optional.
    :type start: int.
    :param stop: The largest coordinate to include in the array, measured from fend midpoints or the end of the last bin. If 'binbounds' is given, this value is ignored. If both 'stop' and 'stopfend' are given, 'stop' will override 'stopfend'. If unspecified, this will be set to the midpoint of the last fend plus one for 'chrom', adjusted to the last multiple of 'start' + 'binsize' if not zero. Optional.
    :type stop: int.
    :param startfend: The first fend to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'start' is not given, this is set to the first valid fend in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given preference. Optional
    :type startfend: int.
    :param stopfend: The first fend not to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'stop' is not given, this is set to the last valid fend in 'chrom' plus one. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given preference. Optional.
    :type stopfend: str.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
    :type datatype: str.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    :type arraytype: str.
    :param maxdistance: This specifies the maximum coordinate distance between bins that will be included in the array. If set to zero, all distances are included.
    :type maxdistance: str.
    :param skipfiltered: If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
    :type skipfiltered: bool.
    :param returnmapping: If 'True', a list containing the data array and a 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that all values are acceptable
    if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
        if not silent:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    elif datatype in ['fend', 'enrichment'] and hic.normalization == 'none':
        if not silent:
            print >> sys.stderr, ("Normalization has not been performed yet on this project. Select either 'raw' or 'distance' for datatype. No data returned\n"),
        return None
    elif datatype in ['distance', 'enrichment'] and hic.distance_parameters is None:
        if not silent:
            print >> sys.stderr, ("No distance-dependence relationship has been calculated for this project yet. Select either 'raw' or 'fend' for datatype. No data returned\n"),
        return None
    if arraytype not in ['full', 'compact', 'upper']:
        if not silent:
            print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine start, stop, startfend, and stopfend
    chrint = hic.chr2int[chrom.strip('chr')]
    if not binbounds is None:
        start = binbounds[0, 0]
        stop = binbounds[-1, 1]
        startfend = _find_fend_from_coord(hic, chrint, start)
        stopfend = _find_fend_from_coord(hic, chrint, stop)
    else:
        if start is None and startfend is None:
            startfend = hic.fends['chr_indices'][chrint]
            while startfend < hic.fends['chr_indices'][chrint + 1] and hic.filter[startfend] == 0:
                startfend += 1
            if startfend == hic.fends['chr_indices'][chrint + 1]:
                if not silent:
                    print >> sys.stderr, ("Insufficient data.\n"),
                return None
            start = hic.fends['fends']['mid'][startfend]
            if binsize > 0:
                start = (start / binsize) * binsize
        elif start is None:
            start = hic.fends['fends']['mid'][startfend]
            if binsize > 0:
                start = (start / binsize) * binsize
        else:
            startfend = _find_fend_from_coord(hic, chrint, start)
        if (stop is None or stop == 0) and stopfend is None:
            stopfend = hic.fends['chr_indices'][chrint + 1]
            while stopfend > hic.fends['chr_indices'][chrint] and hic.filter[stopfend - 1] == 0:
                stopfend -= 1
            stop = hic.fends['fends']['mid'][stopfend - 1]
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        elif stop is None or stop == 0:
            stop = hic.fends['fends']['mid'][stopfend - 1]
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
        else:
            if binsize > 0:
                stop = ((stop - 1 - start) / binsize + 1) * binsize + start
            stopfend = _find_fend_from_coord(hic, chrint, stop)
    if not silent:
        print >> sys.stderr, ("Finding %s %s array for %s:%i-%i...") % (datatype, arraytype, chrom, start, stop),
    # If datatype is not 'expected', pull the needed slice of data
    if datatype != 'expected':
        start_index = hic.data['cis_indices'][startfend]
        stop_index = hic.data['cis_indices'][stopfend]
        if start_index == stop_index:
            if not silent:
                print >> sys.stderr, ("Insufficient data\n"),
            return None
        data_indices = hic.data['cis_indices'][startfend:(stopfend + 1)]
        data_indices -= data_indices[0]
        data = hic.data['cis_data'][start_index:stop_index, :]
        data[:, :2] -= startfend
    else:
        data_indices = None
        data = None
    # Determine mapping of valid fends to bins
    mapping = numpy.zeros(stopfend - startfend, dtype=numpy.int32) - 1
    valid = numpy.where(hic.filter[startfend:stopfend] > 0)[0]
    mids = hic.fends['fends']['mid'][startfend:stopfend]
    if binsize == 0 and binbounds is None:
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
    # Find maximum interaction partner for each fend
    if num_bins < 2:
        if not silent:
            print >> sys.stderr, ("Insufficient data\n"),
        return None
    max_fend = numpy.zeros(mapping.shape[0], dtype=numpy.int32)
    find_max_fend(max_fend, mids, hic.fends['fends']['chr'][startfend:stopfend],
                  hic.fends['chr_indices'][...], startfend, maxdistance)
    max_fend = numpy.minimum(max_fend, mapping.shape[0])
    max_bin = numpy.amax(max_fend - numpy.arange(mapping.shape[0]))
    if max_bin <= 0:
        if not silent:
            print >> sys.stderr, ("Insufficient data.\n"),
        return None
    # If correction is required, determine what type and get appropriate data
    corrections = None
    binning_corrections = None
    correction_indices = None
    binning_num_bins = None
    fend_indices = None
    if datatype in ['fend', 'enrichment', 'expected']:
        if hic.normalization in ['express', 'probability', 'binning-express', 'binning-probability']:
            corrections = hic.corrections[startfend:stopfend]
        if hic.normalization in ['binning', 'binning-express', 'binning-probability']:
            binning_corrections = hic.binning_corrections
            correction_indices = hic.binning_correction_indices
            binning_num_bins = hic.binning_num_bins
            fend_indices = hic.binning_fend_indices
    if datatype in ['distance', 'enrichment', 'expected']:
        distance_parameters = hic.distance_parameters
        chrom_mean = hic.chromosome_means[chrint]
    else:
        distance_parameters = None
        chrom_mean = 0.0
    # Create requested array
    if arraytype == 'compact':
        data_array = numpy.zeros((num_bins, max_bin, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in data values
    if arraytype == 'compact':
        if datatype != 'raw':
            _hic_binning.find_cis_compact_expected(mapping, corrections, binning_corrections, correction_indices,
                                                   binning_num_bins, fend_indices, mids, distance_parameters,
                                                   max_fend, data_array, chrom_mean, startfend)
        if datatype != 'expected':
            _hic_binning.find_cis_compact_observed(data, data_indices, mapping, max_fend, data_array)
        else:
            where = numpy.where(data_array[:, :, 1] > 0.0)
            data_array[where[0], where[1], 0] = 1.0
        if datatype == 'raw':
            where = numpy.where(data_array[:, :, 0] > 0.0)
            data_array[where[0], where[1], 1] = 1.0
    else:
        if datatype != 'raw':
            _hic_binning.find_cis_upper_expected(mapping, corrections, binning_corrections, correction_indices,
                                                 binning_num_bins, fend_indices, mids, distance_parameters,
                                                 max_fend, data_array, chrom_mean, startfend)
        if datatype != 'expected':
            _hic_binning.find_cis_upper_observed(data, data_indices, mapping, max_fend, data_array)
        else:
            where = numpy.where(data_array[:, 1] > 0.0)[0]
            data_array[where, 0] = 1.0
        if datatype == 'raw':
            where = numpy.where(data_array[:, 0] > 0.0)[0]
            data_array[where, 1] = 1.0
    # If requesting 'full' array, convert 'upper' array type to 'full'
    if arraytype == 'full':
        indices = numpy.triu_indices(num_bins, 1)
        full_data_array = numpy.zeros((num_bins, num_bins, 2), dtype=numpy.float32)
        full_data_array[indices[1], indices[0], :] = data_array
        full_data_array[indices[0], indices[1], :] = data_array
        del data_array
        data_array = full_data_array
    if returnmapping:
        bin_mapping = numpy.zeros((num_bins, 4), dtype=numpy.int32)
        if binsize == 0 and binbounds is None:
            if skipfiltered:
                bin_mapping[:, 2] = valid + startfend
            else:
                bin_mapping[:, 2] = numpy.arange(startfend, stopfend)
            bin_mapping[:, 3] = bin_mapping[:, 2] + 1
            bin_mapping[:, 0] = hic.fends['fends']['start'][bin_mapping[:, 2]]
            bin_mapping[:, 1] = hic.fends['fends']['stop'][bin_mapping[:, 2]]
        else:
            if binbounds is None:
                bin_mapping[:, 0] = start + binsize * numpy.arange(num_bins)
                bin_mapping[:, 1] = bin_mapping[:, 0] + binsize
            else:
                bin_mapping[:, :2] = binbounds
            bin_mapping[:, 2] = numpy.searchsorted(mids, bin_mapping[:, 0]) + startfend
            bin_mapping[:, 3] = numpy.searchsorted(mids, bin_mapping[:, 1]) + startfend
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [data_array, bin_mapping]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return data_array

def _find_fend_from_coord(hic, chrint, coord):
    """Find the next fend after the coordinate on chromosome 'chrint'."""
    first_fend = hic.fends['chr_indices'][chrint]
    last_fend = hic.fends['chr_indices'][chrint + 1]
    return numpy.searchsorted(hic.fends['fends']['mid'][first_fend:last_fend], coord) + first_fend

def bin_cis_array(data_array, data_mapping, binsize=10000, binbounds=None, start=None, stop=None, arraytype='full',
                  returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill 'binsize' bins or bins defined by 'binbounds' with data provided in the array passed by 'data_array'.

    :param data_array: A 2d (upper) or 3d (compact) array containing data to be binned. Array format will be determined from the number of dimensions.
    :type data_array: numpy array
    :param data_mapping: An N x 4 2d integer array containing the start and stop coordinates, and start and stop fends for each of the N bin ranges in 'data_array'.
    :type data_mapping: numpy array
    :param binsize: This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any bin from 'data_array' not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The coordinate at the beginning of the first bin of the binned data. If unspecified, 'start' will be the first multiple of 'binsize' below the first coordinate from 'data_mapping'. If 'binbounds' is given, 'start' is ignored. Optional.
    :type start: int.
    :param stop: The coordinate at the end of the last bin of the binned data. If unspecified, 'stop' will be the first multiple of 'binsize' after the last coordinate from 'data_mapping'. If needed, 'stop' is adjusted upward to create a complete last bin. If 'binbounds' is given, 'stop' is ignored. Optional.
    :type stop: int.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    :type arraytype: str.
    :param returnmapping: If 'True', a list containing the data array and a 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing binned data requested with 'datatype' pulled from 'data_array' or list of binned data array and mapping array.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that arraytype value is acceptable
    if arraytype not in ['full', 'compact', 'upper']:
        if not silent:
            print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine input array type
    if len(data_array.shape) == 2 and data_mapping.shape[0] * (data_mapping.shape[0] - 1) / 2 == data_array.shape[0]:
        input_type = 'upper'
    elif len(data_array.shape) == 3 and data_array.shape[0] == data_mapping.shape[0]:
        input_type = 'compact'
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
    fend_ranges = numpy.zeros((binbounds.shape[0], 2), dtype=numpy.int32)
    for i in range(binbounds.shape[0]):
        firstbin = numpy.searchsorted(mids, binbounds[i, 0])
        lastbin = numpy.searchsorted(mids, binbounds[i, 1])
        mapping[firstbin:lastbin] = i
        fend_ranges[i, 0] = data_mapping[firstbin, 2]
        fend_ranges[i, 1] = data_mapping[lastbin, 3]
    # Create requested array
    if arraytype == 'compact':
        max_bin = (stop - start) / binsize + 1
        binned_array = numpy.zeros((num_bins, max_bin, 2), dtype=numpy.float32)
    else:
        binned_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in binned data values
    if arraytype == 'compact':
        if input_type == 'compact':
            _hic_binning.bin_compact_to_compact(binned_array, data_array, mapping)
        else:
            _hic_binning.bin_upper_to_compact(binned_array, data_array, mapping)
        # Trim unused bins
        valid = numpy.where(numpy.sum(binned_array[:, :, 1] > 0, axis=0) > 0)[0][-1]
        binned_array = binned_array[:, :(valid + 1), :]
    else:
        if input_type == 'compact':
            _hic_binning.bin_compact_to_upper(binned_array, data_array, mapping, num_bins)
        else:
            _hic_binning.bin_upper_to_upper(binned_array, data_array, mapping, num_bins)
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
        mapping[:, 2:4] = fend_ranges
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [binned_array, mapping]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return binned_array

def dynamically_bin_cis_array(unbinned, unbinnedpositions, binned, binbounds, minobservations=10,
                              searchdistance=0, removefailed=True, **kwargs):
    """
    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations', or 'searchdistance' criteria.

    :param unbinned: A 2d or 3d array containing data in either compact or upper format to be used for filling expanding bins. Array format will be determined from the number of dimensions.
    :type unbinned: numpy array
    :param unbinnedpositions: A 2d integer array indicating the first and last coordinate of each bin in 'unbinned' array.
    :type unbinnedpositions: numpy array
    :param binned: A 2d or 3d array containing binned data in either compact or upper format to be dynamically binned. Array format will be determined from the number of dimensions. Data in this array will be altered by this function.
    :type binned: numpy array
    :param binbounds: An integer array indicating the start and end position of each bin in 'binned' array. This array should be N x 2, where N is the number of intervals in 'binned'.
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
        unbinned_type = 'upper'
    elif len(unbinned.shape) == 3 and unbinned.shape[0] == unbinnedpositions.shape[0]:
        unbinned_type = 'compact'
    else:
        if not silent:
            print >> sys.stderr, ("Unrecognized unbinned array type. No data returned.\n"),
        return None
    # Determine binned array type
    if len(binned.shape) == 2 and binbounds.shape[0] * (binbounds.shape[0] - 1) / 2 == binned.shape[0]:
        binned_type = 'upper'
    elif len(binned.shape) == 3 and binned.shape[0] == binbounds.shape[0]:
        binned_type = 'compact'
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
    if unbinned_type == 'upper':
        if binned_type == 'upper':
            _hic_binning.dynamically_bin_upper_from_upper(unbinned, unbinnedmids, binned, binedges,
                                                      mids, minobservations, searchdistance, int(removefailed))
        else:
            _hic_binning.dynamically_bin_compact_from_upper(unbinned, unbinnedmids, binned, binedges,
                                                        mids, minobservations, searchdistance, int(removefailed))
    else:
        if binned_type == 'upper':
            _hic_binning.dynamically_bin_upper_from_compact(unbinned, unbinnedmids, binned, binedges,
                                                        mids, minobservations, searchdistance, int(removefailed))
        else:
            _hic_binning.dynamically_bin_compact_from_compact(unbinned, unbinnedmids, binned, binedges,
                                                          mids, minobservations, searchdistance, int(removefailed))
    if not silent:
        print >> sys.stderr, ("Done\n"),
    return None

def find_trans_signal(hic, chrom1, chrom2, binsize=10000, binbounds1=None, binbounds2=None, start1=None, stop1=None,
                      startfend1=None, stopfend1=None, start2=None, stop2=None, startfend2=None, stopfend2=None,
                      datatype='enrichment', skipfiltered=False, returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param chrom: The name of a chromosome contained in 'hic'.
    :type chrom: str.
    :param binsize: This is the coordinate width of each bin. A value of zero indicates unbinned. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds: An array containing start and stop coordinates for a set of user-defined bins. Any fend not falling in a bin is ignored.
    :type binbounds: numpy array
    :param start: The smallest coordinate to include in the array, measured from fend midpoints or the start of the first bin. If 'binbounds' is given, this value is ignored. If both 'start' and 'startfend' are given, 'start' will override 'startfend'. If unspecified, this will be set to the midpoint of the first fend for 'chrom', adjusted to the first multiple of 'binsize' if not zero. Optional.
    :type start: int.
    :param stop: The largest coordinate to include in the array, measured from fend midpoints or the end of the last bin. If 'binbounds' is given, this value is ignored. If both 'stop' and 'stopfend' are given, 'stop' will override 'stopfend'. If unspecified, this will be set to the midpoint of the last fend plus one for 'chrom', adjusted to the last multiple of 'start' + 'binsize' if not zero. Optional.
    :type stop: int.
    :param startfend: The first fend to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'start' is not given, this is set to the first valid fend in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given preference. Optional
    :type startfend: int.
    :param stopfend: The first fend not to include in the array. If 'binbounds' is given, this value is ignored. If unspecified and 'stop' is not given, this is set to the last valid fend in 'chrom' plus one. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given preference. Optional.
    :type stopfend: str.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
    :type datatype: str.
    :param arraytype: This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'. 'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    :type arraytype: str.
    :param maxdistance: This specifies the maximum coordinate distance between bins that will be included in the array. If set to zero, all distances are included.
    :type maxdistance: str.
    :param skipfiltered: If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
    :type skipfiltered: bool.
    :param returnmapping: If 'True', a list containing the data array and two 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin for the first and second axis is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing data requested with 'datatype'.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # check that all values are acceptable
    if datatype not in ['raw', 'fend', 'distance', 'enrichment', 'expected']:
        if not silent:
            print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    elif datatype in ['fend', 'enrichment'] and hic.normalization == 'none':
        if not silent:
            print >> sys.stderr, ("Normalization has not been performed yet on this project. Select either 'raw' or 'distance' for datatype. No data returned\n"),
        return None
    # Determine start, stop, startfend, and stopfend
    chrint1 = hic.chr2int[chrom1.strip('chr')]
    chrint2 = hic.chr2int[chrom2.strip('chr')]
    if not binbounds1 is None:
        start1 = binbounds1[0, 0]
        stop1 = binbounds1[-1, 1]
        startfend1 = _find_fend_from_coord(hic, chrint1, start1)
        stopfend1 = _find_fend_from_coord(hic, chrint1, stop1)
    else:
        if start1 is None and startfend1 is None:
            startfend1 = hic.fends['chr_indices'][chrint1]
            while startfend1 < hic.fends['chr_indices'][chrint1 + 1] and hic.filter[startfend1] == 0:
                startfend1 += 1
            if startfend1 == hic.fends['chr_indices'][chrint1 + 1]:
                if not silent:
                    print >> sys.stderr, ("Insufficient data.\n"),
                return None
            start1 = hic.fends['fends']['mid'][startfend1]
            if binsize > 0:
                start1 = (start1 / binsize) * binsize
        elif start1 is None:
            start1 = hic.fends['fends']['mid'][startfend1]
            if binsize > 0:
                start1 = (start1 / binsize) * binsize
        else:
            startfend1 = _find_fend_from_coord(hic, chrint1, start1)
        if (stop1 is None or stop1 == 0) and stopfend1 is None:
            stopfend1 = hic.fends['chr_indices'][chrint1 + 1]
            while stopfend1 > hic.fends['chr_indices'][chrint1] and hic.filter[stopfend1 - 1] == 0:
                stopfend1 -= 1
            stop1 = hic.fends['fends']['mid'][stopfend1 - 1]
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        elif stop1 is None or stop1 == 0:
            stop1 = hic.fends['fends']['mid'][stopfend1 - 1]
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        else:
            if binsize > 0:
                stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
            stopfend1 = _find_fend_from_coord(hic, chrint1, stop1)
    if not binbounds1 is None:
        start2 = binbounds1[0, 0]
        stop2 = binbounds1[-1, 1]
        startfend2 = _find_fend_from_coord(hic, chrint2, start2)
        stopfend2 = _find_fend_from_coord(hic, chrint2, stop2)
    else:
        if start2 is None and startfend2 is None:
            startfend2 = hic.fends['chr_indices'][chrint2]
            while startfend2 < hic.fends['chr_indices'][chrint2 + 1] and hic.filter[startfend2] == 0:
                startfend2 += 1
            if startfend2 == hic.fends['chr_indices'][chrint2 + 1]:
                if not silent:
                    print >> sys.stderr, ("Insufficient data.\n"),
                return None
            start2 = hic.fends['fends']['mid'][startfend2]
            if binsize > 0:
                start2 = (start2 / binsize) * binsize
        elif start2 is None:
            start2 = hic.fends['fends']['mid'][startfend2]
            if binsize > 0:
                start2 = (start2 / binsize) * binsize
        else:
            startfend2 = _find_fend_from_coord(hic, chrint2, start2)
        if (stop2 is None or stop2 == 0) and stopfend2 is None:
            stopfend2 = hic.fends['chr_indices'][chrint2 + 1]
            while stopfend2 > hic.fends['chr_indices'][chrint2] and hic.filter[stopfend2 - 1] == 0:
                stopfend2 -= 1
            stop2 = hic.fends['fends']['mid'][stopfend2 - 1]
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        elif stop2 is None or stop2 == 0:
            stop2 = hic.fends['fends']['mid'][stopfend2 - 1]
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        else:
            if binsize > 0:
                stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
            stopfend2 = _find_fend_from_coord(hic, chrint2, stop2)
    if not silent:
        print >> sys.stderr, ("Finding %s array for %s:%i-%i by %s:%i-%i...") % (datatype,  chrom1,
                                                                                 start1, stop1, chrom2, start2,
                                                                                 stop2),
    # If datatype is not 'expected', pull the needed slice of data
    if datatype != 'expected':
        if chrint1 < chrint2:
            start_index = hic.data['trans_indices'][startfend1]
            stop_index = hic.data['trans_indices'][stopfend1]
        else:
            start_index = hic.data['trans_indices'][startfend2]
            stop_index = hic.data['trans_indices'][stopfend2]
        if start_index == stop_index:
            if not silent:
                print >> sys.stderr, ("Insufficient data\n"),
            return None
        if chrint1 < chrint2:
            data_indices = hic.data['trans_indices'][startfend1:(stopfend1 + 1)]
        else:
            data_indices = hic.data['trans_indices'][startfend2:(stopfend2 + 1)]
        data_indices -= data_indices[0]
        data = hic.data['trans_data'][start_index:stop_index, :]
        if chrint1 < chrint2:
            data[:, 0] -= startfend1
            data[:, 1] -= startfend2
        else:
            data[:, 0] -= startfend2
            data[:, 1] -= startfend1
    else:
        data_indices = None
        data = None
    # Determine mapping of valid fends to bins
    mapping1 = numpy.zeros(stopfend1 - startfend1, dtype=numpy.int32) - 1
    mapping2 = numpy.zeros(stopfend2 - startfend2, dtype=numpy.int32) - 1
    valid1 = numpy.where(hic.filter[startfend1:stopfend1] > 0)[0]
    valid2 = numpy.where(hic.filter[startfend2:stopfend2] > 0)[0]
    mids1 = hic.fends['fends']['mid'][startfend1:stopfend1]
    mids2 = hic.fends['fends']['mid'][startfend2:stopfend2]
    if binsize == 0 and binbounds1 is None:
        if skipfiltered:
            mapping1[valid1] = numpy.arange(valid1.shape[0])
            num_bins1 = valid1.shape[0]
        else:
            mapping1[valid1] = valid1
            num_bins1 = mapping1.shape[0]
    elif not binbounds1 is None:
        start_indices = numpy.searchsorted(binbounds1[:, 0], mids1[valid1], side='right') - 1
        stop_indices = numpy.searchsorted(binbounds1[:, 1], mids1[valid1], side='right')
        where = numpy.where(start_indices == stop_indices)[0]
        mapping1[valid1[where]] = start_indices[where]
        num_bins1 = binbounds1.shape[0]
    else:
        mapping1[valid1] = (mids1[valid1] - start1) / binsize
        num_bins1 = (stop1 - start1) / binsize
    if binsize == 0 and binbounds2 is None:
        if skipfiltered:
            mapping2[valid2] = numpy.arange(valid2.shape[0])
            num_bins2 = valid2.shape[0]
        else:
            mapping2[valid2] = valid2
            num_bins2 = mapping2.shape[0]
    elif not binbounds2 is None:
        start_indices = numpy.searchsorted(binbounds2[:, 0], mids2[valid2], side='right') - 1
        stop_indices = numpy.searchsorted(binbounds2[:, 1], mids2[valid2], side='right')
        where = numpy.where(start_indices == stop_indices)[0]
        mapping2[valid2[where]] = start_indices[where]
        num_bins2 = binbounds2.shape[0]
    else:
        mapping2[valid2] = (mids2[valid2] - start2) / binsize
        num_bins2 = (stop2 - start2) / binsize
    # Find maximum interaction partner for each fend
    if num_bins1 < 1 or num_bins2 < 1:
        if not silent:
            print >> sys.stderr, ("Insufficient data\n"),
        return None
    # If correction is required, determine what type and get appropriate data
    corrections1 = None
    corrections2 = None
    binning_corrections = None
    correction_indices = None
    binning_num_bins = None
    fend_indices = None
    if datatype in ['fend', 'enrichment', 'expected']:
        if hic.normalization in ['express', 'probability', 'binning-express', 'binning-probability']:
            corrections1 = hic.corrections[startfend1:stopfend1]
            corrections2 = hic.corrections[startfend2:stopfend2]
        if hic.normalization in ['binning', 'binning-express', 'binning-probability']:
            binning_corrections = hic.binning_corrections
            correction_indices = hic.binning_correction_indices
            binning_num_bins = hic.binning_num_bins
            fend_indices = hic.binning_fend_indices
    if datatype in ['distance', 'enrichment', 'expected']:
        if 'trans_means' not in hic.__dict__.keys():
            hic.find_trans_means()
        if chrint1 < chrint2:
            index = chrint1 * (hic.fends['chromosomes'].shape[0] - 1) - chrint1 * (chrint1 + 1) / 2 - 1 + chrint2
        else:
            index = chrint2 * (hic.fends['chromosomes'].shape[0] - 1) - chrint2 * (chrint2 + 1) / 2 - 1 + chrint1
        trans_mean = hic.trans_means[index]
    else:
        trans_mean = 1.0
    # Create data array
    if chrint1 < chrint2:
        data_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins2, num_bins1, 2), dtype=numpy.float32)
    # Fill in data values
    if chrint1 < chrint2:
        if datatype != 'raw':
            _hic_binning.find_trans_expected(mapping1, mapping2, corrections1, corrections2, binning_corrections,
                                             correction_indices, binning_num_bins, fend_indices, data_array,
                                             trans_mean, startfend1, startfend2)
        if datatype != 'expected':
            _hic_binning.find_trans_observed(data, data_indices, mapping1, mapping2, data_array)
    else:
        if datatype != 'raw':
            _hic_binning.find_trans_expected(mapping2, mapping1, corrections2, corrections1, binning_corrections,
                                             correction_indices, binning_num_bins, fend_indices, data_array,
                                             trans_mean, startfend2, startfend1)
        if datatype != 'expected':
            _hic_binning.find_trans_observed(data, data_indices, mapping2, mapping1, data_array)
    if chrint2 < chrint1:
        data_array = numpy.transpose(data_array, (1, 0, 2))
    if datatype == 'expected':
        where = numpy.where(data_array[:, :, 1] > 0.0)
        data_array[where[0], where[1], 0] = 1.0
    if datatype == 'raw':
        where = numpy.where(data_array[:, :, 0] > 0.0)
        data_array[where[0], where[1], 1] = 1.0
    if returnmapping:
        bin_mapping1 = numpy.zeros((num_bins1, 4), dtype=numpy.int32)
        if binsize == 0 and binbounds1 is None:
            if skipfiltered:
                bin_mapping1[:, 2] = valid1 + startfend1
            else:
                bin_mapping1[:, 2] = numpy.arange(startfend1, stopfend1)
            bin_mapping1[:, 3] = bin_mapping1[:, 2] + 1
            bin_mapping1[:, 0] = hic.fends['fends']['start'][bin_mapping1[:, 2]]
            bin_mapping1[:, 1] = hic.fends['fends']['stop'][bin_mapping1[:, 2]]
        else:
            if binbounds1 is None:
                bin_mapping1[:, 0] = start1 + binsize * numpy.arange(num_bins1)
                bin_mapping1[:, 1] = bin_mapping1[:, 0] + binsize
            else:
                bin_mapping1[:, :2] = binbounds1
            bin_mapping1[:, 2] = numpy.searchsorted(mids1, bin_mapping1[:, 0]) + startfend1
            bin_mapping1[:, 3] = numpy.searchsorted(mids1, bin_mapping1[:, 1]) + startfend1
        bin_mapping2 = numpy.zeros((num_bins2, 4), dtype=numpy.int32)
        if binsize == 0 and binbounds2 is None:
            if skipfiltered:
                bin_mapping2[:, 2] = valid2 + startfend2
            else:
                bin_mapping2[:, 2] = numpy.arange(startfend2, stopfend2)
            bin_mapping2[:, 3] = bin_mapping2[:, 2] + 1
            bin_mapping2[:, 0] = hic.fends['fends']['start'][bin_mapping2[:, 2]]
            bin_mapping2[:, 1] = hic.fends['fends']['stop'][bin_mapping2[:, 2]]
        else:
            if binbounds2 is None:
                bin_mapping2[:, 0] = start2 + binsize * numpy.arange(num_bins2)
                bin_mapping2[:, 1] = bin_mapping2[:, 0] + binsize
            else:
                bin_mapping2[:, :2] = binbounds2
            bin_mapping2[:, 2] = numpy.searchsorted(mids2, bin_mapping2[:, 0]) + startfend2
            bin_mapping2[:, 3] = numpy.searchsorted(mids2, bin_mapping2[:, 1]) + startfend2
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [data_array, bin_mapping1, bin_mapping2]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return data_array

def bin_trans_array(data_array, data_mapping1, data_mapping2, binsize=10000, binbounds1=None, start1=None, stop1=None,
                    binbounds2=None, start2=None, stop2=None, returnmapping=False, **kwargs):
    """
    Create an array of format 'arraytype' and fill 'binsize' bins or bins defined by 'binbounds' with data provided in the array passed by 'unbinned'.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
    :param data_array: A 3d array containing data to be binned.
    :type data_array: numpy array
    :param data_mapping1: An N x 4 2d integer array containing the start and stop coordinates, and start and stop fends for each of the N bin ranges along the first axis in 'data_array'.
    :type data_mapping1: numpy array
    :param data_mapping2: An N x 4 2d integer array containing the start and stop coordinates, and start and stop fends for each of the N bin ranges along the second axis in 'data_array'.
    :type data_mapping2: numpy array
    :param binsize: This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
    :type binsize: int.
    :param binbounds1: An array containing start and stop coordinates for a set of user-defined bins along the first axis. Any bin from 'data_array' not falling in a bin is ignored.
    :type binbounds1: numpy array
    :param start1: The coordinate at the beginning of the first bin for the first axis of the binned data. If unspecified, 'start1' will be the first multiple of 'binsize' below the first coordinate from 'data_mapping1'. If 'binbounds1' is given, 'start1' is ignored. Optional.
    :type start1: int.
    :param stop1: The coordinate at the end of the last bin for the first axis of the binned data. If unspecified, 'stop1' will be the first multiple of 'binsize' after the last coordinate from 'data_mapping1'. If needed, 'stop1' is adjusted upward to create a complete last bin. If 'binbounds1' is given, 'stop1' is ignored. Optional.
    :type stop1: int.
    :param binbounds2: An array containing start and stop coordinates for a set of user-defined bins along the second axis. Any bin from 'data_array' not falling in a bin is ignored.
    :type binbounds2: numpy array
    :param start2: The coordinate at the beginning of the first bin for the second axis of the binned data. If unspecified, 'start2' will be the first multiple of 'binsize' below the first coordinate from 'data_mapping2'. If 'binbounds2' is given, 'start2' is ignored. Optional.
    :type start2: int.
    :param stop2: The coordinate at the end of the last bin for the second axis of the binned data. If unspecified, 'stop2' will be the first multiple of 'binsize' after the last coordinate from 'data_mapping2'. If needed, 'stop2' is adjusted upward to create a complete last bin. If 'binbounds2' is given, 'stop2' is ignored. Optional.
    :type stop2: int.
    :param datatype: This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend', 'enrichment', and 'expected'. Observed values are always in the first index along the last axis, except when 'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, unfiltered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both correction and distance mean values.
    :type datatype: str.
    :param returnmapping: If 'True', a list containing the data array and a 2d array containing first coordinate included and excluded from each bin, and the first fend included and excluded from each bin is returned. Otherwise only the data array is returned.
    :type returnmapping: bool.
    :returns: Array in format requested with 'arraytype' containing binned data requested with 'datatype' pulled from 'unbinned'.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    # Determine start and stop, if necessary
    if binbounds1 is None:
        if start1 is None:
            start1 = (data_mapping1[0, 0] / binsize) * binsize
        if stop1 is None:
            stop1 = ((data_mapping1[-1, 1] - 1) / binsize + 1) * binsize
        else:
            stop1 = ((stop1 - 1 - start1) / binsize + 1) * binsize + start1
        num_bins1 = (stop1 - start1) / binsize
        binbounds1 = numpy.zeros((num_bins1, 2), dtype=numpy.int32)
        binbounds1[:, 0] = numpy.arange(num_bins1) * binsize + start1
        binbounds1[:, 1] = binbounds1[:, 0] + binsize
    else:
        num_bins1 = binbounds1.shape[0]
        start1 = binbounds1[0, 0]
        stop1 = binbounds1[0, 1]
    if binbounds2 is None:
        if start2 is None:
            start2 = (data_mapping2[0, 0] / binsize) * binsize
        if stop2 is None:
            stop2 = ((data_mapping2[-1, 1] - 1) / binsize + 1) * binsize
        else:
            stop2 = ((stop2 - 1 - start2) / binsize + 1) * binsize + start2
        num_bins2 = (stop2 - start2) / binsize
        binbounds2 = numpy.zeros((num_bins2, 2), dtype=numpy.int32)
        binbounds2[:, 0] = numpy.arange(num_bins2) * binsize + start2
        binbounds2[:, 1] = binbounds2[:, 0] + binsize
    else:
        num_bins2 = binbounds2.shape[0]
        start2 = binbounds2[0, 0]
        stop2 = binbounds2[0, 1]
    mids1 = (data_mapping1[:, 0] + data_mapping1[:, 1]) / 2
    mids2 = (data_mapping2[:, 0] + data_mapping2[:, 1]) / 2
    if not silent:
        print >> sys.stderr, ("Finding binned trans array..."),
    # Find bin mapping for each fend
    mapping1 = numpy.zeros(mids1.shape[0], dtype=numpy.int32) - 1
    fend_ranges1 = numpy.zeros((binbounds1.shape[0], 2), dtype=numpy.int32)
    for i in range(binbounds1.shape[0]):
        firstbin = numpy.searchsorted(mids1, binbounds1[i, 0])
        lastbin = numpy.searchsorted(mids1, binbounds1[i, 1])
        mapping1[firstbin:lastbin] = i
        fend_ranges1[i, 0] = data_mapping1[firstbin, 2]
        fend_ranges1[i, 1] = data_mapping1[lastbin, 3]
    valid1 = numpy.where(mapping1 >= 0)[0]
    mapping2 = numpy.zeros(mids2.shape[0], dtype=numpy.int32) - 1
    fend_ranges2 = numpy.zeros((binbounds2.shape[0], 2), dtype=numpy.int32)
    for i in range(binbounds2.shape[0]):
        firstbin = numpy.searchsorted(mids2, binbounds2[i, 0])
        lastbin = numpy.searchsorted(mids2, binbounds2[i, 1])
        mapping2[firstbin:lastbin] = i
        fend_ranges2[i, 0] = data_mapping2[firstbin, 2]
        fend_ranges2[i, 1] = data_mapping2[lastbin, 3]
    valid2 = numpy.where(mapping2 >= 0)[0]
    # Create requested array
    binned_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
    # Fill in binned data values
    for i in range(valid1.shape[0]):
        binned_array[i, :, 0] = numpy.bincount(mapping2[valid2], weights=data_array[valid1[i], valid2, 0],
                                            minlength=num_bins2)
        binned_array[i, :, 1] = numpy.bincount(mapping2[valid2], weights=data_array[valid1[i], valid2, 1],
                                            minlength=num_bins2)
    # If mapping requested, calculate bin bounds
    if returnmapping:
        mapping1 = numpy.zeros((num_bins1, 4), dtype=numpy.int32)
        mapping1[:, 0] = binbounds1[:, 0]
        mapping1[:, 1] = binbounds1[:, 1]
        mapping1[:, 2:4] = fend_ranges1
        mapping2 = numpy.zeros((num_bins2, 4), dtype=numpy.int32)
        mapping2[:, 0] = binbounds2[:, 0]
        mapping2[:, 1] = binbounds2[:, 1]
        mapping2[:, 2:4] = fend_ranges2
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return [binned_array, mapping1, mapping2]
    else:
        if not silent:
            print >> sys.stderr, ("Done\n"),
        return binned_array

def dynamically_bin_trans_array(unbinned, unbinnedpositions1, unbinnedpositions2, binned, binbounds1, binbounds2,
                                minobservations=10, searchdistance=0, removefailed=False, **kwargs):
    """
    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations', or 'searchdistance' criteria.

    :param unbinned: A 3d array containing data to be used for filling expanding bins. This array should be  N x M x 2, where N is the number of bins or fends from the first chromosome and M is the number of bins or fends from the second chromosome.
    :type unbinned: numpy array
    :param unbinnedpositions1: A 2d integer array indicating the first and last coordinate of each bin along the first axis in 'unbinned' array.
    :type unbinnedpositions1: numpy array
    :param unbinnedpositions2: A 2d integer array indicating the first and last coordinate of each bin along the first axis in 'unbinned' array.
    :type unbinnedpositions2: numpy array
    :param binned: A 3d array containing binned data to be dynamically binned. This array should be  N x M x 2, where N is the number of bins from the first chromosome and M is the number of bins from the second chromosome. Data in this array will be altered by this function.
    :type binned: numpy array
    :param binbounds1: An integer array indicating the start and end position of each bin from the first chromosome in the 'binned' array. This array should be N x 2, where N is the size of the first dimension of 'binned'.
    :type binbounds1: numpy array
    :param binbounds2: An integer array indicating the start and end position of each bin from the second chromosome in the 'binned' array. This array should be N x 2, where N is the size of the second dimension of 'binned'.
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
    _hic_binning.dynamically_bin_trans(unbinned, unbinnedmids1, unbinnedmids2, binned, binedges1,
                                   binedges2, mids1, mids2, minobservations, searchdistance, int(removefailed))
    if not silent:
        print >> sys.stderr, ("Done\n"),
    return None

def write_heatmap_dict(hic, filename, binsize, includetrans=True, datatype='enrichment', chroms=[], 
                       dynamically_binned=False, minobservations=0, searchdistance=0, expansion_binsize=0,
                       removefailed=False, **kwargs):
    """
    Create an h5dict file containing binned interaction arrays, bin positions, and an index of included chromosomes. This function is MPI compatible.

    :param hic: A :class:`HiC <hifive.hic.HiC>` class object containing fend and count data.
    :type hic: :class:`HiC <hifive.hic.HiC>`
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
    :returns: None
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
    if ('silent' in kwargs and kwargs['silent']) or rank > 0:
        silent = True
    else:
        silent = False
    # Check if trans mean is needed and calculate if not already done
    if includetrans and datatype in ['distance', 'enrichment'] and 'trans_mean' not in hic.__dict__.keys():
        hic.find_trans_means()
    # Check if filename already exists, and remove if it does
    if rank == 0:
        if os.path.exists(filename):
            if not silent:
                print >> sys.stderr, ("%s already exists, overwriting.") % filename
            subprocess.call('rm %s' % filename, shell=True)
        if not silent:
            print >> sys.stderr, ("Creating binned heatmap...\n"),
        output = h5py.File(filename, 'w')
        output.attrs['resolution'] = binsize
        # If chromosomes not specified, fill list
        if len(chroms) == 0:
            chroms = list(hic.fends['chromosomes'][...])
        # Assemble list of requested arrays
        needed = []
        chr_indices = hic.fends['chr_indices'][...]
        for i in range(len(chroms))[::-1]:
            chrom = chroms[i]
            chrint = hic.chr2int[chrom]
            if numpy.sum(hic.filter[chr_indices[chrint]:chr_indices[chrint + 1]]) > 0:
                needed.append((chrom,))
            else:
                del chroms[i]
        if includetrans:
            for i in range(len(chroms)-1):
                for j in range(i + 1, len(chroms)):
                    needed.append((chroms[i],chroms[j]))
        if num_procs == 1:
            node_needed = needed
        else:
            node_ranges = numpy.round(numpy.linspace(0, len(needed), num_procs + 1)).astype(numpy.int32)
            for i in range(1, num_procs):
                comm.send(needed[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
            node_needed = needed[node_ranges[0]:node_ranges[1]]
    else:
        node_needed = comm.recv(source=0, tag=11)
    heatmaps = {}
    # Find heatmaps
    for chrom in node_needed:
        if len(chrom) == 1:
            # Find cis heatmap
            # determine if data is to be dynamically binned
            if not dynamically_binned:
                heatmaps[chrom] = find_cis_signal(hic, chrom[0], binsize=binsize, datatype=datatype,
                                                  arraytype='upper', returnmapping=True, silent=silent,
                                                  skipfiltered=True)
            else:
                temp = find_cis_signal(hic, chrom[0], binsize=expansion_binsize, datatype=datatype, arraytype='upper',
                                       returnmapping=True, silent=silent)
                if temp is None:
                    continue
                expansion, exp_mapping = temp
                binned, mapping = find_cis_signal(hic, chrom[0], binsize=binsize, datatype=datatype,
                                                  arraytype='upper', returnmapping=True, silent=silent)
                dynamically_bin_cis_array(expansion, exp_mapping, binned, mapping, minobservations=minobservations,
                                          searchdistance=searchdistance, removefailed=removefailed, silent=silent)

                heatmaps[chrom] = [binned, mapping]
        else:
            # Find trans heatmap
            # determine if data is to be dynamically binned
            if not dynamically_binned:
                heatmaps[chrom] = find_trans_signal(hic, chrom[0], chrom[1],  binsize=binsize, datatype=datatype,
                                                    returnmapping=False, silent=silent, skipfiltered=True)
            else:
                temp = find_trans_signal(hic, chrom[0], chrom[1], binsize=expansion_binsize, datatype=datatype, 
                                         returnmapping=True, silent=silent)
                if temp is None:
                    continue
                expansion, exp_mapping1, exp_mapping2 = temp
                binned, mapping1, mapping2 = find_trans_signal(hic, chrom[0], chrom[1], binsize=binsize,
                                                               datatype=datatype, returnmapping=True, silent=silent)
                dynamically_bin_trans_array(expansion, exp_mapping1, exp_mapping2, binned, mapping1, mapping2,
                                            minobservations=minobservations, searchdistance=searchdistance,
                                            removefailed=removefailed, silent=silent)

                heatmaps[chrom] = binned
        # Check if array contains data
        if heatmaps[chrom] is None or heatmaps[chrom][0].shape[0] == 0:
            del heatmaps[chrom]
    # Collect heatmaps at node 0 and write to h5dict
    if rank == 0:
        if num_procs > 1:
            for i in range(1, num_procs):
                if node_ranges[i + 1] - node_ranges[i] > 0:
                    temp = comm.recv(source=i, tag=11)
                    heatmaps.update(temp)
            del temp
        for chrom in heatmaps.keys():
            if len(chrom) == 1:
                output.create_dataset('%s.counts' % chrom[0], data=heatmaps[chrom][0][:, 0])
                output.create_dataset('%s.expected' % chrom[0], data=heatmaps[chrom][0][:, 1])
                output.create_dataset('%s.positions' % chrom[0], data=heatmaps[chrom][1][:, :2])
            else:
                output.create_dataset('%s_by_%s.counts' % (chrom[0], chrom[1]), data=heatmaps[chrom][:, :, 0])
                output.create_dataset('%s_by_%s.expected' % (chrom[0], chrom[1]), data=heatmaps[chrom][:, :, 1])
        output.create_dataset('chromosomes', data=numpy.array(chroms))
        if 'history' in kwargs:
            output.attrs['history'] = kwargs['history']
        output.close()
        if not silent:
            print >> sys.stderr, ("Creating binned heatmap...Done\n"),
    else:
        if len(heatmaps) > 0:
            comm.send(heatmaps, dest=0, tag=11)
        del heatmaps
    return None
