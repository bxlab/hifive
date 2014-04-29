#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

"""
This is a module contains scripts for generating compact and full matrices of
interaction data.

Input data
----------

These functions rely on the "HiC" class in conjunction with the "Fend" and
"HiCData" classes.

Concepts
--------

Data can either be arranged in compact, complete, or flattened (row-major)
upper-triangle arrays. Compact arrays are N x M, where N is the number of
fends or bins, and M is the maximum distance between fends or bins. This is
useful for working with sets of short interactions. Data can be raw,
fend-corrected, distance-dependence removed, or enrichment values. Arrays
are 3-dimensional with observed values in the first layer of d3, expected
values in the second layer of d3. The exception to this is upper-triangle
arrays, which are 2d, divinding observed and expected along the second axis.

-----------------------------------------------------------------------------

API documentation
-----------------



"""

import os
import sys
import subprocess
from math import floor, ceil

import numpy
import h5py
try:
    from mpi4py import MPI
except:
    pass

import _distance
import _binning


def unbinned_cis_signal(hic, chrom, start=None, stop=None, startfend=None, stopfend=None, datatype='enrichment',
                        arraytype='compact', maxdistance=0, skipfiltered=False, returnmapping=False):
    """
    unbinned_cis_signal method

    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    Parameters
    ----------
    hic : HiC class object
        A HiC class object containing fend and count data.
    chrom : str
        The name of a chromosome contained in 'hic'.
    start : int, optional
        The smallest coordinate to include in the array, measured from fend midpoints. If both 'start' and 'startfend'
        are given, 'start' will override 'startfend'. If unspecified, this will be set to the midpoint of the first
        fend for 'chrom'.
    stop : int, optional
        The largest coordinate to include in the array, measured from fend midpoints. If both 'stop' and 'stopfend'
        are given, 'stop' will override 'stopfend'. If unspecified, this will be set to the midpoint of the last
        fend plus one for 'chrom'.
    startfend : int, optional
        The first fend to include in the array. If unspecified and 'start' is not given, this is set to the first fend
        in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given preference.
    stopfend : int, optional
        The first fend not to include in the array. If unspecified and 'stop' is not given, this is set to the last
        fend in 'chrom' plus one. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given
        preference.
    datatype : str, optional
        This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend',
        'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when
        'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-
        filtered fends return value of one. Expected values are returned for 'distance', 'fend', 'enrichment', and
        'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the
        expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both
        correction and distance mean values.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'.
        'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number
        of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values
        between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the
        flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    maxdistance : int, optional
        This specifies the maximum coordinate distance between bins that will be included in the array. If set to
        zero, all distances are included.
    skipfiltered : bool, optional
        If 'True', all interaction bins for filtered out fends are removed and a reduced-size array is returned.
    returnmapping : bool, optional
        If 'True', a list containing the data array and a 1d array containing fend numbers included in the data array
        is return. Otherwise only the data array is returned.
    """
    # check that all values are acceptable
    datatypes = {'raw': 0, 'fend': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
    if datatype not in datatypes:
        print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    else:
        datatype_int = datatypes[datatype]
    if arraytype not in ['full', 'compact', 'upper']:
        print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine start, stop, startfend, and stopfend
    chrint = hic.chr2int[chrom.strip('chr')]
    if start is None and startfend is None:
        startfend = hic.fends['chr_indices'][chrint]
        start = hic.fends['fends']['mid'][startfend]
    elif start is None:
        start = hic.fends['fends']['mid'][startfend]
    else:
        startfend = _find_fend_from_coord(hic, chrint, start)
    if stop is None and stopfend is None:
        stopfend = hic.fends['chr_indices'][chrint + 1]
        stop = hic.fends['fends']['mid'][stopfend - 1] + 1
    elif stop is None:
        stop = hic.fends['fends']['mid'][stopfend - 1] + 1
    else:
        stopfend = _find_fend_from_coord(hic, chrint, stop)
    print >> sys.stderr, ("Finding %s %s array for %s:%i-%i...") % (datatype, arraytype, chrom, start, stop),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        start_index = hic.data['cis_indices'][startfend]
        stop_index = hic.data['cis_indices'][stopfend]
        if start_index == stop_index:
            print >> sys.stderr, ("Insufficient data\n"),
            return None
        data_indices = hic.data['cis_indices'][startfend:(stopfend + 1)]
        data_indices -= data_indices[0]
        data = hic.data['cis_data'][start_index:stop_index, :]
    else:
        data_indices = None
        data = None
    # If skipfiltered is True, map only valid fends
    if skipfiltered:
        mapping = numpy.where(hic.filter[startfend:stopfend] > 0)[0].astype(numpy.int32) + startfend
    else:
        mapping = numpy.arange(stopfend - startfend, dtype=numpy.int32) + startfend
    # Find maximum interaction partner for each fend
    num_bins = mapping.shape[0]
    if num_bins < 2:
        print >> sys.stderr, ("Insufficient data\n"),
        return None
    max_fend = numpy.zeros(num_bins, dtype=numpy.int32)
    mids = hic.fends['fends']['mid'][mapping]
    _distance.find_max_fend(max_fend, mids, hic.fends['fends']['chr'][mapping],
                            hic.fends['chr_indices'][...], startfend, maxdistance)
    max_fend = numpy.minimum(max_fend, num_bins)
    max_bin = numpy.amax(max_fend - numpy.arange(num_bins))
    if max_bin <= 0:
        print >> sys.stderr, ("Insufficient data.\n"),
        return None
    # Create requested array
    if arraytype == 'compact':
        data_array = numpy.zeros((num_bins, max_bin, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in data values
    if arraytype == 'compact':
        _binning.unbinned_signal_compact(data, data_indices, hic.filter, mapping,
                                         hic.corrections, mids, hic.distance_mids,
                                         hic.distance_means, max_fend, data_array, datatype_int, startfend)
    else:
        _binning.unbinned_signal_upper(data, data_indices, hic.filter, mapping,
                                       hic.corrections, mids, hic.distance_mids,
                                       hic.distance_means, max_fend, data_array, datatype_int, startfend)
    # If requesting 'full' array, convert 'upper' array type to 'full'
    if arraytype == 'full':
        indices = numpy.triu_indices(num_bins, 1)
        full_data_array = numpy.zeros((num_bins, num_bins, 2), dtype=numpy.float32)
        full_data_array[indices[1], indices[0], :] = data_array
        full_data_array[indices[0], indices[1], :] = data_array
        del data_array
        data_array = full_data_array
    if returnmapping:
        print >> sys.stderr, ("Done\n"),
        return [data_array, mapping]
    else:
        print >> sys.stderr, ("Done\n"),
        return data_array


def _find_fend_from_coord(hic, chrint, coord):
    """Find the next fend after the coordinate on chromosome 'chrint'."""
    first_fend = hic.fends['chr_indices'][chrint]
    last_fend = hic.fends['chr_indices'][chrint + 1]
    return numpy.searchsorted(hic.fends['fends']['mid'][first_fend:last_fend], coord) + first_fend


def bin_cis_signal(hic, chrom, start=None, stop=None, startfend=None, stopfend=None, binsize=10000,
                   binbounds=None, datatype='enrichment', arraytype='full', maxdistance=0, returnmapping=False):
    """
    bin_cis_signal method

    Create an array of format 'arraytype' and fill 'binsize' bins or user-defined 'binbound' bins with data requested
    in 'datatype'.

    Parameters
    ----------
    hic : HiC class object
        A HiC class object containing fend and count data.
    chrom : str
        The name of a chromosome contained in 'hic'.
    start : int, optional
        The coordinate at the beginning of the smallest bin. If unspecified, 'start' will be the first multiple of
        'binsize' below the 'startfend' mid. If there is a conflict between 'start' and 'startfend', 'start' is given
        preference.
    stop : int, optional
        The largest coordinate to include in the array, measured from fend midpoints. If both 'stop' and 'stopfend'
        are given, 'stop' will override 'stopfend'.
    startfend : int, optional
        The first fend to include in the array. If unspecified and 'start' is not given, this is set to the first
        valid fend in 'chrom'. In cases where 'start' is specified and conflicts with 'startfend', 'start' is given
        preference.
    stopfend : int, optional
        The first fend not to include in the array. If unspecified and 'stop' is not given, this is set to the last
        valid fend in 'chrom' + 1. In cases where 'stop' is specified and conflicts with 'stopfend', 'stop' is given
        preference.
    binsize : int, optional
        This is the coordinate width of each bin. If binbounds is not None, this value is ignored.
    binbounds : 2D numpy array
        An array containing start and stop coordinates for a set of user-defined bins. Any fends not falling in a bin
        is ignored. 
    datatype : str, optional
        This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend',
        'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when
        'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-
        filtered bins return value of 1. Expected values are returned for 'distance', 'fend', 'enrichment', and
        'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the
        expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both
        correction and distance mean values.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'.
        'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number
        of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values
        between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only
        theflattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    maxdistance : int, optional
        This specifies the maximum coordinate distance between bins that will be included in the array. If set to
        zero, all distances are included.
    returnmapping : bool, optional
        If 'True', a list containing the data array and a 2d array of N x 4 containing the first fend and last fend
        plus one included in each bin and first and last coordinates is return. Otherwise only the data array is
        returned.
    """
    # check that all values are acceptable
    datatypes = {'raw': 0, 'fend': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
    if datatype not in datatypes:
        print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    else:
        datatype_int = datatypes[datatype]
    if arraytype not in ['full', 'compact', 'upper']:
        print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine start, stop, startfend, and stopfend
    chrint = hic.chr2int[chrom.strip('chr')]
    if start is None and startfend is None:
        startfend = hic.fends['chr_indices'][chrint]
        while startfend < hic.fends['chr_indices'][chrint + 1] and hic.filter[startfend] == 0:
            startfend += 1
        if startfend == hic.fends['chr_indices'][chrint + 1]:
            print >> sys.stderr, ("Insufficient data.\n"),
            return None
        start = (hic.fends['fends']['mid'][startfend] / binsize) * binsize
    elif start is None:
        start = (hic.fends['fends']['mid'][startfend] / binsize) * binsize
    else:
        startfend = _find_fend_from_coord(hic, chrint, start)
    if stop is None and stopfend is None:
        stopfend = hic.fends['chr_indices'][chrint + 1]
        while stopfend > hic.fends['chr_indices'][chrint] and hic.filter[stopfend - 1] == 0:
            stopfend -= 1
        stop = int(ceil((hic.fends['fends']['mid'][stopfend - 1] + 1 - start) / float(binsize))) * binsize + start
    elif stop is None:
        stop = int(ceil((hic.fends['fends']['mid'][stopfend - 1] + 1 - start) / float(binsize))) * binsize + start
    else:
        stopfend = _find_fend_from_coord(hic, chrint, stop)
    if binbounds is None:
        num_bins = (stop - start) / binsize
    else:
        num_bins = binbounds.shape[0]
    num_fends = stopfend - startfend
    print >> sys.stderr, ("Finding %s %s array for %s:%i-%i...") % (datatype, arraytype, chrom, start, stop),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        start_index = hic.data['cis_indices'][startfend]
        stop_index = hic.data['cis_indices'][stopfend]
        data_indices = hic.data['cis_indices'][startfend:(stopfend + 1)]
        data_indices -= data_indices[0]
        data = hic.data['cis_data'][start_index:stop_index, :]
        data[:, :2] -= startfend
    else:
        data_indices = None
        data = None
    # Find fend ranges for each bin
    mids = hic.fends['fends']['mid'][startfend:stopfend]
    if binbounds is None:
        mapping = numpy.searchsorted(numpy.arange(1, num_bins + 1) * binsize + start, mids).astype(numpy.int32)
    else:
        mapping = numpy.zeros(mids.shape[0], dtype=numpy.int32) - 1
        starts = numpy.searchsorted(binbounds[:, 0], mids, side='right') - 1
        stops = numpy.searchsorted(binbounds[:, 1], mids)
        where = numpy.where(starts == stops)[0]
        mapping[where] = starts[where].astype(numpy.int32)
    # Find maximum interaction partner for each fend
    max_fend = numpy.zeros(num_fends, dtype=numpy.int32)
    _distance.find_max_fend(max_fend, mids, hic.fends['fends']['chr'][startfend:stopfend],
                            hic.fends['chr_indices'][...], startfend, maxdistance)
    max_fend = numpy.minimum(max_fend, num_fends)
    if maxdistance == 0:
        max_bin = num_bins
    else:
        max_bin = int(ceil(maxdistance / float(binsize)))
    if max_bin <= 0:
        print >> sys.stderr, ("Insufficient data.\n"),
        return None
    # Create requested array
    if arraytype == 'compact':
        data_array = numpy.zeros((num_bins, max_bin, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in data values
    if arraytype == 'compact':
        _binning.binned_signal_compact(data, data_indices, hic.filter[startfend:stopfend], mapping,
                                       hic.corrections[startfend:stopfend], mids, hic.distance_mids,
                                       hic.distance_means, max_fend, data_array, datatype_int)
    else:
        _binning.binned_signal_upper(data, data_indices, hic.filter[startfend:stopfend], mapping,
                                     hic.corrections[startfend:stopfend], mids, hic.distance_mids,
                                     hic.distance_means, max_fend, data_array, datatype_int, num_bins)
    # If requesting 'full' array, convert 'upper' array type to 'full'
    if arraytype == 'full':
        indices = numpy.triu_indices(num_bins, 1)
        full_data_array = numpy.zeros((num_bins, num_bins, 2), dtype=numpy.float32)
        full_data_array[indices[1], indices[0], :] = data_array
        full_data_array[indices[0], indices[1], :] = data_array
        del data_array
        data_array = full_data_array
    # If mapping requested, calculate bin bounds
    if returnmapping:
        mapping = numpy.zeros((num_bins, 4), dtype=numpy.int32)
        mapping[:, 2] = numpy.arange(num_bins) * binsize + start
        mapping[:, 3] = mapping[:, 2] + binsize
        mapping[:, 0] = numpy.searchsorted(mids, mapping[:, 2]) + startfend
        mapping[:, 1] = numpy.searchsorted(mids, mapping[:, 3]) + startfend
        print >> sys.stderr, ("Done\n"),
        return [data_array, mapping]
    else:
        print >> sys.stderr, ("Done\n"),
        return data_array


def bin_cis_array(hic, unbinned, fends, start=None, stop=None, binsize=10000, binbounds=None, arraytype='full',
                  returnmapping=False):
    """
    bin_cis_array method

    Create an array of format 'arraytype' and fill 'binsize' bins or bins defined by 'binbounds' with data provided in
    'unbinned'.

    Parameters
    ----------
    hic : HiC class object
        A HiC class object containing fend and count data.
    unbinned : numpy array of either 'compact' or 'upper' format
        A 2d or 3d array containing data to be binned. Array format will be determined from the number of dimensions.
    fends : 1d numpy array
        An integer array indicating which position corresponds to which fend in the 'unbinned' array.
    start : int, optional
        The coordinate at the beginning of the smallest bin. If unspecified, 'start' will be the first multiple of
        'binsize' below the first mid from 'fends'. If 'binbounds' is given, 'start' is ignored.
    stop : int, optional
        The coordinate at the end of the last bin. If unspecified, 'stop' will be the first multiple of 'binsize'
        above the last mid from 'fends'. If needed, 'stop' is adjusted upward to create a complete last bin. If
        'binbounds' is given, 'stop' is ignored.
    binsize : int, optional
        This is the coordinate width of each bin. This is ignored if 'binbounds' is given.
    binbounds : numpy array
        This is an optional 2D array containing a row for each bin and start and stop coordinates.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'.
        'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number
        of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values
        between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2. 'upper' returns only the
        flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    returnmapping : bool, optional
        If 'True', a list containing the data array and a 2d array of N x 4 containing the first fend and last fend
        plus one included in each bin and first and last coordinates is return. Otherwise only the data array is
        returned.
    """
    # check that arraytype value is acceptable
    if arraytype not in ['full', 'compact', 'upper']:
        print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine input array type
    if len(unbinned.shape) == 2 and fends.shape[0] * (fends.shape[0] - 1) / 2 == unbinned.shape[0]:
        input_type = 'upper'
    elif len(unbinned.shape) == 3 and unbinned.shape[0] == fends.shape[0]:
        input_type = 'compact'
    else:
        print >> sys.stderr, ("Unrecognized input array type. No data returned.\n"),
        return None
    # Determine start and stop, if necessary
    if binbounds is None:
        if start is None:
            start = (hic.fends['fends']['mid'][fends[0]] / binsize) * binsize
        if stop is None:
            stop = int(ceil((hic.fends['fends']['mid'][fends[-1]] + 1 - start) / float(binsize))) * binsize + start
        else:
            stop = int(ceil((stop - start) / float(binsize))) * binsize + start
        num_bins = (stop - start) / binsize
        binbounds = numpy.zeros((num_bins, 2), dtype=numpy.int32)
        binbounds[:, 0] = numpy.arange(num_bins) * binsize + start
        binbounds[:, 1] = binbounds[:, 0] + binsize
    else:
        num_bins = binbounds.shape[0]
        start = binbounds[0, 0]
        stop = binbounds[0, 1]
    print >> sys.stderr, ("Finding binned %s array...") % (arraytype),
    # Find bin mapping for each fend
    mapping = numpy.zeros(fends.shape[0], dtype=numpy.int32) - 1
    mids = hic.fends['fends']['mid'][fends]
    for i in range(binbounds.shape[0]):
        firstfend = numpy.searchsorted(mids, binbounds[i, 0])
        lastfend = numpy.searchsorted(mids, binbounds[i, 1])
        mapping[firstfend:lastfend] = i
    # Create requested array
    if arraytype == 'compact':
        max_bin = (stop - start) / binsize + 1
        binned_array = numpy.zeros((num_bins, max_bin, 2), dtype=numpy.float32)
    else:
        binned_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in binned data values
    if arraytype == 'compact':
        if input_type == 'compact':
            _binning.bin_compact_to_compact(binned_array, unbinned, mapping)
        else:
            _binning.bin_upper_to_compact(binned_array, unbinned, mapping)
        # Trim unused bins
        valid = numpy.where(numpy.sum(binned_array[:, :, 1] > 0, axis=0) > 0)[0][-1]
        binned_array = binned_array[:, :(valid + 1), :]
    else:
        if input_type == 'compact':
            _binning.bin_compact_to_upper(binned_array, unbinned, mapping, num_bins)
        else:
            _binning.bin_upper_to_upper(binned_array, unbinned, mapping, num_bins)
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
        mapping[:, 2] = binbounds[:, 0]
        mapping[:, 3] = binbounds[:, 1]
        mapping[:, 0] = numpy.r_[fends, fends[-1] + 1][numpy.searchsorted(mids, mapping[:, 2])]
        mapping[:, 1] = numpy.r_[fends, fends[-1] + 1][numpy.searchsorted(mids, mapping[:, 3])]
        print >> sys.stderr, ("Done\n"),
        return [binned_array, mapping]
    else:
        print >> sys.stderr, ("Done\n"),
        return binned_array


def dynamically_bin_cis_array(unbinned, unbinnedpositions, binned, binbounds, minobservations=10,
                              searchdistance=0, removefailed=True):
    """
    dynamically_bin_cis_array method

    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations',
    or 'searchdistance' criteria.

    Parameters
    ----------
    unbinned : numpy array of either 'compact' or 'upper' format
        A 2d or 3d array containing data to be binned. Array format will be determined from the number of dimensions.
    unbinnedpositions : 1d numpy array
        An integer array indicating the mid-point of each bin in 'unbinned' array.
    binned : numpy array of either 'compact' or 'upper' format
        A 2d or 3d array containing binned data to be dynamically binned. Array format will be determined from the
        number of dimensions. Data in this array will be altered by this function.
    binbounds : 2d numpy array
        An integer array indicating the start and end position of each bin in 'binned' array.
    minobservations : int, optional
        The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
    searchdistance : int, optional
        The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on
        expansion distance.
    """
    # Determine unbinned array type
    if len(unbinned.shape) == 2 and (unbinnedpositions.shape[0] * (unbinnedpositions.shape[0] - 1) / 2 ==
                                     unbinned.shape[0]):
        unbinned_type = 'upper'
    elif len(unbinned.shape) == 3 and unbinned.shape[0] == unbinnedpositions.shape[0]:
        unbinned_type = 'compact'
    else:
        print >> sys.stderr, ("Unrecognized unbinned array type. No data returned.\n"),
        return None
    # Determine binned array type
    if len(binned.shape) == 2 and binbounds.shape[0] * (binbounds.shape[0] - 1) / 2 == binned.shape[0]:
        binned_type = 'upper'
    elif len(binned.shape) == 3 and binned.shape[0] == binbounds.shape[0]:
        binned_type = 'compact'
    else:
        print >> sys.stderr, ("Unrecognized binned array type. No data returned.\n"),
        return None
    print >> sys.stderr, ("Dynamically binning data..."),
    # Determine bin edges relative to unbinned positions
    binedges = numpy.zeros(binbounds.shape, dtype=numpy.int32)
    binedges[:, 0] = numpy.searchsorted(unbinnedpositions, binbounds[:, 0])
    binedges[:, 1] = numpy.searchsorted(unbinnedpositions, binbounds[:, 1])
    # Determine bin midpoints
    mids = (binbounds[:, 0] + binbounds[:, 1]) / 2
    # Dynamically bin using appropriate array type combination
    if unbinned_type == 'upper':
        if binned_type == 'upper':
            _binning.dynamically_bin_upper_from_upper(unbinned, unbinnedpositions, binned, binedges,
                                                      mids, minobservations, searchdistance, int(removefailed))
        else:
            _binning.dynamically_bin_compact_from_upper(unbinned, unbinnedpositions, binned, binedges,
                                                        mids, minobservations, searchdistance, int(removefailed))
    else:
        if binned_type == 'upper':
            _binning.dynamically_bin_upper_from_compact(unbinned, unbinnedpositions, binned, binedges,
                                                        mids, minobservations, searchdistance, int(removefailed))
        else:
            _binning.dynamically_bin_compact_from_compact(unbinned, unbinnedpositions, binned, binedges,
                                                          mids, minobservations, searchdistance, int(removefailed))
    print >> sys.stderr, ("Done\n"),
    return None


def dynamically_bin_trans_array(unbinned, unbinnedpositions1, unbinnedpositions2, binned, binbounds1, binbounds2,
                                minobservations=10, searchdistance=0):
    """
    dynamically_bin_trans_array method

    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations'
    or 'searchdistance' criteria.

    Parameters
    ----------
    unbinned : full numpy array
        A 3d array containing data to be binned.
    unbinnedpositions1 : 1d numpy array
        An integer array indicating the mid-point of each bin in 'unbinned' array for axis zero.
    unbinnedpositions2 : 1d numpy array
        An integer array indicating the mid-point of each bin in 'unbinned' array for axis one.
    binned : full numpy array
        A 3d array containing binned data to be dynamically binned. Data in this array will be
        altered by this function.
    binbounds1 : 2d numpy array
        An integer array indicating the start and end position of each bin in 'binned' array along axis zero.
    binbounds2 : 2d numpy array
        An integer array indicating the start and end position of each bin in 'binned' array along axis one.
    minobservations : int, optional
        The fewest number of observed reads needed for a bin to counted as valid and stop expanding.
    searchdistance : int, optional
        The furthest distance from the bin minpoint to expand bounds. If this is set to zero, there is no limit on
        expansion distance.
    """
    print >> sys.stderr, ("Dynamically binning data..."),
    # Determine bin edges relative to unbinned positions
    binedges1 = numpy.zeros(binbounds1.shape, dtype=numpy.int32)
    binedges1[:, 0] = numpy.searchsorted(unbinnedpositions1, binbounds1[:, 0])
    binedges1[:, 1] = numpy.searchsorted(unbinnedpositions1, binbounds1[:, 1])
    binedges2 = numpy.zeros(binbounds2.shape, dtype=numpy.int32)
    binedges2[:, 0] = numpy.searchsorted(unbinnedpositions2, binbounds2[:, 0])
    binedges2[:, 1] = numpy.searchsorted(unbinnedpositions2, binbounds2[:, 1])
    # Determine bin midpoints
    mids1 = (binbounds1[:, 0] + binbounds1[:, 1]) / 2
    mids2 = (binbounds2[:, 0] + binbounds2[:, 1]) / 2
    # Dynamically bin using appropriate array type combination
    _binning.dynamically_bin_trans(unbinned, unbinnedpositions1, unbinnedpositions2, binned, binedges1,
                                   binedges2, mids1, mids2, minobservations, searchdistance)
    print >> sys.stderr, ("Done\n"),
    return None


def bin_trans_signal(hic, chrom1, chrom2, start1=None, stop1=None, startfend1=None, stopfend1=None, binbounds1=None,
                     start2=None, stop2=None, startfend2=None, stopfend2=None, binbounds2=None, binsize=1000000,
                     datatype='enrichment', returnmapping=False):
    """
    bin_cis_signal method

    Create an array and fill 'binsize' bins with data requested in 'datatype'.

    Parameters
    ----------
    hic : HiC class object
        A HiC class object containing fend and count data.
    chrom1 : str
        The name of the first chromosome to pull data from.
    chrom2 : str
        The name of the second chromosome to pull data from.
    start1 : int, optional
        The coordinate at the beginning of the smallest bin from 'chrom1'. If unspecified, 'start1' will be the first
        multiple of 'binsize' below the 'startfend1' mid. If there is a conflict between 'start1' and 'startfend1',
        'start1' is given preference.
    stop1 : int, optional
        The largest coordinate to include in the array from 'chrom1', measured from fend midpoints. If both 'stop1'
        and 'stopfend1' are given, 'stop1' will override 'stopfend1'. 'stop1' will be shifted higher as needed to make
        the last bin of size 'binsize'.
    startfend1 : int, optional
        The first fend from 'chrom1' to include in the array. If unspecified and 'start1' is not given, this is set to
        the first valid fend in 'chrom1'. In cases where 'start1' is specified and conflicts with 'startfend1',
        'start1' is given preference.
    stopfend1 : int, optional
        The first fend not to include in the array from 'chrom1'. If unspecified and 'stop1' is not given, this is set
        to the last valid fend in 'chrom1' + 1. In cases where 'stop1' is specified and conflicts with 'stopfend1',
        'stop1' is given preference.
    binbounds1 : 2D numpy array
        An array containing start and stop coordinates for a set of user-defined bins for chrom1. Any fends not
        falling in a bin is ignored. If binbounds1 is specified, start1, stop1, startfend1, and stopfend1 are ignored.
    start2 : int, optional
        The coordinate at the beginning of the smallest bin from 'chrom2'. If unspecified, 'start2' will be the first
        multiple of 'binsize' below the 'startfend2' mid. If there is a conflict between 'start2' and 'startfend2',
        'start2' is given preference.
    stop2 : int, optional
        The largest coordinate to include in the array from 'chrom2', measured from fend midpoints. If both 'stop2'
        and 'stopfend2' are given, 'stop2' will override 'stopfend2'. 'stop2' will be shifted higher as needed to make
        the last bin of size 'binsize'.
    startfend2 : int, optional
        The first fend from 'chrom2' to include in the array. If unspecified and 'start2' is not given, this is set to
        the first valid fend in 'chrom2'. In cases where 'start2' is specified and conflicts with 'startfend2',
        'start2' is given preference.
    stopfend2 : int, optional
        The first fend not to include in the array from 'chrom2'. If unspecified and 'stop2' is not given, this is set
        to the last valid fend in 'chrom2' + 1. In cases where 'stop2' is specified and conflicts with 'stopfend2',
        'stop2' is given preference.
    binbounds2 : 2D numpy array
        An array containing start and stop coordinates for a set of user-defined bins for chrom2. Any fends not
        falling in a bin is ignored. If binbounds1 is specified, start2, stop2, startfend2, and stopfend2 are ignored.
    binsize : int, optional
        This is the coordinate width of each bin.
    datatype : str, optional
        This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fend',
        'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when
        'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-
        filtered bins return value of 1. Expected values are returned for 'distance', 'fend', 'enrichment', and
        'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the
        expected values, 'fend' uses only fend correction values, and both 'enrichment' and 'expected' use both
        correction and distance mean values.
    returnmapping : bool, optional
        If 'True', a list containing the data array and two 2d arrasw of N x 4 containing the first fend and last fend
        plus one included in each bin and first and last coordinates for 'chrom1' and 'chrom2' is return. Otherwise
        only the data array is returned.
    """
    # check that all values are acceptable
    datatypes = {'raw': 0, 'fend': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
    if datatype not in datatypes:
        print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    else:
        datatype_int = datatypes[datatype]
    # Determine start, stop, startfend, and stopfend
    chrint1 = hic.chr2int[chrom1.strip('chr')]
    chrint2 = hic.chr2int[chrom2.strip('chr')]
    if binbounds1 is None:
        if start1 is None and startfend1 is None:
            startfend1 = hic.fends['chr_indices'][chrint1]
            while startfend1 < hic.fends['chr_indices'][chrint1 + 1] and hic.filter[startfend1] == 0:
                startfend1 += 1
            start1 = (hic.fends['fends']['mid'][startfend1] / binsize) * binsize
        elif start1 is None:
            start1 = (hic.fends['fends']['mid'][startfend1] / binsize) * binsize
        else:
            startfend1 = _find_fend_from_coord(hic, chrint1, start1)
        if stop1 is None and stopfend1 is None:
            stopfend1 = hic.fends['chr_indices'][chrint1 + 1]
            while stopfend1 > hic.fends['chr_indices'][chrint1] and hic.filter[stopfend1 - 1] == 0:
                stopfend1 -= 1
            stop1 = int(ceil((hic.fends['fends']['mid'][stopfend1 - 1] + 1 - start1) / float(binsize))) * binsize + start1
        elif stop1 is None:
            stop1 = int(ceil((hic.fends['fends']['mid'][stopfend1 - 1] + 1 - start1) / float(binsize))) * binsize + start1
        else:
            stopfend1 = _find_fend_from_coord(hic, chrint1, stop1)
        num_bins1 = (stop1 - start1) / binsize
    else:
        num_bins1 = binbounds1.shape[0]
        start1 = binbounds1[0, 0]
        stop1 = binbounds1[-1, 1]
        startfend1 = _find_fend_from_coord(hic, chrint1, start1)
        stopfend1 = _find_fend_from_coord(hic, chrint1, stop1)
    if binbounds2 is None:
        if start2 is None and startfend2 is None:
            startfend2 = hic.fends['chr_indices'][chrint2]
            while startfend2 < hic.fends['chr_indices'][chrint2 + 1] and hic.filter[startfend2] == 0:
                startfend2 += 1
            start2 = (hic.fends['fends']['mid'][startfend2] / binsize) * binsize
        elif start2 is None:
            start2 = (hic.fends['fends']['mid'][startfend2] / binsize) * binsize
        else:
            startfend2 = _find_fend_from_coord(hic, chrint2, start2)
        if stop2 is None and stopfend2 is None:
            stopfend2 = hic.fends['chr_indices'][chrint2 + 1]
            while stopfend2 > hic.fends['chr_indices'][chrint2] and hic.filter[stopfend2 - 1] == 0:
                stopfend2 -= 1
            stop2 = int(ceil((hic.fends['fends']['mid'][stopfend2 - 1] + 1 - start2) / float(binsize))) * binsize + start2
        elif stop2 is None:
            stop2 = int(ceil((hic.fends['fends']['mid'][stopfend2 - 1] + 1 - start2) / float(binsize))) * binsize + start2
        else:
            stopfend2 = _find_fend_from_coord(hic, chrint2, stop2)
        num_bins2 = (stop2 - start2) / binsize
    else:
        num_bins2 = binbounds2.shape[0]
        start2 = binbounds2[0, 0]
        stop2 = binbounds2[-1, 1]
        startfend2 = _find_fend_from_coord(hic, chrint2, start2)
        stopfend2 = _find_fend_from_coord(hic, chrint2, stop2)
    # If trans mean not already in hic and 'distance', 'enrichment', or 'expected' is requested, calculate
    if datatype in ['raw', 'fend']:
        trans_mean = 0.0
    elif 'trans_mean' not in hic.__dict__.keys():
        hic.find_trans_mean()
        trans_mean = hic.trans_mean
    else:
        trans_mean = hic.trans_mean
    print >> sys.stderr, ("Finding %s array for %s:%i-%i by %s:%i-%i...") %\
                         (datatype, chrom1, start1, stop1, chrom2, start2, stop2),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        if startfend1 < startfend2:
            start_index = hic.data['trans_indices'][startfend1]
            stop_index = hic.data['trans_indices'][stopfend1]
            data_indices = hic.data['trans_indices'][startfend1:(stopfend1 + 1)]
        else:
            start_index = hic.data['trans_indices'][startfend2]
            stop_index = hic.data['trans_indices'][stopfend2]
            data_indices = hic.data['trans_indices'][startfend2:(stopfend2 + 1)]
        data_indices -= data_indices[0]
        data = hic.data['trans_data'][start_index:stop_index, :]
    else:
        data_indices = None
        data = None
    # Find fend ranges for each bin
    mids1 = hic.fends['fends']['mid'][startfend1:stopfend1]
    if binbounds1 is None:
        mapping1 = numpy.searchsorted(numpy.arange(1, num_bins1 + 1) * binsize + start1, mids1).astype(numpy.int32)
    else:
        mapping1 = numpy.zeros(mids1.shape[0], dtype=numpy.int32) - 1
        starts = numpy.searchsorted(binbounds1[:, 0], mids1, side='right') - 1
        stops = numpy.searchsorted(binbounds1[:, 1], mids1, side='right')
        where = numpy.where(starts == stops)[0]
        mapping1[where] = starts[where].astype(numpy.int32)
    mids2 = hic.fends['fends']['mid'][startfend2:stopfend2]
    if binbounds2 is None:
        mapping2 = numpy.searchsorted(numpy.arange(1, num_bins2 + 1) * binsize + start2, mids2).astype(numpy.int32)
    else:
        mapping2 = numpy.zeros(mids2.shape[0], dtype=numpy.int32) - 1
        starts = numpy.searchsorted(binbounds2[:, 0], mids2, side='right') - 1
        stops = numpy.searchsorted(binbounds2[:, 1], mids2, side='right')
        where = numpy.where(starts == stops)[0]
        mapping2[where] = starts[where].astype(numpy.int32)
    # Create requested array
    if startfend1 < startfend2:
        data_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins2, num_bins1, 2), dtype=numpy.float32)
    # Fill in data values
    if startfend1 < startfend2:
        _binning.binned_signal_trans(data, data_indices, hic.filter, mapping1, mapping2, hic.corrections, data_array,
                                     trans_mean, startfend1, startfend2, datatype_int)
    else:
        _binning.binned_signal_trans(data, data_indices, hic.filter, mapping2, mapping1, hic.corrections, data_array,
                                     trans_mean, startfend2, startfend1, datatype_int)
        data_array = numpy.transpose(data_array, axes=[1, 0, 2])
    # If mapping requested, calculate bin bounds
    if returnmapping:
        mapping1 = numpy.zeros((num_bins1, 4), dtype=numpy.int32)
        if binbounds1 is None:
            mapping1[:, 2] = numpy.arange(num_bins1) * binsize + start1
            mapping1[:, 3] = mapping1[:, 2] + binsize
        else:
            mapping1[:, 2:4] = binbounds1[:, :]
        mapping1[:, 0] = numpy.searchsorted(mids1, mapping1[:, 2]) + startfend1
        mapping1[:, 1] = numpy.searchsorted(mids1, mapping1[:, 3]) + startfend1
        mapping2 = numpy.zeros((num_bins2, 4), dtype=numpy.int32)
        if binbounds2 is None:
            mapping2[:, 2] = numpy.arange(num_bins2) * binsize + start2
            mapping2[:, 3] = mapping2[:, 2] + binsize
        else:
            mapping2[:, 2:4] = binbounds2[:, :]
        mapping2[:, 0] = numpy.searchsorted(mids2, mapping2[:, 2]) + startfend2
        mapping2[:, 1] = numpy.searchsorted(mids2, mapping2[:, 3]) + startfend2
        print >> sys.stderr, ("Done\n"),
        return [data_array, mapping1, mapping2]
    else:
        print >> sys.stderr, ("Done\n"),
        return data_array


def write_heatmap_dict(hic, filename, binsize, includetrans=True, removedistance=False, chroms=[]):
    """
    write_heatmap_dict method

    Create an h5dict file containing binned interaction arrays.

    Parameters
    ----------
    hic : HiC class object
        A HiC class object containing fend and count data.
    filename : string
        Location to write h5dict object to.
    binsize : int
        Size of bins for interaction arrays.
    includetrans : bool, optional
        Indicates whether trans interaction arrays should be calculated and saved.
    removedistance : bool, optional
        If 'True', the expected value is calculated including the expected distance mean. Otherwise, only fend
        corrections are used.
    chroms : list, optional
        If given, indicates which chromosomes should be included. If left empty, all chromosomes are included.
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
    # Check if trans mean is needed and calculate if not already done
    if includetrans and removedistance and 'trans_mean' not in hic.__dict__.keys():
        hic.find_trans_mean()
    # Check if filename already exists, and remove if it does
    if rank == 0:
        if os.path.exists(filename):
            print >> sys.stderr, ("%s already exists, overwriting.") % filename
            subprocess.call('rm %s' % filename, shell=True)
        print >> sys.stderr, ("Creating binned heatmap...\n"),
        output = h5py.File(filename, 'w')
        output['resolution'] = binsize
        # If chromosomes not specified, fill list
        if len(chroms) == 0:
            chroms = list(hic.fends['chromosomes'][...])
        # Assemble list of requested arrays
        needed = []
        for chrom in chroms:
            needed.append((chrom,))
        if includetrans:
            for i in range(len(chroms)-1):
                for j in range(i + 1, len(chroms)):
                    needed.append((chroms[i],chroms[j]))
        if num_procs == 1:
            node_needed = needed
        else:
            worker_size = int(ceil(len(needed) / float(num_procs)))
            for i in range(1, num_procs):
                comm.send(needed[(worker_size * (i - 1)):(worker_size * i)], dest=i, tag=11)
            node_needed = needed[(worker_size * (num_procs - 1)):]
    else:
        node_needed = comm.recv(source=0, tag=11)
    # Determine the requested data type
    if removedistance:
        datatype = 'enrichment'
    else:
        datatype = 'fend'
    heatmaps = {}
    # Find heatmaps
    for chrom in node_needed:
        if len(chrom) == 1:
            # Find cis heatmap
            heatmaps[chrom] = bin_cis_signal(hic, chrom[0], binsize=binsize, datatype=datatype, arraytype='upper',
                                             returnmapping=True)
            # Check if array contains data
            if heatmaps[chrom] is None or heatmaps[chrom][0].shape[0] == 0:
                del heatmaps[chrom]
        else:
            # Find trans heatmap
                heatmaps[chrom] = bin_trans_signal(hic, chrom[0], chrom[1], binsize=binsize, datatype=datatype)
    # Collect heatmaps at node 0 and write to h5dict
    if rank == 0:
        if num_procs > 1:
            for i in range(1, num_procs):
                temp = comm.recv(source=i, tag=11)
                heatmaps.update(temp)
            del temp
        for chrom in heatmaps.keys():
            if len(chrom) == 1:
                output.create_dataset('%s.counts' % chrom[0], data=heatmaps[chrom][0][:, 0])
                output.create_dataset('%s.expected' % chrom[0], data=heatmaps[chrom][0][:, 1])
                output.create_dataset('%s.positions' % chrom[0], data=heatmaps[chrom][1][:, 2:4])
            else:
                output.create_dataset('%s_by_%s.counts' % (chrom[0], chrom[1]), data=heatmaps[chrom][:, :, 0])
                output.create_dataset('%s_by_%s.expected' % (chrom[0], chrom[1]), data=heatmaps[chrom][:, :, 1])
        output.create_dataset('chromosomes', data=numpy.array(chroms))
        output.close()
        print >> sys.stderr, ("Creating binned heatmap...Done\n"),
    else:
        comm.send(heatmaps, dest=0, tag=11)
        del heatmaps
    return None
