#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

"""
This is a module contains scripts for generating compact and full matrices of
interaction data.

Input data
----------

These functions rely on the "FiveC" class in conjunction with the "Fragment"
and "FiveCData" classes.

Concepts
--------

Data can either be arranged in compact, complete, or flattened (row-major)
upper-triangle arrays. Compact arrays are N x M, where N is the number of
forward probe fragments and M is the number of reverse probe fragments. Data
can be raw, fragment-corrected, distance-dependence removed, or enrichment
values. Arrays are 3-dimensional with observed values in the first layer of
d3, expected values in the second layer of d3. The exception to this is
upper-triangle arrays, which are 2d, divinding observed and expected along
the second axis.

-----------------------------------------------------------------------------

API documentation
-----------------



"""

import os
import sys
import subprocess

import numpy
import h5py
from math import floor, ceil

import _distance
import _binning


def unbinned_cis_signal(fivec, region, start=None, stop=None, startfrag=None, stopfrag=None, datatype='enrichment',
                        arraytype='compact', skipfiltered=False, returnmapping=False):
    """
    unbinned_cis_signal method

    Create an array of format 'arraytype' and fill with data requested in 'datatype'.

    Parameters
    ----------
    fivec : FiveC class object
        A FiveC class object containing fragment and count data.
    region : int
        The index of the region to return.
    start : int, optional
        The smallest coordinate to include in the array, measured from fragment midpoints. If both 'start' and
        'startfrag' are given, 'start' will override 'startfrag'. If unspecified, this will be set to the midpoint of
        the first fragment for 'region'.
    stop : int, optional
        The largest coordinate to include in the array, measured from fragment midpoints. If both 'stop' and 'stopfrag'
        are given, 'stop' will override 'stopfrag'. If unspecified, this will be set to the midpoint of the last
        fragment plus one for 'region'.
    startfrag : int, optional
        The first fragment to include in the array. If unspecified and 'start' is not given, this is set to the first
        fragment in 'region'. In cases where 'start' is specified and conflicts with 'startfrag', 'start' is given
        preference.
    stopfrag : int, optional
        The first fragment not to include in the array. If unspecified and 'stop' is not given, this is set to the last
        fragment in 'region' plus one. In cases where 'stop' is specified and conflicts with 'stopfrag', 'stop' is
        given preference.
    datatype : str, optional
        This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment',
        'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when
        'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-
        filtered fragments return value of one. Expected values are returned for 'distance', 'fragment', 'enrichment',
        and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating
        the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use
        both correction and distance mean values. 'enrichment' also scales both observed and expected by the standard
        deviation, giving a completely normalized set of values.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'compact', 'full', and 'upper'.
        'compact' means data are arranged in a N x M x 2 array where N and M are the number of forward and reverse
        probe fragments, respectively. 'full' returns a square, symmetric array of size N x N x 2 where N is the total
        number of fragments. 'upper' returns only the flattened upper triangle of a full array, excluding the diagonal
        of size (N * (N - 1) / 2) x 2, where N is the total number of fragments.
    skipfiltered : bool, optional
        If 'True', all interaction bins for filtered out fragments are removed and a reduced-size array is returned.
    returnmapping : bool, optional
        If 'True', a list containing the data array and a 1d array containing fragment numbers included in the data
        array is return. Otherwise only the data array is returned.
    """
    # check that all values are acceptable
    datatypes = {'raw': 0, 'fragment': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
    if datatype not in datatypes:
        print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    else:
        datatype_int = datatypes[datatype]
    if arraytype not in ['full', 'compact', 'upper']:
        print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # if parameters are no already found, calculate now
    if 'gamma' not in fivec.__dict__.keys():
        fivec.find_distance_parameters()
    # Determine start, stop, startfrag, and stopfrag
    chrint = fivec.frags['regions']['chromosome'][region]
    chrom = fivec.frags['chromosomes'][chrint]
    if start is None and startfrag is None:
        startfrag = fivec.frags['regions']['start_frag'][region]
        start = fivec.frags['fragments']['mid'][startfrag]
    elif start is None:
        start = fivec.frags['fragments']['mid'][startfrag]
    else:
        startfrag = _find_frag_from_coord(fivec, chrint, start)
    if stop is None and stopfrag is None:
        stopfrag = fivec.frags['regions']['stop_frag'][region]
        stop = fivec.frags['fragments']['mid'][stopfrag - 1] + 1
    elif stop is None:
        stop = fivec.frags['fragments']['mid'][stopfrag - 1] + 1
    else:
        stopfrag = _find_frag_from_coord(fivec, chrint, stop)
    print >> sys.stderr, ("Finding %s %s array for %s:%i-%i...") % (datatype, arraytype, chrom, start, stop),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        start_index = fivec.data['cis_indices'][startfrag]
        stop_index = fivec.data['cis_indices'][stopfrag]
        data_indices = fivec.data['cis_indices'][startfrag:(stopfrag + 1)]
        data_indices -= data_indices[0]
        data = fivec.data['cis_data'][start_index:stop_index, :]
        data[:, :2] -= startfrag
    else:
        data_indices = None
        data = None
    mids = fivec.frags['fragments']['mid'][startfrag:stopfrag]
    strands = fivec.frags['fragments']['strand'][startfrag:stopfrag]
    temp_strands = fivec.frags['fragments']['strand'][...]
    temp_data = fivec.data['cis_data'][...]
    # Create requested array and corresponding mapping
    if arraytype == 'compact':
        mapping = numpy.zeros(stopfrag - startfrag, dtype=numpy.int32)
        forward = 0
        reverse = 0
        for i in range(mapping.shape[0]):
            if skipfiltered and fivec.filter[i + startfrag] == 0:
                continue
            if fivec.frags['fragments']['strand'][i + startfrag] == 0:
                mapping[i] = forward + 1
                forward += 1
            else:
                mapping[i] = reverse - 1
                reverse -= 1
        if forward == 0 or reverse == 0:
            print >> sys.stderr, ("Insufficient data.\n"),
            return None
        data_array = numpy.zeros((forward, -reverse, 2), dtype=numpy.float32)
    else:
        if skipfiltered:
            num_bins = numpy.sum(fivec.filter[startfrag:stopfrag])
        else:
            num_bins = stopfrag - startfrag
        if num_bins < 2:
            print >> sys.stderr, ("Insufficient data.\n"),
            return None
        mapping = numpy.zeros(stopfrag - startfrag, dtype=numpy.int32) - 1
        i = 0
        for j in range(stopfrag - startfrag):
            if skipfiltered and fivec.filter[j + startfrag] == 0:
                continue
            mapping[j] = i
            i += 1
        data_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    if datatype_int <= 1:
        mu = 0.0
    else:
        mu = fivec.mu
    # Fill in data values
    if arraytype == 'compact':
        _binning.unbinned_signal_compact(data, fivec.filter[startfrag:stopfrag], mapping,
                                         fivec.corrections[startfrag:stopfrag], mids, data_array, mu,
                                         fivec.gamma, fivec.sigma, datatype_int)
    else:
        _binning.unbinned_signal_upper(data, data_indices, fivec.filter[startfrag:stopfrag], strands, mapping,
                                       fivec.corrections[startfrag:stopfrag], mids, data_array, mu, fivec.gamma,
                                       fivec.sigma, num_bins, datatype_int)
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
        if arraytype == 'compact':
            xmapping = numpy.where(mapping > 0)[0] + startfrag
            ymapping = numpy.where(mapping < 0)[0] + startfrag
            return [data_array, xmapping, ymapping]
        else:
            return [data_array, numpy.where(mapping >= 0)[0] + startfrag]
    else:
        print >> sys.stderr, ("Done\n"),
        return data_array


def _find_frag_from_coord(fivec, chrint, coord):
    """Find the next fragment after the coordinate on chromosome 'chrint'."""
    first_frag = fivec.frags['chr_indices'][chrint]
    last_frag = fivec.frags['chr_indices'][chrint + 1]
    return numpy.searchsorted(fivec.frags['fragments']['mid'][first_frag:last_frag], coord) + first_frag


def bin_cis_signal(fivec, region, start=None, stop=None, startfrag=None, stopfrag=None, binsize=10000,
                   datatype='enrichment', arraytype='full', returnmapping=False):
    """
    bin_cis_signal method

    Create an array of format 'arraytype' and fill 'binsize' bins with data requested in 'datatype'.

    Parameters
    ----------
    fivec : FiveC class object
        A FiveC class object containing fragment and count data.
    region : int
        The index of the region to return.
    start : int, optional
        The coordinate at the beginning of the smallest bin. If unspecified, 'start' will be the first multiple of
        'binsize' below the 'startfend' mid. If there is a conflict between 'start' and 'startfrag', 'start' is given
        preference.
    stop : int, optional
        The largest coordinate to include in the array, measured from fend midpoints. If both 'stop' and 'stopfrag'
        are given, 'stop' will override 'stopfrag'.
    startfrag : int, optional
        The first fragment to include in the array. If unspecified and 'start' is not given, this is set to the first
        valid fragment in 'region'. In cases where 'start' is specified and conflicts with 'startfrag', 'start' is
        given preference.
    stopfrag : int, optional
        The first fragment not to include in the array. If unspecified and 'stop' is not given, this is set to the last
        valid fragment in 'region' + 1. In cases where 'stop' is specified and conflicts with 'stopfrad', 'stop' is
        given preference.
    binsize : int, optional
        This is the coordinate width of each bin.
    datatype : str, optional
        This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment',
        'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when
        'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-
        filtered expected bins return value of 1. Expected values are returned for 'distance', 'fragment',
        'enrichment', and 'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for
        calculating the expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and
        'expected' use both correction and distance values.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'full' and 'upper'. 'full'
        returns a square, symmetric array of size N x N x 2 where N is the total number of fragments. 'upper' returns
        only the flattened upper triangle of a full array, excluding the diagonal of size (N * (N - 1) / 2) x 2, where
        N is the total number of fragments.
    returnmapping : bool, optional
        If 'True', a list containing the data array and a 2d array of N x 4 containing the first fend and last fend
        plus one included in each bin and first and last coordinates is return. Otherwise only the data array is
        returned.
    """
    # check that all values are acceptable
    datatypes = {'raw': 0, 'fragment': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
    if datatype not in datatypes:
        print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    else:
        datatype_int = datatypes[datatype]
    if arraytype not in ['full', 'upper']:
        print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # if parameters are no already found, calculate now
    if 'gamma' not in fivec.__dict__.keys():
        fivec.find_distance_parameters()
    # Determine start, stop, startfrag, and stopfrag
    chrint = fivec.frags['regions']['chromosome'][region]
    chrom = fivec.frags['chromosomes'][chrint]
    if start is None and startfrag is None:
        startfrag = fivec.frags['regions']['start_frag'][region]
        while startfrag < fivec.frags['regions']['stop_frag'][region] and fivec.filter[startfrag] == 0:
            startfrag += 1
        if startfrag == fivec.frags['regions']['stop_frag'][region]:
            print >> sys.stderr, ("Insufficient data.\n"),
            return None
        start = (fivec.frags['fragments']['mid'][startfrag] / binsize) * binsize
    elif start is None:
        start = (fivec.frags['fragments']['mid'][startfrag] / binsize) * binsize
    else:
        startfrag = _find_frag_from_coord(fivec, chrint, start)
    if stop is None and stopfrag is None:
        stopfrag = fivec.frags['regions']['stop_frag'][region]
        while stopfrag > startfrag and fivec.filter[stopfrag - 1] == 0:
            stopfrag -= 1
        stop = (int(ceil((fivec.frags['fragments']['mid'][stopfrag - 1] + 1 - start) / float(binsize))) *
                binsize + start)
    elif stop is None:
        stop = (int(ceil((fivec.frags['fragments']['mid'][stopfrag - 1] + 1 - start) / float(binsize))) *
                binsize + start)
    else:
        stopfrag = _find_frag_from_coord(fivec, chrint, stop)
    num_bins = (stop - start) / binsize
    num_frags = stopfrag - startfrag
    print >> sys.stderr, ("Finding %s %s array for %s:%i-%i...") % (datatype, arraytype, chrom, start, stop),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        start_index = fivec.data['cis_indices'][startfrag]
        stop_index = fivec.data['cis_indices'][stopfrag]
        data_indices = fivec.data['cis_indices'][startfrag:(stopfrag + 1)]
        data_indices -= data_indices[0]
        data = fivec.data['cis_data'][start_index:stop_index, :]
        data[:, :2] -= startfrag
    else:
        data_indices = None
        data = None
    mids = fivec.frags['fragments']['mid'][startfrag:stopfrag]
    strands = fivec.frags['fragments']['strand'][startfrag:stopfrag]    
    # Find fragment ranges for each bin
    mapping = numpy.searchsorted(numpy.arange(1, num_bins + 1) * binsize + start, mids).astype(numpy.int32)
    # Create requested array
    data_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in data values
    _binning.binned_signal_upper(data, data_indices, fivec.filter[startfrag:stopfrag], mapping,
                                 fivec.corrections[startfrag:stopfrag], mids, strands, data_array,
                                 fivec.mu, fivec.gamma, fivec.sigma, datatype_int, num_bins)
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
        mapping[:, 0] = numpy.searchsorted(mids, mapping[:, 2]) + startfrag
        mapping[:, 1] = numpy.searchsorted(mids, mapping[:, 3]) + startfrag
        print >> sys.stderr, ("Done\n"),
        return [data_array, mapping]
    else:
        print >> sys.stderr, ("Done\n"),
        return data_array


def bin_cis_array(fivec, unbinned, fragments, start=None, stop=None, binsize=10000, binbounds=None, arraytype='full',
                  returnmapping=False):
    """
    bin_cis_array method

    Create an array of format 'arraytype' and fill 'binsize' bins or bins defined by 'binbounds' with data provided in
    'unbinned'.

    Parameters
    ----------
    fivec : FiveC class object
        A FiveC class object containing fragment and count data.
    unbinned : numpy array of 'upper' or 'full' format
        A 2d or 3d array containing data to be binned. Array format will be determined from the number of dimensions.
    fragments : 1d numpy array
        An integer array indicating which position corresponds to which fragment in the 'unbinned' array.
    start : int, optional
        The coordinate at the beginning of the smallest bin. If unspecified, 'start' will be the first multiple of
        'binsize' below the first mid from 'fragments'. If 'binbounds' is given, 'start' is ignored.
    stop : int, optional
        The coordinate at the end of the last bin. If unspecified, 'stop' will be the first multiple of 'binsize'
        above the last mid from 'fragments'. If needed, 'stop' is adjusted upward to create a complete last bin. If
        'binbounds' is given, 'stop' is ignored.
    binsize : int, optional
        This is the coordinate width of each bin. This is ignored if 'binbounds' is given.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'full' and 'upper'. 'full'
        returns a square, symmetric array of size N x N x 2. 'upper' returns only the flattened upper triangle of a
        full array, excluding the diagonal of size (N * (N - 1) / 2) x 2.
    returnmapping : bool, optional
        If 'True', a list containing the data array and a 2d array of N x 4 containing the first fragment and last
        fragment plus one included in each bin and first and last coordinates is return. Otherwise only the data array
        is returned.
    """
    # check that arraytype value is acceptable
    if arraytype not in ['full', 'upper']:
        print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # Determine input array type
    if len(unbinned.shape) == 2 and fragments.shape[0] * (fragments.shape[0] - 1) / 2 == unbinned.shape[0]:
        input_type = 'upper'
        ub_signal = unbinned
    elif len(unbinned.shape) == 3 and unbinned.shape[0] == fragments.shape[0]:
        input_type = 'full'
        ub_signal = numpy.zeros((unbinned.shape[0] * (unbinned.shape[0] - 1) / 2, 2), dtype=numpy.float32)
        indices = numpy.triu_indices(unbinned.shape[0], 1)
        ub_signal[:, 0] = unbinned[indices[0], indices[1], 0]
        ub_signal[:, 1] = unbinned[indices[0], indices[1], 1]
    else:
        print >> sys.stderr, ("Unrecognized input array type. No data returned.\n"),
        return None
    # Determine start and stop, if necessary
    if binbounds is None:
        if start is None:
            start = (fivec.frags['fragments']['mid'][fragments[0]] / binsize) * binsize
        if stop is None:
            stop = (int(ceil((fivec.frags['fragments']['mid'][fragments[-1]] + 1 - start) / float(binsize))) *
                    binsize + start)
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
    # Find bin mapping for each fragment
    mapping = numpy.zeros(fragments.shape[0], dtype=numpy.int32) - 1
    mids = fivec.frags['fragments']['mid'][fragments]
    for i in range(binbounds.shape[0]):
        firstfrag = numpy.searchsorted(mids, binbounds[i, 0])
        lastfrag = numpy.searchsorted(mids, binbounds[i, 1])
        mapping[firstfrag:lastfrag] = i
    # Create requested array
    binned_array = numpy.zeros((num_bins * (num_bins - 1) / 2, 2), dtype=numpy.float32)
    # Fill in binned data values
    _binning.bin_upper_to_upper(binned_array, ub_signal, mapping, num_bins)
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
        mapping[:, 0] = numpy.r_[fragments, fragments[-1] + 1][numpy.searchsorted(mids, mapping[:, 2])]
        mapping[:, 1] = numpy.r_[fragments, fragments[-1] + 1][numpy.searchsorted(mids, mapping[:, 3])]
        print >> sys.stderr, ("Done\n"),
        return [binned_array, mapping]
    else:
        print >> sys.stderr, ("Done\n"),
        return binned_array


def dynamically_bin_cis_array(unbinned, unbinnedpositions, binned, binbounds, minobservations=50,
                              searchdistance=0):
    """
    dynamically_bin_cis_array method

    Expand bins in 'binned' to include additional data provided in 'unbinned' as necessary to meet 'minobservations',
    or 'searchdistance' criteria.

    Parameters
    ----------
    unbinned : numpy array of either 'full' or 'upper' format
        A 2d or 3d array containing data to be binned. Array format will be determined from the number of dimensions.
    unbinnedpositions : 1d numpy array
        An integer array indicating the mid-point of each bin in 'unbinned' array.
    binned : numpy array of either 'full' or 'upper' format
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
        ub_signal = unbinned
    elif len(unbinned.shape) == 3 and unbinned.shape[0] == unbinnedpositions.shape[0]:
        unbinned_type = 'full'
        ub_signal = numpy.zeros((unbinned.shape[0] * (unbinned.shape[0] - 1) / 2, 2), dtype=numpy.float32)
        indices = numpy.triu_indices(unbinned.shape[0], 1)
        ub_signal[:, 0] = unbinned[indices[0], indices[1], 0]
        ub_signal[:, 1] = unbinned[indices[0], indices[1], 1]
    else:
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
    _binning.dynamically_bin_upper_from_upper(ub_signal, unbinnedpositions, b_signal, binedges,
                                              mids, minobservations, searchdistance)
    if binned_type == 'full':
        binned[indices[0], indices[1], 0] = b_signal[:, 0]
        binned[indices[0], indices[1], 1] = b_signal[:, 1]
        binned[indices[1], indices[0], 0] = b_signal[:, 0]
        binned[indices[1], indices[0], 1] = b_signal[:, 1]
    print >> sys.stderr, ("Done\n"),
    return None


def unbinned_trans_signal(fivec, region1, region2, start1=None, stop1=None, startfrag1=None, stopfrag1=None,
                          start2=None, stop2=None, startfrag2=None, stopfrag2=None, datatype='enrichment',
                          arraytype='full', skipfiltered=False, returnmapping=False):
    """
    unbinned_cis_signal method

    Create an array of format 'arraytype' and fill 'binsize' bins with data requested in 'datatype'.

    Parameters
    ----------
    fivec : FiveC class object
        A FiveC class object containing fragment and count data.
    region1 : int
        The index of the first region to pull data from.
    region2 : int
        The index of the second region to pull data from.
    start1 : int, optional
        The coordinate at the beginning of the smallest bin from 'region1'. If unspecified, 'start1' will be the first
        multiple of 'binsize' below the 'startfrag1' mid. If there is a conflict between 'start1' and 'startfrag1',
        'start1' is given preference.
    stop1 : int, optional
        The largest coordinate to include in the array from 'region1', measured from fragment midpoints. If both
        'stop1' and 'stopfrag1' are given, 'stop1' will override 'stopfrag1'. 'stop1' will be shifted higher as needed
        to make the last bin of size 'binsize'.
    startfrag1 : int, optional
        The first fragment from 'region1' to include in the array. If unspecified and 'start1' is not given, this is
        set to the first valid fend in 'region1'. In cases where 'start1' is specified and conflicts with
        'startfrag1', 'start1' is given preference.
    stopfrag1 : int, optional
        The first fragment not to include in the array from 'region1'. If unspecified and 'stop1' is not given, this
        is set to the last valid fragment in 'region1' + 1. In cases where 'stop1' is specified and conflicts with
        'stopfrag1', 'stop1' is given preference.
    start2 : int, optional
        The coordinate at the beginning of the smallest bin from 'region2'. If unspecified, 'start2' will be the first
        multiple of 'binsize' below the 'startfrag2' mid. If there is a conflict between 'start2' and 'startfrag2',
        'start2' is given preference.
    stop2 : int, optional
        The largest coordinate to include in the array from 'region2', measured from fragment midpoints. If both
        'stop2' and 'stopfrag2' are given, 'stop2' will override 'stopfrag2'. 'stop2' will be shifted higher as needed
        to make the last bin of size 'binsize'.
    startfrag2 : int, optional
        The first fragment from 'region2' to include in the array. If unspecified and 'start2' is not given, this is
        set to the first valid fragment in 'region2'. In cases where 'start2' is specified and conflicts with
        'startfrag2', 'start2' is given preference.
    stopfrag2 : int, optional
        The first fragment not to include in the array from 'region2'. If unspecified and 'stop2' is not given, this
        is set to the last valid fragment in 'region2' + 1. In cases where 'stop2' is specified and conflicts with
        'stopfrag2', 'stop2' is given preference.
    datatype : str, optional
        This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment',
        'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when
        'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-
        filtered bins return value of 1. Expected values are returned for 'distance', 'fragment', 'enrichment', and
        'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the
        expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both
        correction and distance mean values.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'compact' and 'full'.
        'compact' means data are arranged in a N x M x 2 array where N and M are the number of forward and reverse
        probe fragments, respectively. Two arrays will be returned for this format, the first with forward probe
        fragments from region1 and reverse probe fragments from region2. The second is the compliment of the first.
        'full' returns a square, symmetric array of size N x N x 2 where N is the total number of fragments.
    skipfiltered : bool, optional
        If 'True', all interaction bins for filtered out fragments are removed and a reduced-size array is returned.
    returnmapping : bool, optional
        If 'True', a list containing the data array and two 2d arrasw of N x 4 containing the first fragment and last
        fragment plus one included in each bin and first and last coordinates for 'region1' and 'region2' is return.
        Otherwise only the data array is returned.
    """
    # check that all values are acceptable
    datatypes = {'raw': 0, 'fragment': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
    if datatype not in datatypes:
        print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    else:
        datatype_int = datatypes[datatype]
    if arraytype not in ['full', 'compact']:
        print >> sys.stderr, ("Unrecognized array type. No data returned.\n"),
        return None
    # if parameters are no already found, calculate now
    if 'gamma' not in fivec.__dict__.keys():
        fivec.find_distance_parameters()
    # Determine start, stop, startfend, and stopfend
    chrint1 = fivec.frags['regions']['chromosome'][region1]
    chrom1 = fivec.frags['chromosomes'][chrint1]
    chrint2 = fivec.frags['regions']['chromosome'][region2]
    chrom2 = fivec.frags['chromosomes'][chrint2]
    if start1 is None and startfrag1 is None:
        startfrag1 = fivec.frags['regions']['start_frag'][region1]
        start1 = fivec.frags['fragments']['mid'][startfrag1]
    elif start1 is None:
        start1 = fivec.frags['fragments']['mid'][startfrag1]
    else:
        startfrag1 = _find_frag_from_coord(fivec, chrint1, start1)
    if stop1 is None and stopfrag1 is None:
        stopfrag1 = fivec.frags['regions']['stop_frag'][region1]
        stop1 = fivec.frags['fragments']['mid'][stopfrag1 - 1] + 1
    elif stop1 is None:
        stop1 = fivec.frags['fragments']['mid'][stopfrag1 - 1] + 1
    else:
        stopfrag1 = _find_frag_from_coord(fivec, chrint1, stop1)
    if start2 is None and startfrag2 is None:
        startfrag2 = fivec.frags['regions']['start_frag'][region2]
        start2 = fivec.frags['fragments']['mid'][startfrag2]
    elif start2 is None:
        start2 = fivec.frags['fragments']['mid'][startfrag2]
    else:
        startfrag2 = _find_frag_from_coord(fivec, chrint2, start2)
    if stop2 is None and stopfrag2 is None:
        stopfrag2 = fivec.frags['regions']['stop_frag'][region2]
        stop2 = fivec.frags['fragments']['mid'][stopfrag2 - 1] + 1
    elif stop2 is None:
        stop2 = fivec.frags['fragments']['mid'][stopfrag2 - 1] + 1
    else:
        stopfrag2 = _find_frag_from_coord(fivec, chrint2, stop2)
    # If trans mean not already in fivec and 'distance', 'enrichment', or 'expected' is requested, calculate
    if datatype in ['raw', 'fragment']:
        trans_mean = 0.0
    elif 'trans_mean' not in fivec.__dict__.keys():
        fivec.find_trans_mean()
        trans_mean = fivec.trans_mean
    else:
        trans_mean = fivec.trans_mean
    print >> sys.stderr, ("Finding %s array for %s:%i-%i by %s:%i-%i...") %\
                         (datatype, chrom1, start1, stop1, chrom2, start2, stop2),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        if startfrag1 < startfrag2:
            start_index = fivec.data['trans_indices'][startfrag1]
            stop_index = fivec.data['trans_indices'][stopfrag1]
        else:
            start_index = fivec.data['trans_indices'][startfrag2]
            stop_index = fivec.data['trans_indices'][stopfrag2]
        data = fivec.data['trans_data'][start_index:stop_index, :]
    else:
        data = None
    strands = fivec.frags['fragments']['strand'][...]
    # Create requested array and corresponding mapping
    if arraytype == 'compact':
        mapping1 = numpy.zeros(stopfrag1 - startfrag1, dtype=numpy.int32)
        mapping2 = numpy.zeros(stopfrag2 - startfrag2, dtype=numpy.int32)
        forward1 = 0
        reverse1 = 0
        forward2 = 0
        reverse2 = 0
        for i in range(mapping1.shape[0]):
            if skipfiltered and fivec.filter[i + startfrag1] == 0:
                continue
            if fivec.frags['fragments']['strand'][i + startfrag1] == 0:
                mapping1[i] = forward1 + 1
                forward1 += 1
            else:
                mapping1[i] = reverse1 - 1
                reverse1 -= 1
        for i in range(mapping2.shape[0]):
            if skipfiltered and fivec.filter[i + startfrag2] == 0:
                continue
            if fivec.frags['fragments']['strand'][i + startfrag2] == 0:
                mapping2[i] = forward2 + 1
                forward2 += 1
            else:
                mapping2[i] = reverse2 - 1
                reverse2 -= 1
        if forward1 == 0 or reverse1 == 0 or forward2 == 0 or reverse2 == 0:
            print >> sys.stderr, ("Insufficient data.\n"),
            return None
        data_array1 = numpy.zeros((forward1, -reverse2, 2), dtype=numpy.float32)
        data_array2 = numpy.zeros((forward2, -reverse1, 2), dtype=numpy.float32)
    else:
        if skipfiltered:
            num_bins1 = numpy.sum(fivec.filter[startfrag1:stopfrag1])
            num_bins2 = numpy.sum(fivec.filter[startfrag2:stopfrag2])
        else:
            num_bins1 = stopfrag1 - startfrag1
            num_bins2 = stopfrag2 - startfrag2
        if num_bins1 < 1 or num_bins2 < 1:
            print >> sys.stderr, ("Insufficient data.\n"),
            return None
        mapping1 = numpy.zeros(stopfrag1 - startfrag1, dtype=numpy.int32) - 1
        mapping2 = numpy.zeros(stopfrag2 - startfrag2, dtype=numpy.int32) - 1
        i = 0
        for j in range(stopfrag1 - startfrag1):
            if skipfiltered and fivec.filter[j + startfrag1] == 0:
                continue
            mapping1[j] = i
            i += 1
        i = 0
        for j in range(stopfrag2 - startfrag2):
            if skipfiltered and fivec.filter[j + startfrag2] == 0:
                continue
            mapping2[j] = i
            i += 1
        if startfrag1 < startfrag2:
            data_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
        else:
            data_array = numpy.zeros((num_bins2, num_bins1, 2), dtype=numpy.float32)
    # Fill in data values
    if arraytype == 'full':
        if startfrag1 < startfrag2:
            _binning.unbinned_signal_trans_full(data, fivec.filter, fivec.corrections, strands, mapping1, mapping2,
                                                data_array, trans_mean, fivec.sigma, startfrag1, startfrag2,
                                                datatype_int)
        else:
            _binning.unbinned_signal_trans_full(data, fivec.filter, fivec.corrections, strands, mapping2, mapping1,
                                                data_array, trans_mean, fivec.sigma, startfrag2, startfrag1,
                                                datatype_int)
            data_array = numpy.transpose(data_array, axes=[1, 0, 2])
    else:
        if startfrag1 < startfrag2:
            _binning.unbinned_signal_trans_compact(data, fivec.filter, fivec.corrections, strands, mapping1, mapping2,
                                                   data_array1, data_array2, trans_mean, fivec.sigma, startfrag1,
                                                   startfrag2, datatype_int)
        else:
            _binning.unbinned_signal_trans_compact(data, fivec.filter, fivec.corrections, strands, mapping2, mapping1,
                                                   data_array2, data_array1, trans_mean, fivec.sigma, startfrag2,
                                                   startfrag1, datatype_int)
    # If mapping requested, calculate bin bounds
    print >> sys.stderr, ("Done\n"),
    if returnmapping:
        if arraytype == 'full':
            if startfrag1 < startfrag2:
                return [data_array, numpy.where(mapping1 >= 0)[0] + startfrag1,
                        numpy.where(mapping2 >= 0)[0] + startfrag2]
            else:
                return [data_array, numpy.where(mapping2 >= 0)[0] + startfrag2,
                        numpy.where(mapping1 >= 0)[0] + startfrag1]
        else:
            xmapping1 = numpy.where(mapping1 > 0)[0] + startfrag1
            ymapping2 = numpy.where(mapping1 < 0)[0] + startfrag1
            xmapping2 = numpy.where(mapping2 > 0)[0] + startfrag2
            ymapping1 = numpy.where(mapping2 < 0)[0] + startfrag2
            if startfrag1 < startfrag2:
                return [data_array1, data_array2, xmapping1, ymapping1, ymapping2, xmapping2]
            else:
                return [data_array2, data_array1, xmapping2, ymapping2, ymapping1, xmapping1]
    else:
        if arraytype == 'full':
            return data_array
        elif startfrag1 < startfrag2:
            return [data_array1, data_array2]
        else:
            return [data_array2, data_array1]


def bin_trans_signal(fivec, region1, region2, start1=None, stop1=None, startfrag1=None, stopfrag1=None, start2=None,
                     stop2=None, startfrag2=None, stopfrag2=None, binsize=1000000, datatype='enrichment',
                     returnmapping=False):
    """
    bin_cis_signal method

    Create an array and fill 'binsize' bins with data requested in 'datatype'.

    Parameters
    ----------
    fivec : FiveC class object
        A FiveC class object containing fragment and count data.
    region1 : int
        The index of the first region to pull data from.
    region2 : int
        The index of the second region to pull data from.
    start1 : int, optional
        The coordinate at the beginning of the smallest bin from 'region1'. If unspecified, 'start1' will be the first
        multiple of 'binsize' below the 'startfrag1' mid. If there is a conflict between 'start1' and 'startfrag1',
        'start1' is given preference.
    stop1 : int, optional
        The largest coordinate to include in the array from 'region1', measured from fragment midpoints. If both
        'stop1' and 'stopfrag1' are given, 'stop1' will override 'stopfrag1'. 'stop1' will be shifted higher as needed
        to make the last bin of size 'binsize'.
    startfrag1 : int, optional
        The first fragment from 'region1' to include in the array. If unspecified and 'start1' is not given, this is
        set to the first valid fend in 'region1'. In cases where 'start1' is specified and conflicts with
        'startfrag1', 'start1' is given preference.
    stopfrag1 : int, optional
        The first fragment not to include in the array from 'region1'. If unspecified and 'stop1' is not given, this
        is set to the last valid fragment in 'region1' + 1. In cases where 'stop1' is specified and conflicts with
        'stopfrag1', 'stop1' is given preference.
    start2 : int, optional
        The coordinate at the beginning of the smallest bin from 'region2'. If unspecified, 'start2' will be the first
        multiple of 'binsize' below the 'startfrag2' mid. If there is a conflict between 'start2' and 'startfrag2',
        'start2' is given preference.
    stop2 : int, optional
        The largest coordinate to include in the array from 'region2', measured from fragment midpoints. If both
        'stop2' and 'stopfrag2' are given, 'stop2' will override 'stopfrag2'. 'stop2' will be shifted higher as needed
        to make the last bin of size 'binsize'.
    startfrag2 : int, optional
        The first fragment from 'region2' to include in the array. If unspecified and 'start2' is not given, this is
        set to the first valid fragment in 'region2'. In cases where 'start2' is specified and conflicts with
        'startfrag2', 'start2' is given preference.
    stopfrag2 : int, optional
        The first fragment not to include in the array from 'region2'. If unspecified and 'stop2' is not given, this
        is set to the last valid fragment in 'region2' + 1. In cases where 'stop2' is specified and conflicts with
        'stopfrag2', 'stop2' is given preference.
    binsize : int, optional
        This is the coordinate width of each bin.
    datatype : str, optional
        This specifies the type of data that is processed and returned. Options are 'raw', 'distance', 'fragment',
        'enrichment', and 'expected'. Observed values are aways in the first index along the last axis, except when
        'datatype' is 'expected'. In this case, filter values replace counts. Conversely, if 'raw' is specified, non-
        filtered bins return value of 1. Expected values are returned for 'distance', 'fragment', 'enrichment', and
        'expected' values of 'datatype'. 'distance' uses only the expected signal given distance for calculating the
        expected values, 'fragment' uses only fragment correction values, and both 'enrichment' and 'expected' use both
        correction and distance mean values.
    returnmapping : bool, optional
        If 'True', a list containing the data array and two 2d arrasw of N x 4 containing the first fragment and last
        fragment plus one included in each bin and first and last coordinates for 'region1' and 'region2' is return.
        Otherwise only the data array is returned.
    """
    # check that all values are acceptable
    datatypes = {'raw': 0, 'fragment': 1, 'distance': 2, 'enrichment': 3, 'expected': 4}
    if datatype not in datatypes:
        print >> sys.stderr, ("Datatype given is not recognized. No data returned\n"),
        return None
    else:
        datatype_int = datatypes[datatype]
    # if parameters are no already found, calculate now
    if 'gamma' not in fivec.__dict__.keys():
        fivec.find_distance_parameters()
    # Determine start, stop, startfend, and stopfend
    chrint1 = fivec.frags['regions']['chromosome'][region1]
    chrom1 = fivec.frags['chromosomes'][chrint1]
    chrint2 = fivec.frags['regions']['chromosome'][region2]
    chrom2 = fivec.frags['chromosomes'][chrint2]
    if start1 is None and startfrag1 is None:
        startfrag1 = fivec.frags['regions']['start_frag'][region1]
        while startfrag1 < fivec.frags['regions']['stop_frag'][region1] and fivec.filter[startfrag1] == 0:
            startfrag1 += 1
        start1 = (fivec.frags['fragments']['mid'][startfrag1] / binsize) * binsize
    elif start1 is None:
        start1 = (fivec.frags['fragments']['mid'][startfrag1] / binsize) * binsize
    else:
        startfrag1 = _find_frag_from_coord(fivec, chrint1, start1)
    if stop1 is None and stopfrag1 is None:
        stopfrag1 = fivec.frags['regions']['stop_frag'][region1]
        while stopfrag1 > startfrag1 and fivec.filter[stopfrag1 - 1] == 0:
            stopfrag1 -= 1
        stop1 = (int(ceil((fivec.frags['fragments']['mid'][stopfrag1 - 1] + 1 - start1) / float(binsize))) *
                 binsize + start1)
    elif stop1 is None:
        stop1 = (int(ceil((fivec.frags['fragments']['mid'][stopfrag1-1] + 1 - start1) / float(binsize))) *
                 binsize + start1)
    else:
        stopfrag1 = _find_frag_from_coord(fivec, chrint1, stop1)
    num_bins1 = (stop1 - start1) / binsize
    if start2 is None and startfrag2 is None:
        startfrag2 = fivec.frags['regions']['start_frag'][region2]
        while startfrag2 < fivec.frags['regions']['stop_frag'][region2] and fivec.filter[startfrag2] == 0:
            startfrag2 += 1
        start2 = (fivec.frags['fragments']['mid'][startfrag2] / binsize) * binsize
    elif start2 is None:
        start2 = (fivec.frags['fragments']['mid'][startfrag2] / binsize) * binsize
    else:
        startfrag2 = _find_frag_from_coord(fivec, chrint2, start2)
    if stop2 is None and stopfrag2 is None:
        stopfrag2 = fivec.frags['regions']['stop_frag'][region2]
        while stopfrag2 > startfrag2 and fivec.filter[stopfrag2 - 1] == 0:
            stopfrag2 -= 1
        stop2 = (int(ceil((fivec.frags['fragments']['mid'][stopfrag2 - 1] + 1 - start2) / float(binsize))) *
                 binsize + start2)
    elif stop2 is None:
        stop2 = (int(ceil((fivec.frags['fragments']['mid'][stopfrag2-1] + 1 - start2) / float(binsize))) *
                 binsize + start2)
    else:
        stopfrag2 = _find_frag_from_coord(fivec, chrint2, stop2)
    num_bins2 = (stop2 - start2) / binsize
    # If trans mean not already in fivec and 'distance', 'enrichment', or 'expected' is requested, calculate
    if datatype in ['raw', 'fragment']:
        trans_mean = 0.0
    elif 'trans_mean' not in fivec.__dict__.keys():
        fivec.find_trans_mean()
        trans_mean = fivec.trans_mean
    else:
        trans_mean = fivec.trans_mean
    print >> sys.stderr, ("Finding %s array for %s:%i-%i by %s:%i-%i...") %\
                         (datatype, chrom1, start1, stop1, chrom2, start2, stop2),
    # Copy needed data from h5dict for faster access
    if datatype != 'expected':
        if startfrag1 < startfrag2:
            start_index = fivec.data['trans_indices'][startfrag1]
            stop_index = fivec.data['trans_indices'][stopfrag1]
        else:
            start_index = fivec.data['trans_indices'][startfrag2]
            stop_index = fivec.data['trans_indices'][stopfrag2]
        data = fivec.data['trans_data'][start_index:stop_index, :]
    else:
        data = None
    strands = fivec.frags['fragments']['strand'][...]
    # Find fragment ranges for each bin
    mids1 = fivec.frags['fragments']['mid'][startfrag1:stopfrag1]
    mapping1 = numpy.searchsorted(numpy.arange(1, num_bins1 + 1) * binsize + start1, mids1).astype(numpy.int32)
    mids2 = fivec.frags['fragments']['mid'][startfrag2:stopfrag2]
    mapping2 = numpy.searchsorted(numpy.arange(1, num_bins2 + 1) * binsize + start2, mids2).astype(numpy.int32)
    # Create requested array
    if startfrag1 < startfrag2:
        data_array = numpy.zeros((num_bins1, num_bins2, 2), dtype=numpy.float32)
    else:
        data_array = numpy.zeros((num_bins2, num_bins1, 2), dtype=numpy.float32)
    # Fill in data values
    if startfrag1 < startfrag2:
        _binning.binned_signal_trans(data, fivec.filter, fivec.corrections, strands, mapping1, mapping2, data_array,
                                     trans_mean, fivec.sigma, startfrag1, startfrag2, datatype_int)
    else:
        _binning.binned_signal_trans(data, fivec.filter, fivec.corrections, strands, mapping2, mapping1, data_array,
                                     trans_mean, fivec.sigma, startfrag2, startfrag1, datatype_int)
        data_array = numpy.transpose(data_array, axes=[1, 0, 2])
    # If mapping requested, calculate bin bounds
    if returnmapping:
        mapping1 = numpy.zeros((num_bins1, 4), dtype=numpy.int32)
        mapping1[:, 2] = numpy.arange(num_bins1) * binsize + start1
        mapping1[:, 3] = mapping1[:, 2] + binsize
        mapping1[:, 0] = numpy.searchsorted(mids1, mapping1[:, 2]) + startfrag1
        mapping1[:, 1] = numpy.searchsorted(mids1, mapping1[:, 3]) + startfrag1
        mapping2 = numpy.zeros((num_bins2, 4), dtype=numpy.int32)
        mapping2[:, 2] = numpy.arange(num_bins2) * binsize + start2
        mapping2[:, 3] = mapping2[:, 2] + binsize
        mapping2[:, 0] = numpy.searchsorted(mids2, mapping2[:, 2]) + startfrag2
        mapping2[:, 1] = numpy.searchsorted(mids2, mapping2[:, 3]) + startfrag2
        print >> sys.stderr, ("Done\n"),
        return [data_array, mapping1, mapping2]
    else:
        print >> sys.stderr, ("Done\n"),
        return data_array


def write_heatmap_dict(fivec, filename, binsize, includetrans=True, removedistance=False, arraytype='full',
                       regions=[]):
    """
    write_heatmap_dict method

    Create an h5dict file containing binned interaction arrays.

    Parameters
    ----------
    fivec : FiveC class object
        A FiveC class object containing fragment and count data.
    filename : string
        Location to write h5dict object to.
    binsize : int
        Size of bins for interaction arrays. If "binsize" is zero, fragment interactions are returned without
        binning.
    includetrans : bool, optional
        Indicates whether trans interaction arrays should be calculated and saved.
    removedistance : bool, optional
        If 'True', the expected value is calculated including the expected distance mean. Otherwise, only fragment
        corrections are used.
    arraytype : str, optional
        This determines what shape of array data are returned in. Acceptable values are 'compact' and 'full'.
        'compact' means data are arranged in a N x M x 2 array where N is the number of bins, M is the maximum number
        of steps between included bin pairs, and data are stored such that bin n,m contains the interaction values
        between n and n + m + 1. 'full' returns a square, symmetric array of size N x N x 2.
    regions : list, optional
        If given, indicates which regions should be included. If left empty, all regions are included.
    """
    # Check if trans mean is needed and calculate if not already done
    if includetrans and removedistance and 'trans_mean' not in fivec.__dict__.keys():
        fivec.find_trans_mean()
    # Check if filename already exists, and remove if it does
    if os.path.exists(filename):
        print >> sys.stderr, ("%s already exists, overwriting.") % filename
        subprocess.call('rm %s' % filename, shell=True)
    print >> sys.stderr, ("Creating binned heatmap...\n"),
    output = h5py.File(filename, 'w')
    # Determine the requested data type
    if removedistance:
        datatype = 'enrichment'
    else:
        datatype = 'fragment'
    # If regions not specified, fill list
    if len(regions) == 0:
        regions = list(numpy.arange(fivec.frags['regions'].shape[0]))
    if binsize > 0:
        output['resolution'] = binsize
    else:
        output['resolution'] = 'fragment'
    # Find cis heatmaps
    remove = []
    for region in regions:
        if binsize == 0:
            results = unbinned_cis_signal(fivec, region, datatype=datatype, arraytype=arraytype, returnmapping=True)
        else:
            results = bin_cis_signal(fivec, region, binsize=binsize, datatype=datatype, returnmapping=True)
        # Check if array contains data
        if results is None or results[0].shape[0] == 0:
            remove.append(region)
            continue
        output.create_dataset('%i.counts' % region, data=results[0][:, :, 0])
        output.create_dataset('%i.expected' % region, data=results[0][:, :, 1])
        if binsize > 0:
            output.create_dataset('%i.positions' % region, data=results[1][:, 2])
        elif arraytype == 'full':
            output.create_dataset('%i.fragments' % region, data=numpy.arange(
                                  fivec.frags['regions']['start_frag'][region],
                                  fivec.frags['regions']['stop_frag'][region]))
        else:
            strands = fivec.frags['fragments']['strand'][fivec.frags['regions']['start_frag'][region]:
                                  fivec.frags['regions']['stop_frag'][region]]
            output.create_dataset('%i.forward_fragments' % region, data=(
                                  numpy.where(strands == 0)[0] + fivec.frags['regions']['start_frag'][region]))
            output.create_dataset('%i.reverse_fragments' % region, data=(
                                  numpy.where(strands == 1)[0] + fivec.frags['regions']['start_frag'][region]))
    for region in remove:
        del regions[regions.index(region)]
    all_regions = fivec.frags['regions'][...]
    output.create_dataset('regions', data=all_regions[regions][...])
    # If requested, find trans heatmaps
    if includetrans:
        for i in range(len(regions)-1):
            for j in range(i + 1, len(regions)):
                if binsize == 0:
                    results = unbinned_trans_signal(fivec, regions[i], regions[j], datatype=datatype,
                                                    arraytype=arraytype)
                else:
                    results = bin_trans_signal(fivec, regions[i], regions[j], binsize=binsize, datatype=datatype)
                if binsize > 0 or arraytype == 'full':
                    output.create_dataset('%s_by_%s.counts' % (regions[i], regions[j]), data=results[:, :, 0])
                    output.create_dataset('%s_by_%s.expected' % (regions[i], regions[j]), data=results[:, :, 1])
                else:
                    output.create_dataset('%s_by_%s.counts' % (regions[i], regions[j]), data=results[0][:, :, 0])
                    output.create_dataset('%s_by_%s.expected' % (regions[i], regions[j]), data=results[0][:, :, 1])
                    output.create_dataset('%s_by_%s.counts' % (regions[j], regions[i]), data=results[1][:, :, 0])
                    output.create_dataset('%s_by_%s.expected' % (regions[j], regions[i]), data=results[1][:, :, 1])
    output.close
    print >> sys.stderr, ("Creating binned heatmap...Done\n"),
    return None