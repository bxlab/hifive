# distutils: language = c++

"""These functions provide increased speed in handling the signal-binning
functions necessary for supporting hic analysis using HiFive.
"""

import cython
cimport numpy as np
import numpy

ctypedef np.float32_t DTYPE_t
ctypedef np.float64_t DTYPE_64_t
ctypedef np.int32_t DTYPE_int_t
ctypedef np.int64_t DTYPE_int64_t
ctypedef np.uint32_t DTYPE_uint_t
ctypedef np.int8_t DTYPE_int8_t
cdef double Inf = numpy.inf

cdef extern from "math.h":
    double exp(double x) nogil
    double log(double x) nogil
    double log2(double x) nogil
    double log10(double x) nogil
    double sqrt(double x) nogil
    double pow(double x, double x) nogil
    double abs(double x) nogil
    double round(double x) nogil
    double floor(double x) nogil
    double ceil(double x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_compact_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        np.ndarray[DTYPE_int_t, ndim=2] ranges,
        np.ndarray[DTYPE_t, ndim=2] overlap,
        int maxdistance,
        int diag):
    cdef long long int fend1, fend2, i, j, k, map1, map2, start1, start2, stop1, stop2
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        for fend1 in range(num_fends - 1 + diag):
            map1 = mapping[fend1]
            if map1 == -1:
                continue
            for i in range(indices[fend1], indices[fend1 + 1]):
                fend2 = data[i, 1]
                if fend2 >= num_fends:
                    continue
                map2 = mapping[fend2]
                if map2 == -1:
                    continue
                if mids[fend2] - mids[fend1] > maxdistance:
                    continue
                if diag == 0 and map1 == map2:
                    continue
                if ranges is None:
                    signal[map1, map2 - map1 - 1 + diag, 0] += data[i, 2]
                else:
                    start1 = ranges[fend1, 0]
                    start2 = ranges[fend2, 0]
                    stop1 = ranges[fend1, 1]
                    stop2 = ranges[fend2, 1]
                    if start1 < start2:
                        signal[start1, start2 - start1 - 1 + diag, 0] += data[i, 2] * overlap[fend1, 0] * overlap[fend2, 0]
                    if start1 < stop2:
                        signal[start1, stop2 - start1 - 1 + diag, 0] += data[i, 2] * overlap[fend1, 0] * overlap[fend2, 1]
                    for j in range(max(start1, start2) + 1, stop2):
                        signal[start1, j - start1 - 1 + diag, 0] += data[i, 2] * overlap[fend1, 0]
                    if stop1 < start2:
                        signal[stop1, start2 - stop1 - 1 + diag, 0] += data[i, 2] * overlap[fend1, 1] * overlap[fend2, 0]
                    if stop1 < stop2:
                        signal[stop1, stop2 - stop1 - 1 + diag, 0] += data[i, 2] * overlap[fend1, 1] * overlap[fend2, 1]
                    for j in range(max(stop1, start2) + 1, stop2):
                        signal[stop1, j - stop1 - 1 + diag, 0] += data[i, 2] * overlap[fend1, 1]
                    for j in range(start1 + 1, stop1):
                        if j < start2:
                            signal[j, start2 - j - 1 + diag, 0] += data[i, 2] * overlap[fend2, 0]
                        if j < stop2:
                            signal[j, stop2 - j - 1 + diag, 0] += data[i, 2] * overlap[fend2, 1]
                        for k in range(max(j, start2) + 1, stop2):
                            signal[j, k - j - 1 + diag, 0] += data[i, 2]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_upper_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        np.ndarray[DTYPE_int_t, ndim=2] ranges,
        np.ndarray[DTYPE_t, ndim=2] overlap,
        int maxdistance,
        int diag):
    cdef long long int fend1, fend2, i, j, k, index, map1, map2, start1, start2, stop1, stop2, index1
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5)) - diag
    cdef int diag2 = diag * 2
    with nogil:
        for fend1 in range(num_fends - 1 + diag):
            map1 = mapping[fend1]
            if map1 == -1:
                continue
            index = map1 * (num_bins - 1) - map1 * (map1 + 1 - diag2) / 2 - 1 + diag
            for i in range(indices[fend1], indices[fend1 + 1]):
                fend2 = data[i, 1]
                if fend2 >= num_fends:
                    continue
                map2 = mapping[fend2]
                if map2 == -1:
                    continue
                if mids[fend2] - mids[fend1] > maxdistance:
                    continue
                if diag == 0 and map1 == map2:
                    continue
                if ranges is None:
                    signal[index + map2, 0] += data[i, 2]
                else:
                    start1 = ranges[fend1, 0]
                    start2 = ranges[fend2, 0]
                    stop1 = ranges[fend1, 1]
                    stop2 = ranges[fend2, 1]
                    index1 = start1 * (num_bins - 1) - start1 * (start1 + 1 - diag2) / 2 - 1 + diag
                    if start1 < start2 + diag:
                        signal[index1 + start2, 0] += data[i, 2] * overlap[fend1, 0] * overlap[fend2, 0]
                    if start1 < stop2 + diag:
                        signal[index1 + stop2, 0] += data[i, 2] * overlap[fend1, 0] * overlap[fend2, 1]
                    for j in range(max(start1, start2) + 1, stop2):
                        signal[index1 + j, 0] += data[i, 2] * overlap[fend1, 0]
                    index1 = stop1 * (num_bins - 1) - stop1 * (stop1 + 1 - diag2) / 2 - 1 + diag
                    if stop1 < start2 + diag:
                        signal[index1 + start2, 0] += data[i, 2] * overlap[fend1, 1] * overlap[fend2, 0]
                    if stop1 < stop2 + diag:
                        signal[index1 + stop2, 0] += data[i, 2] * overlap[fend1, 1] * overlap[fend2, 1]
                    for j in range(max(stop1, start2) + 1, stop2):
                        signal[index1 + j, 0] += data[i, 2] * overlap[fend1, 1]
                    for j in range(start1 + 1, stop1):
                        index1 = j * (num_bins - 1) - j * (j + 1 - diag2) / 2 - 1 + diag
                        if j < start2 + diag:
                            signal[index1 + start2, 0] += data[i, 2] * overlap[fend2, 0]
                        if j < stop2 + diag:
                            signal[index1 + stop2, 0] += data[i, 2] * overlap[fend2, 1]
                        for k in range(max(j, start2) + 1, stop2):
                            signal[index1 + k, 0] += data[i, 2]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_compact_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        np.ndarray[DTYPE_int_t, ndim=2] ranges,
        np.ndarray[DTYPE_t, ndim=2] overlap,
        double chrom_mean,
        int startfend,
        int maxdistance,
        int diag):
    cdef long long int fend1, fend2, afend1, afend2, j, k, map1, map2, index, num_parameters, num_bins, max_bin
    cdef long long int start1, start2, stop1, stop2, l, m
    cdef double distance, value
    cdef long long int num_fends = mapping.shape[0]
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    else:
        num_parameters = 0
    if not correction_sums is None:
        num_bins = correction_sums.shape[0]
        max_bin = signal.shape[1]
    with nogil:
        if correction_sums is None or not ranges is None:
            for fend1 in range(num_fends - 2):
                map1 = mapping[fend1]
                if map1 == -1:
                    continue
                k = 0
                # find opposite strand adjacents, skipping same fragment and same strand adjacents
                fend2 = fend1 + 2
                map2 = mapping[fend2]
                if map2 >= 0 and (diag == 1 or map2 != map1) and mids[fend2] - mids[fend1] <= maxdistance:
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections is None:
                        value *= corrections[fend1] * corrections[fend2]
                    # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                    if not binning_corrections is None:
                        for j in range(num_parameters):
                            afend1 = fend1 + startfend
                            afend2 = fend2 + startfend
                            if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                                value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                            else:
                                value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                    # if finding distance, enrichment, or expected, correct for distance
                    if not parameters is None:
                        distance = log(<double>(mids[fend2] - mids[fend1]))
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    if ranges is None:
                        signal[map1, map2 - map1 - 1 + diag, 1] += value
                    else:
                        start1 = ranges[fend1, 0]
                        start2 = ranges[fend2, 0]
                        stop1 = ranges[fend1, 1]
                        stop2 = ranges[fend2, 1]
                        if start1 < start2 + diag:
                            signal[start1, start2 - start1 - 1 + diag, 1] += value * overlap[fend1, 0] * overlap[fend2, 0]
                        if start1 < stop2 + diag:
                            signal[start1, stop2 - start1 - 1 + diag, 1] += value * overlap[fend1, 0] * overlap[fend2, 1]
                        for l in range(max(start1, start2) + 1, stop2):
                            signal[start1, l - start1 - 1 + diag, 1] += value * overlap[fend1, 0]
                        if stop1 < start2 + diag:
                            signal[stop1, start2 - stop1 - 1 + diag, 1] += value * overlap[fend1, 1] * overlap[fend2, 0]
                        if stop1 < stop2 + diag:
                            signal[stop1, stop2 - stop1 - 1 + diag, 1] += value * overlap[fend1, 1] * overlap[fend2, 1]
                        for l in range(max(stop1, start2) + 1, stop2):
                            signal[stop1, l - stop1 - 1 + diag, 1] += value * overlap[fend1, 1]
                        for l in range(start1 + 1, stop1):
                            if l < start2 + diag:
                                signal[l, start2 - l - 1 + diag, 1] += value * overlap[fend2, 0]
                            if l < stop2 + diag:
                                signal[l, stop2 - l - 1 + diag, 1] += value * overlap[fend2, 1]
                            for m in range(max(l, start2) + 1, stop2):
                                signal[l, m - l - 1 + diag, 1] += value
                for fend2 in range(((fend1 + startfend) / 2) * 2 + 4 - startfend, num_fends):
                    map2 = mapping[fend2]
                    if map2 == -1 or mids[fend2] - mids[fend1] > maxdistance or (diag == 0 and map2 == map1):
                        continue
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections is None:
                        value *= corrections[fend1] * corrections[fend2]
                    # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                    if not binning_corrections is None:
                        for j in range(num_parameters):
                            afend1 = fend1 + startfend
                            afend2 = fend2 + startfend
                            if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                                value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                            else:
                                value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                    # if finding distance, enrichment, or expected, correct for distance
                    if not parameters is None:
                        distance = log(<double>(mids[fend2] - mids[fend1]))
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    if ranges is None:
                        signal[map1, map2 - map1 - 1 + diag, 1] += value
                    else:
                        start1 = ranges[fend1, 0]
                        start2 = ranges[fend2, 0]
                        stop1 = ranges[fend1, 1]
                        stop2 = ranges[fend2, 1]
                        if start1 < start2 + diag:
                            signal[start1, start2 - start1 - 1 + diag, 1] += value * overlap[fend1, 0] * overlap[fend2, 0]
                        if start1 < stop2 + diag:
                            signal[start1, stop2 - start1 - 1 + diag, 1] += value * overlap[fend1, 0] * overlap[fend2, 1]
                        for l in range(max(start1, start2) + 1, stop2):
                            signal[start1, l - start1 - 1 + diag, 1] += value * overlap[fend1, 0]
                        if stop1 < start2 + diag:
                            signal[stop1, start2 - stop1 - 1 + diag, 1] += value * overlap[fend1, 1] * overlap[fend2, 0]
                        if stop1 < stop2 + diag:
                            signal[stop1, stop2 - stop1 - 1 + diag, 1] += value * overlap[fend1, 1] * overlap[fend2, 1]
                        for l in range(max(stop1, start2) + 1, stop2):
                            signal[stop1, l - stop1 - 1 + diag, 1] += value * overlap[fend1, 1]
                        for l in range(start1 + 1, stop1):
                            if l < start2 + diag:
                                signal[l, start2 - l - 1 + diag, 1] += value * overlap[fend2, 0]
                            if l < stop2 + diag:
                                signal[l, stop2 - l - 1 + diag, 1] += value * overlap[fend2, 1]
                            for m in range(max(l, start2) + 1, stop2):
                                signal[l, m - l - 1 + diag, 1] += value
        else:
            for j in range(num_bins - 1 + diag):
                for k in range(j + 1 - diag, min(num_bins, j + max_bin + 1 - diag)):
                    signal[j, k - j - 1 + diag, 1] += correction_sums[j] * correction_sums[k]
            if diag == 1:
                for j in range(num_bins):
                    signal[j, 0, 1] /= 2.0
            for fend1 in range(num_fends - 1):
                map1 = mapping[fend1]
                if map1 == -1:
                    continue
                map2 = mapping[fend1 + 1]
                if diag == 1 or map2 > map1:
                    signal[map1, map2 - map1 - 1 + diag, 1] -= corrections[fend1] * corrections[fend1 + 1]
                if (fend1 + startfend) % 2 == 0 and fend1 + 3 < num_fends:
                    map2 = mapping[fend1 + 3]
                    if diag == 1 or map2 > map1:
                        signal[map1, map2 - map1 - 1 + diag, 1] -= corrections[fend1] * corrections[fend1 + 3]
            if diag == 1:
                for fend1 in range(num_fends):
                    map1 = mapping[fend1]
                    if map1 != -1:
                        signal[map1, 0, 1] -= corrections[fend1] * corrections[fend1] / 2.0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_upper_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        np.ndarray[DTYPE_int_t, ndim=2] ranges,
        np.ndarray[DTYPE_t, ndim=2] overlap,
        double chrom_mean,
        int startfend,
        int maxdistance,
        int diag):
    cdef long long int fend1, fend2, afend1, afend2, j, k, index, map1, map2, index2, num_parameters
    cdef long long int start1, start2, stop1, stop2, l, m, index1
    cdef double distance, value
    cdef long long int num_fends = mapping.shape[0]
    cdef int diag2 = diag * 2
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5)) - diag
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    else:
        num_parameters = 0
    with nogil:
        if correction_sums is None or not ranges is None:
            for fend1 in range(num_fends - 2):
                map1 = mapping[fend1]
                if map1 == -1:
                    continue
                k = 0
                index = map1 * (num_bins - 1) - map1 * (map1 + 1 - diag2) / 2 - 1 + diag
                # find opposite strand adjacents, skipping same fragment and same strand adjacents
                fend2 = fend1 + 2
                map2 = mapping[fend2]
                if map2 >= 0 and (diag == 1 or map2 != map1) and mids[fend2] - mids[fend1] <= maxdistance:
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections is None:
                        value *= corrections[fend1] * corrections[fend2]
                    # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                    if not binning_corrections is None:
                        for j in range(num_parameters):
                            afend1 = fend1 + startfend
                            afend2 = fend2 + startfend
                            if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                                value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                            else:
                                value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                    # if finding distance, enrichment, or expected, correct for distance
                    if not parameters is None:
                        distance = log(<double>(mids[fend2] - mids[fend1]))
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    if ranges is None:
                        signal[index + map2, 1] += value
                    else:
                        start1 = ranges[fend1, 0]
                        start2 = ranges[fend2, 0]
                        stop1 = ranges[fend1, 1]
                        stop2 = ranges[fend2, 1]
                        index1 = start1 * (num_bins - 1) - start1 * (start1 + 1 - diag2) / 2 - 1 + diag
                        if start1 < start2 + diag:
                            signal[index1 + start2, 1] += value * overlap[fend1, 0] * overlap[fend2, 0]
                        if start1 < stop2 + diag:
                            signal[index1 + stop2, 1] += value * overlap[fend1, 0] * overlap[fend2, 1]
                        for l in range(max(start1, start2) + 1, stop2):
                            signal[index1 + l, 1] += value * overlap[fend1, 0]
                        index1 = stop1 * (num_bins - 1) - stop1 * (stop1 + 1) / 2 - 1
                        if stop1 < start2 + diag:
                            signal[index1 + start2, 1] += value * overlap[fend1, 1] * overlap[fend2, 0]
                        if stop1 < stop2 + diag:
                            signal[index1 + stop2, 1] += value * overlap[fend1, 1] * overlap[fend2, 1]
                        for l in range(max(stop1, start2) + 1, stop2):
                            signal[index1 + l, 1] += value * overlap[fend1, 1]
                        for l in range(start1 + 1, stop1):
                            index1 = l * (num_bins - 1) - l * (l + 1) / 2 - 1
                            if l < start2 + diag:
                                signal[index1 + start2, 1] += value * overlap[fend2, 0]
                            if l < stop2 + diag:
                                signal[index1 + stop2, 1] += value * overlap[fend2, 1]
                            for m in range(max(l, start2) + 1, stop2):
                                signal[index1 + m, 1] += value
                for fend2 in range(((fend1 + startfend) / 2) * 2 + 4 - startfend, num_fends):
                    map2 = mapping[fend2]
                    if map2 == -1 or mids[fend2] - mids[fend1] > maxdistance or (diag == 0 and map2 == map1):
                        continue
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections is None:
                        value *= corrections[fend1] * corrections[fend2]
                    # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                    if not binning_corrections is None:
                        for j in range(num_parameters):
                            afend1 = fend1 + startfend
                            afend2 = fend2 + startfend
                            if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                                value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                            else:
                                value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                    # if finding distance, enrichment, or expected, correct for distance
                    if not parameters is None:
                        distance = log(<double>(mids[fend2] - mids[fend1]))
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    if ranges is None:
                        signal[index + map2, 1] += value
                    else:
                        start1 = ranges[fend1, 0]
                        start2 = ranges[fend2, 0]
                        stop1 = ranges[fend1, 1]
                        stop2 = ranges[fend2, 1]
                        index1 = start1 * (num_bins - 1) - start1 * (start1 + 1) / 2 - 1
                        if start1 < start2 + diag:
                            signal[index1 + start2, 1] += value * overlap[fend1, 0] * overlap[fend2, 0]
                        if start1 < stop2 + diag:
                            signal[index1 + stop2, 1] += value * overlap[fend1, 0] * overlap[fend2, 1]
                        for l in range(max(start1, start2) + 1, stop2):
                            signal[index1 + l, 1] += value * overlap[fend1, 0]
                        index1 = stop1 * (num_bins - 1) - stop1 * (stop1 + 1) / 2 - 1
                        if stop1 < start2 + diag:
                            signal[index1 + start2, 1] += value * overlap[fend1, 1] * overlap[fend2, 0]
                        if stop1 < stop2 + diag:
                            signal[index1 + stop2, 1] += value * overlap[fend1, 1] * overlap[fend2, 1]
                        for l in range(max(stop1, start2) + 1, stop2):
                            signal[index1 + l, 1] += value * overlap[fend1, 1]
                        for l in range(start1 + 1, stop1):
                            index1 = l * (num_bins - 1) - l * (l + 1) / 2 - 1
                            if l < start2 + diag:
                                signal[index1 + start2, 1] += value * overlap[fend2, 0]
                            if l < stop2 + diag:
                                signal[index1 + stop2, 1] += value * overlap[fend2, 1]
                            for m in range(max(l, start2) + 1, stop2):
                                signal[index1 + m, 1] += value
        else:
            for j in range(num_bins - 1 + diag):
                index = j * (num_bins - 1) - (j * (j + 1 - diag2)) / 2 - 1 + diag
                for k in range(j + 1 - diag, num_bins):
                    signal[index + k, 1] += correction_sums[j] * correction_sums[k]
            if diag == 1:
                for j in range(num_bins):
                    index = j * num_bins - (j * (j - 1)) / 2
                    signal[index, 1] /= 2.0
            for fend1 in range(num_fends - 1):
                map1 = mapping[fend1]
                if map1 == -1:
                    continue
                index = map1 * (num_bins - 1) - (map1 * (map1 + 1 - diag2)) / 2 - 1 + diag
                map2 = mapping[fend1 + 1]
                if map2 >= 0 and (diag == 1 or map2 > map1):
                    signal[index + map2, 1] -= corrections[fend1] * corrections[fend1 + 1]
                if (fend1 + startfend) % 2 == 0 and fend1 + 3 < num_fends:
                    map2 = mapping[fend1 + 3]
                    if map2 >= 0 and (diag == 1 or map2 > map1):
                        signal[index + map2, 1] -= corrections[fend1] * corrections[fend1 + 3]
            if diag == 1:
                for fend1 in range(num_fends):
                    map1 = mapping[fend1]
                    if map1 != -1:
                        index = map1 * num_bins - map1 * (map1 - 1) / 2
                        signal[index, 1] -= corrections[fend1] * corrections[fend1] / 2.0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_binned_cis_compact_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        double chrom_mean,
        int startfend,
        int maxdistance):
    cdef long long int fend1, fend2, j, k, map1, map2, num_bins, max_bin
    cdef double distance, value
    cdef long long int num_fends = mapping.shape[0]
    num_bins = signal.shape[0]
    max_bin = signal.shape[1]
    with nogil:
        if correction_sums is None:
            for fend1 in range(num_fends):
                map1 = mapping[fend1]
                if map1 == -1:
                    continue
                k = 0
                for fend2 in range(fend1, num_fends):
                    map2 = mapping[fend2]
                    if map2 == -1 or mids[fend2] - mids[fend1] > maxdistance:
                        continue
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections is None:
                        value *= corrections[fend1] * corrections[fend2]
                    # if finding distance, enrichment, or expected, correct for distance
                    if not parameters is None:
                        distance = log(<double>(max(1, mids[fend2] - mids[fend1])))
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    signal[map1, map2 - map1, 1] += value
                signal[map1, 0, 1] /= 2
        else:
            for j in range(num_bins):
                for k in range(j, min(num_bins, j + max_bin)):
                    signal[j, k - j, 1] += correction_sums[j] * correction_sums[k]
                signal[j, 0, 1] /= 2.0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_binned_cis_upper_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        double chrom_mean,
        int startfend,
        int maxdistance):
    cdef long long int fend1, fend2, j, k, map1, map2, index, num_bins, max_bin
    cdef double distance, value
    cdef long long int num_fends = mapping.shape[0]
    num_bins = int(-0.5 + pow(0.25 + 2 * signal.shape[0], 0.5))
    with nogil:
        if correction_sums is None:
            for fend1 in range(num_fends):
                map1 = mapping[fend1]
                if map1 == -1:
                    continue
                index = map1 * (num_bins - 1) - map1 * (map1 - 1) / 2
                k = 0
                for fend2 in range(fend1, num_fends):
                    map2 = mapping[fend2]
                    if map2 == -1 or mids[fend2] - mids[fend1] > maxdistance:
                        continue
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections is None:
                        value *= corrections[fend1] * corrections[fend2]
                    # if finding distance, enrichment, or expected, correct for distance
                    if not parameters is None:
                        distance = log(<double>(max(1, mids[fend2] - mids[fend1])))
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    signal[index + map2, 1] += value
            for j in range(num_bins):
                signal[j * num_bins - j * (j - 1) / 2, 1] /= 2.0
        else:
            for j in range(num_bins):
                index = j * (num_bins - 1) - j * (j - 1) / 2
                for k in range(j, num_bins):
                    signal[index + k, 1] = correction_sums[j] * correction_sums[k]
                signal[index + j, 1] /= 2.0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_subregion_observed(
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] indices not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        int startfend1,
        int startfend2):
    cdef long long int fend1, fend2, i, j, map1, map2, map1b, map2b, valid1, valid2
    cdef long long int num_fends = indices.shape[0] - 1
    cdef long long int startfend = min(startfend1, startfend2)
    cdef long long int stopfend1 = startfend1 + mapping1.shape[0]
    cdef long long int stopfend2 = startfend2 + mapping2.shape[0]
    with nogil:
        for i in range(num_fends):
            valid1 = 0
            valid2 = 0
            fend1 = i + startfend
            if fend1 >= startfend1 and fend1 < stopfend1:
                map1 = mapping1[fend1 - startfend1]
                if map1 >= 0:
                    valid1 = 1
            if fend1 >= startfend2 and fend1 < stopfend2:
                map2 = mapping2[fend1 - startfend2]
                if map2 >= 0:
                    valid2 = 1
            if valid1 == 0 and valid2 == 0:
                continue
            for j in range(indices[i], indices[i + 1]):
                fend2 = data[j, 1]
                if valid1 == 1 and fend2 >= startfend2 and fend2 < stopfend2:
                    map2b = mapping2[fend2 - startfend2]
                    if map2b >= 0:
                        signal[map1, map2b, 0] += data[j, 2]
                if valid2 == 1 and fend2 >= startfend1 and fend2 < stopfend1:
                    map1b = mapping1[fend2 - startfend1]
                    if map1b >= 0:
                        signal[map1b, map2, 0] += data[j, 2]
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_subregion_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=1] corrections1,
        np.ndarray[DTYPE_t, ndim=1] corrections2,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        np.ndarray[DTYPE_64_t, ndim=1] correction_sums1,
        np.ndarray[DTYPE_64_t, ndim=1] correction_sums2,
        np.ndarray[DTYPE_int_t, ndim=1] mids1,
        np.ndarray[DTYPE_int_t, ndim=1] mids2,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        double chrom_mean,
        int startfend1,
        int startfend2):
    cdef long long int fend1, fend2, afend1, afend2, i, j, k, map1, map2, index, num_parameters, num_bins1, num_bins2
    cdef long long int start1, start2, stop1, stop2, l, m
    cdef double distance, value
    cdef long long int num_fends1 = mapping1.shape[0]
    cdef long long int num_fends2 = mapping2.shape[0]
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    else:
        num_parameters = 0
    if not correction_sums1 is None:
        num_bins1 = correction_sums1.shape[0]
        num_bins2 = correction_sums2.shape[0]
    with nogil:
        if correction_sums1 is None:
            for fend1 in range(num_fends1):
                map1 = mapping1[fend1]
                if map1 == -1:
                    continue
                afend1 = fend1 + startfend1
                for fend2 in range(num_fends2):
                    map2 = mapping2[fend2]
                    if map2 == -1:
                        continue
                    afend2 = fend2 + startfend2
                    diff = afend1 - afend2
                    if (diff < 2 and diff > -2) or (afend1 % 2 == 0 and diff == 3) or (afend2 % 2 == 0 and diff == -3):
                        continue
                     # give starting expected value
                    value = 1.0
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections1 is None:
                        value *= corrections1[fend1] * corrections2[fend2]
                    # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                    if not binning_corrections is None:
                        for j in range(num_parameters):
                            if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                                value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                            else:
                                value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                    # if finding distance, enrichment, or expected, correct for distance
                    if not parameters is None:
                        distance = mids2[fend2] - mids1[fend1]
                        if distance > 0:
                            distance = log(<double>distance)
                        else:
                            distance = log(<double>(-distance))
                        k = 0
                        while distance > parameters[k, 0]:
                            k += 1
                        value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                    signal[map1, map2, 1] += value
        else:
            for i in range(num_bins1):
                for j in range(num_bins2):
                    signal[i, j, 1] += correction_sums1[i] * correction_sums2[j]
            for i in range(num_fends1):
                map1 = mapping1[i]
                if map1 >= 0:
                    fend2 = i + startfend1 - startfend2
                    if fend2 >= 0 and fend2 < num_fends2:
                        map2 = mapping2[fend2]
                        if map2 >= 0:
                            signal[map1, map2, 1] -= corrections1[i] * corrections2[fend2]
                    fend2 = i + startfend1 - startfend2 - 1
                    if fend2 >= 0 and fend2 < num_fends2:
                        map2 = mapping2[fend2]
                        if map2 >= 0:
                            signal[map1, map2, 1] -= corrections1[i] * corrections2[fend2]
                    fend2 = i + startfend1 - startfend2 + 1
                    if fend2 >= 0 and fend2 < num_fends2:
                        map2 = mapping2[fend2]
                        if map2 >= 0:
                            signal[map1, map2, 1] -= corrections1[i] * corrections2[fend2]
                    if (i + startfend1) % 2 == 0:
                        fend2 = i + startfend1 - startfend2 + 3
                        if fend2 >= 0 and fend2 < num_fends2:
                            map2 = mapping2[fend2]
                            if map2 >= 0:
                                signal[map1, map2, 1] -= corrections1[i] * corrections2[fend2]
                    else:
                        fend2 = i + startfend1 - startfend2 - 3
                        if fend2 >= 0 and fend2 < num_fends2:
                            map2 = mapping2[fend2]
                            if map2 >= 0:
                                signal[map1, map2, 1] -= corrections1[i] * corrections2[fend2]
    return None



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_compact_to_compact(
        np.ndarray[DTYPE_t, ndim=3] binned,
        np.ndarray[DTYPE_t, ndim=3] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        int diag):
    cdef long long int i, j, fend2
    cdef long long int num_bins = binned.shape[0]
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    with nogil:
        for i in range(num_fends - 1 + diag):
            if mapping[i] == -1:
                continue
            for j in range(min(num_fends - i - 1 + diag, max_fend)):
                fend2 = j + i + 1 - diag
                if mapping[fend2] == -1 or mapping[fend2] == mapping[i]:
                    continue
                else:
                    binned[mapping[i], mapping[fend2] - mapping[i] - 1 + diag, 0] += unbinned[i, j, 0]
                    binned[mapping[i], mapping[fend2] - mapping[i] - 1 + diag, 1] += unbinned[i, j, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_compact_to_upper(
        np.ndarray[DTYPE_t, ndim=2] binned,
        np.ndarray[DTYPE_t, ndim=3] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        int num_bins,
        int diag):
    cdef long long int i, j, index
    cdef long long int diag2 = diag * 2
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    with nogil:
        for i in range(num_fends - 1 + diag):
            if mapping[i] == -1:
                continue
            index = mapping[i] * (num_bins - 1) - mapping[i] * (mapping[i] + 1 - diag2) / 2 - 1 + diag
            for j in range(i + 1 - diag, num_fends):
                if mapping[j] == -1 or mapping[j] == mapping[i]:
                    continue
                else:
                    binned[index + mapping[j], 0] += unbinned[i, j - i - 1 + diag, 0]
                    binned[index + mapping[j], 1] += unbinned[i, j - i - 1 + diag, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_upper_to_compact(
        np.ndarray[DTYPE_t, ndim=3] binned,
        np.ndarray[DTYPE_t, ndim=2] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        int diag):
    cdef long long int i, j, index
    cdef long long int diag2 = diag * 2
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        for i in range(num_fends - 1 + diag):
            if mapping[i] == -1:
                continue
            index = i * (num_fends - 1) - i * (i + 1 - diag2) / 2 - 1 + diag
            for j in range(i + 1 - diag, num_fends):
                if mapping[j] == -1 or mapping[j] == mapping[i]:
                    continue
                else:
                    binned[mapping[i], mapping[j] - mapping[i] - 1 + diag, 0] += unbinned[index + j, 0]
                    binned[mapping[i], mapping[j] - mapping[i] - 1 + diag, 1] += unbinned[index + j, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_upper_to_upper(
        np.ndarray[DTYPE_t, ndim=2] binned,
        np.ndarray[DTYPE_t, ndim=2] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        int num_bins,
        int diag):
    cdef long long int i, j, index, index2
    cdef long long int diag2 = diag * 2
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        for i in range(num_fends - 1 + diag):
            if mapping[i] == -1:
                continue
            index = i * (num_fends - 1) - i * (i + 1 - diag2) / 2 - 1 + diag
            index2 = mapping[i] * (num_bins - 1) - mapping[i] * (mapping[i] + 1 - diag2) / 2 - 1 + diag
            for j in range(i + 1 - diag, num_fends):
                if mapping[j] == -1 or mapping[j] == mapping[i]:
                    continue
                else:
                    binned[index2 + mapping[j], 0] += unbinned[index + j, 0]
                    binned[index2 + mapping[j], 1] += unbinned[index + j, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin_upper_from_upper(
        np.ndarray[DTYPE_t, ndim=2] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] ub_mids not None,
        np.ndarray[DTYPE_t, ndim=2] binned not None,
        np.ndarray[DTYPE_int_t, ndim=2] bounds not None,
        np.ndarray[DTYPE_int_t, ndim=1] b_mids not None,
        int minobservations,
        int maxsearch,
        int removefailed,
        int skipinvalid,
        int diag):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist, index, index2
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int num_fends = ub_mids.shape[0]
    cdef long long int diag2 = diag * 2
    with nogil:
        for x in range(num_bins - 1 + diag):
            index = x * (num_bins - 1) - x * (x + 1 - diag2) / 2 - 1 + diag
            for y in range(x + 1 - diag, num_bins):
                if skipinvalid == 1 and binned[index + y, 1] == 0.0:
                    continue
                # if bin already meets our criteria, skip
                if binned[index + y, 0] >= minobservations:
                    continue
                # otherwise, set boarding unbinned positions according to bounds
                lX = bounds[x, 0]
                uX = bounds[x, 1] - 1
                lY = bounds[y, 0]
                uY = bounds[y, 1] - 1
                # find distance in each direction
                if lX > 0:
                    lX_dist = b_mids[x] - ub_mids[lX - 1]
                else:
                    lX_dist = 1000000000
                if uX < num_fends - 1:
                    uX_dist = ub_mids[uX + 1] - b_mids[x]
                else:
                    uX_dist = 1000000000
                if lY > 0:
                    lY_dist = b_mids[y] - ub_mids[lY - 1]
                else:
                    lY_dist = 1000000000
                if uY < num_fends - 1:
                    uY_dist = ub_mids[uY + 1] - b_mids[y]
                else:
                    uY_dist = 1000000000
                # determine min distance
                min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                # keep searching while less than maxsearch and minobservations
                while (maxsearch == 0 or min_dist < maxsearch) and binned[index + y, 0] < minobservations:
                    # find min dist, update distance, and add new row or col of observations
                    if min_dist == lX_dist:
                        lX -= 1
                        if lX > 0:
                            lX_dist = b_mids[x] - ub_mids[lX - 1]
                        else:
                            lX_dist = 1000000000
                        index2 = lX * (num_fends - 1) - lX * (lX + 1 - diag2) / 2 - 1 + diag
                        for i in range(max(lX + 1 - diag, lY), uY + 1):
                            binned[index + y, 0] += unbinned[index2 + i, 0]
                            binned[index + y, 1] += unbinned[index2 + i, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        index2 = uX * (num_fends - 1) - uX * (uX + 1 - diag2) / 2 - 1 + diag
                        for i in range(max(uX + 1 - diag, lY), uY + 1):
                            binned[index + y, 0] += unbinned[index2 + i, 0]
                            binned[index + y, 1] += unbinned[index2 + i, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(lX, min(uX + 1, lY + diag)):
                            index2 = i * (num_fends - 1) - i * (i + 1 - diag2) / 2 - 1 + diag
                            binned[index + y, 0] += unbinned[index2 + lY, 0]
                            binned[index + y, 1] += unbinned[index2 + lY, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(lX, min(uX + 1, uY + diag)):
                            index2 = i * (num_fends - 1) - i * (i + 1 - diag2) / 2 - 1 + diag
                            binned[index + y, 0] += unbinned[index2 + uY, 0]
                            binned[index + y, 1] += unbinned[index2 + uY, 1]
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                if binned[index + y, 0] < minobservations and removefailed == 1:
                    binned[index + y, 0] = 0
                    binned[index + y, 1] = 0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin_compact_from_upper(
        np.ndarray[DTYPE_t, ndim=2] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] ub_mids not None,
        np.ndarray[DTYPE_t, ndim=3] binned not None,
        np.ndarray[DTYPE_int_t, ndim=2] bounds not None,
        np.ndarray[DTYPE_int_t, ndim=1] b_mids not None,
        int minobservations,
        int maxsearch,
        int removefailed,
        int skipinvalid,
        int diag):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist, index
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int max_bin = binned.shape[1]
    cdef long long int num_fends = ub_mids.shape[0]
    cdef long long int diag2 = diag * 2
    with nogil:
        for x in range(num_bins - 1 + diag):
            for y in range(x + 1 - diag, min(x + max_bin + 1 - diag, num_bins)):
                if skipinvalid == 1 and binned[x, y - x - 1 + diag, 1] == 0.0:
                    continue
                # if bin already meets our criteria, skip
                if binned[x, y - x - 1 + diag, 0] >= minobservations:
                    continue
                # otherwise, set boarding unbinned positions according to bounds
                lX = bounds[x, 0]
                uX = bounds[x, 1] - 1
                lY = bounds[y, 0]
                uY = bounds[y, 1] - 1
                # find distance in each direction
                if lX > 0:
                    lX_dist = b_mids[x] - ub_mids[lX - 1]
                else:
                    lX_dist = 1000000000
                if uX < num_fends - 1:
                    uX_dist = ub_mids[uX + 1] - b_mids[x]
                else:
                    uX_dist = 1000000000
                if lY > 0:
                    lY_dist = b_mids[y] - ub_mids[lY - 1]
                else:
                    lY_dist = 1000000000
                if uY < num_fends - 1:
                    uY_dist = ub_mids[uY + 1] - b_mids[y]
                else:
                    uY_dist = 1000000000
                # determine min distance
                min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                # keep searching while less than maxsearch and minobservations
                while (maxsearch == 0 or min_dist < maxsearch) and binned[x, y - x - 1 + diag, 0] < minobservations:
                    # find min dist, update distance, and add new row or col of observations
                    if min_dist == lX_dist:
                        lX -= 1
                        if lX > 0:
                            lX_dist = b_mids[x] - ub_mids[lX - 1]
                        else:
                            lX_dist = 1000000000
                        index = lX * (num_fends - 1) - lX * (lX + 1 - diag2) / 2 - 1 + diag
                        for i in range(max(lX + 1 - diag, lY), uY + 1):
                            binned[x, y - x - 1 + diag, 0] += unbinned[index + i, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[index + i, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        index = uX * (num_fends - 1) - uX * (uX + 1 - diag2) / 2 - 1 + diag
                        for i in range(max(uX + 1 - diag, lY), uY + 1):
                            binned[x, y - x - 1 + diag, 0] += unbinned[index + i, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[index + i, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(lX, min(uX + 1, lY + diag)):
                            index = i * (num_fends - 1) - i * (i + 1 - diag2) / 2 - 1 + diag
                            binned[x, y - x - 1 + diag, 0] += unbinned[index + lY, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[index + lY, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(lX, min(uX + 1, uY + diag)):
                            index = i * (num_fends - 1) - i * (i + 1 - diag2) / 2 - 1 + diag
                            binned[x, y - x - 1 + diag, 0] += unbinned[index + uY, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[index + uY, 1]
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                if binned[x, y - x - 1 + diag, 0] < minobservations and removefailed == 1:
                    binned[x, y - x - 1 + diag, 0] = 0
                    binned[x, y - x - 1 + diag, 1] = 0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin_upper_from_compact(
        np.ndarray[DTYPE_t, ndim=3] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] ub_mids not None,
        np.ndarray[DTYPE_t, ndim=2] binned not None,
        np.ndarray[DTYPE_int_t, ndim=2] bounds not None,
        np.ndarray[DTYPE_int_t, ndim=1] b_mids not None,
        int minobservations,
        int maxsearch,
        int removefailed,
        int skipinvalid,
        int diag):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist, index
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int num_fends = ub_mids.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    cdef long long int diag2 = diag * 2
    with nogil:
        for x in range(num_bins - 1 + diag):
            index = x * (num_bins - 1) - x * (x + 1 - diag2) / 2 - 1 + diag
            for y in range(x + 1 - diag, num_bins):
                if skipinvalid == 1 and binned[index + y, 1] == 0.0:
                    continue
                # if bin already meets our criteria, skip
                if binned[index + y, 0] >= minobservations:
                    continue
                # otherwise, set boarding unbinned positions according to bounds
                lX = bounds[x, 0]
                uX = bounds[x, 1] - 1
                lY = bounds[y, 0]
                uY = bounds[y, 1] - 1
                # find distance in each direction
                if lX > 0:
                    lX_dist = b_mids[x] - ub_mids[lX - 1]
                else:
                    lX_dist = 1000000000
                if uX < num_fends - 1:
                    uX_dist = ub_mids[uX + 1] - b_mids[x]
                else:
                    uX_dist = 1000000000
                if lY > 0:
                    lY_dist = b_mids[y] - ub_mids[lY - 1]
                else:
                    lY_dist = 1000000000
                if uY < num_fends - 1:
                    uY_dist = ub_mids[uY + 1] - b_mids[y]
                else:
                    uY_dist = 1000000000
                # determine min distance
                min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                # keep searching while less than maxsearch and minobservations
                while (maxsearch == 0 or min_dist < maxsearch) and binned[index + y, 0] < minobservations:
                    # find min dist, update distance, and add new row or col of observations
                    if min_dist == lX_dist:
                        lX -= 1
                        if lX > 0:
                            lX_dist = b_mids[x] - ub_mids[lX - 1]
                        else:
                            lX_dist = 1000000000
                        for i in range(max(lX + 1 - diag, lY), min(uY + 1, lX + max_fend + 1 - diag)):
                            binned[index + y, 0] += unbinned[lX, i - lX - 1 + diag, 0]
                            binned[index + y, 1] += unbinned[lX, i - lX - 1 + diag, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        for i in range(max(uX + 1 - diag, lY), min(uY + 1, uX + max_fend + 1 - diag)):
                            binned[index + y, 0] += unbinned[uX, i - uX - 1 + diag, 0]
                            binned[index + y, 1] += unbinned[uX, i - uX - 1 + diag, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(max(lX, lY - max_fend - 1 + diag), min(uX + 1, lY + diag)):
                            binned[index + y, 0] += unbinned[i, lY - i - 1 + diag, 0]
                            binned[index + y, 1] += unbinned[i, lY - i - 1 + diag, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(max(lX, uY - max_fend - 1 + diag), min(uX + 1, uY + diag)):
                            binned[index + y, 0] += unbinned[i, uY - i - 1 + diag, 0]
                            binned[index + y, 1] += unbinned[i, uY - i - 1 + diag, 1]
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                if binned[index + y, 0] < minobservations and removefailed == 1:
                    binned[index + y, 0] = 0
                    binned[index + y, 1] = 0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin_compact_from_compact(
        np.ndarray[DTYPE_t, ndim=3] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] ub_mids not None,
        np.ndarray[DTYPE_t, ndim=3] binned not None,
        np.ndarray[DTYPE_int_t, ndim=2] bounds not None,
        np.ndarray[DTYPE_int_t, ndim=1] b_mids not None,
        int minobservations,
        int maxsearch,
        int removefailed,
        int skipinvalid,
        int diag):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int max_bin = binned.shape[1]
    cdef long long int num_fends = ub_mids.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    with nogil:
        for x in range(num_bins - 1 + diag):
            for y in range(x + 1 - diag, min(x + max_bin + 1 - diag, num_bins)):
                if skipinvalid == 1 and binned[x, y - x - 1 + diag, 1] == 0.0:
                    continue
                # if bin already meets our criteria, skip
                if binned[x, y - x - 1 + diag, 0] >= minobservations:
                    continue
                # otherwise, set boarding unbinned positions according to bounds
                lX = bounds[x, 0]
                uX = bounds[x, 1] - 1
                lY = bounds[y, 0]
                uY = bounds[y, 1] - 1
                # find distance in each direction
                if lX > 0:
                    lX_dist = b_mids[x] - ub_mids[lX - 1]
                else:
                    lX_dist = 1000000000
                if uX < num_fends - 1:
                    uX_dist = ub_mids[uX + 1] - b_mids[x]
                else:
                    uX_dist = 1000000000
                if lY > 0:
                    lY_dist = b_mids[y] - ub_mids[lY - 1]
                else:
                    lY_dist = 1000000000
                if uY < num_fends - 1:
                    uY_dist = ub_mids[uY + 1] - b_mids[y]
                else:
                    uY_dist = 1000000000
                # determine min distance
                min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                # keep searching while less than maxsearch and minobservations
                while (maxsearch == 0 or min_dist < maxsearch) and binned[x, y - x - 1 + diag, 0] < minobservations:
                    # find min dist, update distance, and add new row or col of observations
                    if min_dist == lX_dist:
                        lX -= 1
                        if lX > 0:
                            lX_dist = b_mids[x] - ub_mids[lX - 1]
                        else:
                            lX_dist = 1000000000
                        for i in range(max(lX + 1 - diag, lY), min(uY + 1, lX + max_fend + 1 - diag)):
                            binned[x, y - x - 1 + diag, 0] += unbinned[lX, i - lX - 1 + diag, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[lX, i - lX - 1 + diag, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        for i in range(max(uX + 1 - diag, lY), min(uY + 1, uX + max_fend + 1 - diag)):
                            binned[x, y - x - 1 + diag, 0] += unbinned[uX, i - uX - 1 + diag, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[uX, i - uX - 1 + diag, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(max(lX, lY - max_fend - 1 + diag), min(uX + 1, lY + diag)):
                            binned[x, y - x - 1 + diag, 0] += unbinned[i, lY - i - 1 + diag, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[i, lY - i - 1 + diag, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(max(lX, uY - max_fend - 1 + diag), min(uX + 1, uY + diag)):
                            binned[x, y - x - 1 + diag, 0] += unbinned[i, uY - i - 1 + diag, 0]
                            binned[x, y - x - 1 + diag, 1] += unbinned[i, uY - i - 1 + diag, 1]
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                if binned[x, y - x - 1 + diag, 0] < minobservations and removefailed == 1:
                    binned[x, y - x - 1 + diag, 0] = 0
                    binned[x, y - x - 1 + diag, 1] = 0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin_unbinned_upper(
        np.ndarray[DTYPE_t, ndim=2] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=2] binned not None,
        int minobservations):
    cdef long long int x, y, i, j, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, 
    cdef long long int old_lX, old_uX, old_lY, old_uY, max_dist, min_dist, index, index2
    cdef double temp_count, temp_exp
    cdef long long int num_fends = mids.shape[0]
    with nogil:
        for x in range(num_fends - 1):
            index = x * num_fends - x * (x + 1) / 2 - x - 1
            lX = x
            uX = x
            lY = x
            uY = x
            for y in range(x + 1, num_fends):
                # if bin already meets our criteria, skip
                if binned[index + y, 0] >= minobservations:
                    lX = x
                    uX = x
                    lY = y
                    uY = y
                    continue
                # otherwise shift lower y bound forward 1
                elif y == x + 1:
                    lY = y
                    uY = y
                else:
                    old_lX = lX
                    old_uX = uX
                    old_lY = lY
                    old_uY = uY
                    lY += 1
                    uY = max(uY, y)
                    # determine if any bounds are too far away
                    lX_dist = mids[x] - mids[lX]
                    uX_dist = mids[uX] - mids[x]
                    lY_dist = mids[y] - mids[lY]
                    uY_dist = mids[uY] - mids[y]
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                    if min_dist == 0:
                        lX = x
                        uX = x
                        lY = y
                        uY = y
                    else:
                        while lX < x and mids[x] - mids[lX] > min_dist:
                            lX += 1
                        while uX > x and mids[uX] - mids[x] > min_dist:
                            uX -= 1
                        while lY < y and mids[y] - mids[lY] > min_dist:
                            lY += 1
                        while uY > y and mids[uY] - mids[y] > min_dist:
                            uY -= 1
                        # check to see total enclosed in rectangle
                        binned[index + y, 0] = binned[index + y - 1, 0]
                        binned[index + y, 1] = binned[index + y - 1, 1]
                        if lX > old_lX:
                            for i in range(old_lX, lX):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                for j in range(max(i + 1, old_lY), old_uY + 1):
                                    binned[index + y, 0] -= unbinned[index2 + j, 0]
                                    binned[index + y, 1] -= unbinned[index2 + j, 1]
                        if uX < old_uX:
                            for i in range(uX + 1, old_uX + 1):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                for j in range(max(i + 1, old_lY), old_uY + 1):
                                    binned[index + y, 0] -= unbinned[index2 + j, 0]
                                    binned[index + y, 1] -= unbinned[index2 + j, 1]
                        if lY > old_lY:
                            for i in range(lX, uX + 1):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                for j in range(max(i + 1, old_lY), lY):
                                    binned[index + y, 0] -= unbinned[index2 + j, 0]
                                    binned[index + y, 1] -= unbinned[index2 + j, 1]
                        if uY < old_uY:
                            for i in range(lX, uX + 1):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                for j in range(max(i + 1, uY + 1), old_uY + 1):
                                    binned[index + y, 0] -= unbinned[index2 + j, 0]
                                    binned[index + y, 1] -= unbinned[index2 + j, 1]
                        elif old_uY < uY:
                            for i in range(lX, uX + 1):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                for j in range(max(i + 1, old_uY + 1), uY + 1):
                                    binned[index + y, 0] += unbinned[index2 + j, 0]
                                    binned[index + y, 1] += unbinned[index2 + j, 1]
                # check if we've met our criterion
                if binned[index + y, 0] < minobservations:
                    # find distance in each direction
                    if lX > 0:
                        lX_dist = mids[x] - mids[lX - 1]
                    else:
                        lX_dist = 1000000000
                    if uX < num_fends - 1:
                        uX_dist = mids[uX + 1] - mids[x]
                    else:
                        uX_dist = 1000000000
                    if lY > 0:
                        lY_dist = mids[y] - mids[lY - 1]
                    else:
                        lY_dist = 1000000000
                    if uY < num_fends - 1:
                        uY_dist = mids[uY + 1] - mids[y]
                    else:
                        uY_dist = 1000000000
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                    # keep searching while less minobservations
                    while binned[index + y, 0] < minobservations:
                        # find min dist, update distance, and add new row or col of observations
                        if min_dist == lX_dist:
                            lX -= 1
                            if lX > 0:
                                lX_dist = mids[x] - mids[lX - 1]
                            else:
                                lX_dist = 1000000000
                            index2 = lX * num_fends - lX * (lX + 1) / 2 - lX - 1
                            for i in range(max(lX + 1, lY), uY + 1):
                                binned[index + y, 0] += unbinned[index2 + i, 0]
                                binned[index + y, 1] += unbinned[index2 + i, 1]
                        elif min_dist == uX_dist:
                            uX += 1
                            if uX < num_fends - 1:
                                uX_dist = mids[uX + 1] - mids[x]
                            else:
                                uX_dist = 1000000000
                            index2 = uX * num_fends - uX * (uX + 1) / 2 - uX - 1
                            for i in range(max(uX + 1, lY), uY + 1):
                                binned[index + y, 0] += unbinned[index2 + i, 0]
                                binned[index + y, 1] += unbinned[index2 + i, 1]
                        elif min_dist == lY_dist:
                            lY -= 1
                            if lY > 0:
                                lY_dist = mids[y] - mids[lY - 1]
                            else:
                                lY_dist = 1000000000
                            for i in range(lX, min(uX + 1, lY)):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                binned[index + y, 0] += unbinned[index2 + lY, 0]
                                binned[index + y, 1] += unbinned[index2 + lY, 1]
                        elif min_dist == uY_dist:
                            uY += 1
                            if uY < num_fends - 1:
                                uY_dist = mids[uY + 1] - mids[y]
                            else:
                                uY_dist = 1000000000
                            for i in range(lX, min(uX + 1, uY)):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                binned[index + y, 0] += unbinned[index2 + uY, 0]
                                binned[index + y, 1] += unbinned[index2 + uY, 1]
                        # determine min distance
                        min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                # else if we've overshot our goal, remove as many rows/cols as needed
                elif binned[index + y, 0] > minobservations:
                    # find furthest row or col
                    lX_dist = mids[x] - mids[lX]
                    uX_dist = mids[uX] - mids[x]
                    lY_dist = mids[y] - mids[lY]
                    uY_dist = mids[uY] - mids[y]
                    max_dist = max(max(lX_dist, uX_dist), max(lY_dist, uY_dist))
                    # if bin only includes itself, continue
                    if max_dist == 0:
                        continue
                    continue_search = 1
                    while continue_search == 1:
                        continue_search = 0
                        # find sum of furthest col or row
                        temp_count = 0.0
                        temp_exp = 0.0
                        if lX_dist == max_dist:
                            index2 = lX * num_fends - lX * (lX + 1) / 2 - lX - 1
                            for i in range(max(lX + 1, lY), uY + 1):
                                temp_count += unbinned[index2 + i, 0]
                                temp_exp += unbinned[index2 + i, 1]
                            if binned[index + y, 0] - temp_count >= minobservations:
                                binned[index + y, 0] -= temp_count
                                binned[index + y, 1] -= temp_exp
                                lX += 1
                                if lX < x:
                                    lX_dist = mids[x] - mids[lX]
                                else:
                                    lX_dist = 0
                                continue_search = 1
                        elif uX_dist == max_dist:
                            index2 = uX * num_fends - uX * (uX + 1) / 2 - uX - 1
                            for i in range(max(uX + 1, lY), uY + 1):
                                temp_count += unbinned[index2 + i, 0]
                                temp_exp += unbinned[index2 + i, 1]
                            if binned[index + y, 0] - temp_count >= minobservations:
                                binned[index + y, 0] -= temp_count
                                binned[index + y, 1] -= temp_exp
                                uX -= 1
                                if uX > x:
                                    uX_dist = mids[uX] - mids[x]
                                else:
                                    uX_dist = 0
                                continue_search = 1
                        elif lY_dist == max_dist:
                            for i in range(lX, min(uX + 1, lY)):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                temp_count += unbinned[index2 + lY, 0]
                                temp_exp += unbinned[index2 + lY, 1]
                            if binned[index + y, 0] - temp_count >= minobservations:
                                binned[index + y, 0] -= temp_count
                                binned[index + y, 1] -= temp_exp
                                lY += 1
                                if lY < y:
                                    lY_dist = mids[y] - mids[lY]
                                else:
                                    lY_dist = 0
                                continue_search = 1
                        elif uY_dist == max_dist:
                            for i in range(lX, min(uX + 1, uY)):
                                index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                                temp_count += unbinned[index2 + uY, 0]
                                temp_exp += unbinned[index2 + uY, 1]
                            if binned[index + y, 0] - temp_count >= minobservations:
                                binned[index + y, 0] -= temp_count
                                binned[index + y, 1] -= temp_exp
                                uY -= 1
                                if uY > y:
                                    uY_dist = mids[uY] - mids[y]
                                else:
                                    uY_dist = 0
                                continue_search = 1
                        max_dist = max(max(lX_dist, uX_dist), max(lY_dist, uY_dist))
                        # if bin only includes itself, continue
                        if max_dist == 0:
                            continue_search = 0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_trans_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None):
    cdef long long int fend1, fend2, i, stop
    cdef long long int num_fends1 = mapping1.shape[0]
    cdef long long int num_fends2 = mapping2.shape[0]
    with nogil:
        for fend1 in range(num_fends1 - 1):
            if mapping1[fend1] == -1:
                continue
            i = indices[fend1]
            stop = indices[fend1 + 1]
            while i < stop and data[i, 1] < 0:
                i += 1 
            while i < stop and data[i, 1] < num_fends2:
                fend2 = data[i, 1]
                if mapping2[fend2] != -1:
                    signal[mapping1[fend1], mapping2[fend2], 0] += data[i, 2]
                i += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_trans_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=1] corrections1,
        np.ndarray[DTYPE_t, ndim=1] corrections2,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        np.ndarray[DTYPE_64_t, ndim=1] correction_sums1,
        np.ndarray[DTYPE_64_t, ndim=1] correction_sums2,
        double trans_mean,
        int startfend1,
        int startfend2):
    cdef long long int j, fend1, fend2, afend1, afend2, map1, map2, index, num_parameters, num_bins1, num_bins2
    cdef double value
    cdef long long int num_fends1 = mapping1.shape[0]
    cdef long long int num_fends2 = mapping2.shape[0]
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    else:
        num_parameters = 0
    if not correction_sums1 is None:
        num_bins1 = correction_sums1.shape[0]
        num_bins2 = correction_sums2.shape[0]
    with nogil:
        if correction_sums1 is None:
            for fend1 in range(num_fends1 - 1):
                map1 = mapping1[fend1]
                if map1 == -1:
                    continue
                for fend2 in range(num_fends2):
                    map2 = mapping2[fend2]
                    if map2 == -1:
                        continue
                     # give starting expected value
                    value = trans_mean
                    # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                    if not corrections1 is None:
                        value *= corrections1[fend1] * corrections2[fend2]
                    # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                    if not binning_corrections is None:
                        for j in range(num_parameters):
                            afend1 = fend1 + startfend1
                            afend2 = fend2 + startfend2
                            if fend_indices[afend1, j, 0] < fend_indices[afend2, j, 0]:
                                value *= binning_corrections[fend_indices[afend1, j, 1] + fend_indices[afend2, j, 0]]
                            else:
                                value *= binning_corrections[fend_indices[afend2, j, 1] + fend_indices[afend1, j, 0]]
                    signal[map1, map2, 1] += value
        else:
            for fend1 in range(num_bins1):
                for fend2 in range(num_bins2):
                    signal[fend1, fend2, 1] += trans_mean * correction_sums1[fend1] * correction_sums2[fend2]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin_trans(
        np.ndarray[DTYPE_t, ndim=3] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] ub_mids1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] ub_mids2 not None,
        np.ndarray[DTYPE_t, ndim=3] binned not None,
        np.ndarray[DTYPE_int_t, ndim=2] bounds1 not None,
        np.ndarray[DTYPE_int_t, ndim=2] bounds2 not None,
        np.ndarray[DTYPE_int_t, ndim=1] b_mids1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] b_mids2 not None,
        int minobservations,
        int maxsearch,
        int removefailed):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist
    cdef long long int num_bins1 = bounds1.shape[0]
    cdef long long int num_bins2 = bounds2.shape[0]
    cdef long long int num_fends1 = ub_mids1.shape[0]
    cdef long long int num_fends2 = ub_mids2.shape[0]
    with nogil:
        for x in range(num_bins1):
            for y in range(num_bins2):
                # if bin already meets our criteria, skip
                if binned[x, y, 0] >= minobservations:
                    continue
                # otherwise, set bordering unbinned positions according to bounds
                lX = bounds1[x, 0]
                uX = bounds1[x, 1] - 1
                lY = bounds2[y, 0]
                uY = bounds2[y, 1] - 1
                # find distance in each direction
                if lX > 0:
                    lX_dist = b_mids1[x] - ub_mids1[lX - 1]
                else:
                    lX_dist = 1000000000
                if uX < num_fends1 - 1:
                    uX_dist = ub_mids1[uX + 1] - b_mids1[x]
                else:
                    uX_dist = 1000000000
                if lY > 0:
                    lY_dist = b_mids2[y] - ub_mids2[lY - 1]
                else:
                    lY_dist = 1000000000
                if uY < num_fends2 - 1:
                    uY_dist = ub_mids2[uY + 1] - b_mids2[y]
                else:
                    uY_dist = 1000000000
                # determine min distance
                min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                # keep searching while less than maxsearch and minobservations
                while (maxsearch == 0 or min_dist < maxsearch) and binned[x, y, 0] < minobservations:
                    # find min dist, update distance, and add new row or col of observations
                    if min_dist == lX_dist:
                        lX -= 1
                        if lX > 0:
                            lX_dist = b_mids1[x] - ub_mids1[lX - 1]
                        else:
                            lX_dist = 1000000000
                        for i in range(lY, uY + 1):
                            binned[x, y, 0] += unbinned[lX, i, 0]
                            binned[x, y, 1] += unbinned[lX, i, 1]
                    elif min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends1 - 1:
                            uX_dist = ub_mids1[uX + 1] - b_mids1[x]
                        else:
                            uX_dist = 1000000000
                        for i in range(lY, uY + 1):
                            binned[x, y, 0] += unbinned[uX, i, 0]
                            binned[x, y, 1] += unbinned[uX, i, 1]
                    elif min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids2[y] - ub_mids2[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(lX, uX + 1):
                            binned[x, y, 0] += unbinned[i, lY, 0]
                            binned[x, y, 1] += unbinned[i, lY, 1]
                    elif min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends2 - 1:
                            uY_dist = ub_mids2[uY + 1] - b_mids2[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(lX, uX + 1):
                            binned[x, y, 0] += unbinned[i, uY, 0]
                            binned[x, y, 1] += unbinned[i, uY, 1]
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                if binned[x, y, 0] < minobservations and removefailed == 1:
                    binned[x, y, 0] = 0
                    binned[x, y, 1] = 0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binning_bin_cis_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] filt,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_int64_t, ndim=2] counts,
        np.ndarray[DTYPE_int_t, ndim=3] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] distance_cutoffs,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        int distance_div,
        int distance_bins,
        int mindistance,
        int maxdistance):
    cdef long long int i, j, k, distance, fend1, fend2, index, prev_fend
    cdef long long int num_data = data.shape[0]
    cdef long long int num_features = all_indices.shape[1]
    with nogil:
        prev_fend = -1
        for i in range(num_data):
            fend1 = data[i, 0]
            if filt[fend1] == 0:
                continue
            if fend1 != prev_fend:
                prev_fend = fend1
                k = 0
            fend2 = data[i, 1]
            if filt[fend2] == 0:
                continue
            distance = mids[fend2] - mids[fend1]
            if distance < mindistance or distance > maxdistance:
                continue
            index = 0
            for j in range(num_features):
                if all_indices[fend1, j, 0] < all_indices[fend2, j, 0]:
                    index += (all_indices[fend1, j, 1] + all_indices[fend2, j, 0]) * bin_divs[j]
                else:
                    index += (all_indices[fend2, j, 1] + all_indices[fend1, j, 0]) * bin_divs[j]
            if distance_div > 0:
                while distance > distance_cutoffs[k]:
                    k += 1
                index += k * distance_div
            counts[index, 0] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binning_bin_cis_expected(
        np.ndarray[DTYPE_int_t, ndim=1] filt,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_int64_t, ndim=2] counts,
        np.ndarray[DTYPE_int_t, ndim=3] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] distance_cutoffs,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        int distance_div,
        int distance_bins,
        int mindistance,
        int maxdistance,
        int startfend,
        int stopfend,
        int maxfend):
    cdef long long int i, j, k, distance, fend1, fend2, index
    cdef long long int num_features = all_indices.shape[1]
    with nogil:
        for fend1 in range(startfend, min(stopfend, maxfend - 2)):
            if filt[fend1] == 0:
                continue
            k = 0
            fend2 = fend1 + 2
            distance = mids[fend2] - mids[fend1]
            if filt[fend2] == 1 and distance >= mindistance and distance <= maxdistance:
                index = 0
                for j in range(num_features):
                    if all_indices[fend1, j, 0] < all_indices[fend2, j, 0]:
                        index += (all_indices[fend1, j, 1] + all_indices[fend2, j, 0]) * bin_divs[j]
                    else:
                        index += (all_indices[fend2, j, 1] + all_indices[fend1, j, 0]) * bin_divs[j]
                if distance_div > 0:
                    while distance > distance_cutoffs[k]:
                        k += 1
                    index += k * distance_div
                counts[index, 1] += 1
            for fend2 in range((fend1 / 2 + 2) * 2, maxfend):
                if filt[fend2] == 0:
                    continue
                distance = mids[fend2] - mids[fend1]
                if distance < mindistance or distance > maxdistance:
                    continue
                index = 0
                for j in range(num_features):
                    if all_indices[fend1, j, 0] < all_indices[fend2, j, 0]:
                        index += (all_indices[fend1, j, 1] + all_indices[fend2, j, 0]) * bin_divs[j]
                    else:
                        index += (all_indices[fend2, j, 1] + all_indices[fend1, j, 0]) * bin_divs[j]
                if distance_div > 0:
                    while distance > distance_cutoffs[k]:
                        k += 1
                    index += k * distance_div
                counts[index, 1] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binning_bin_trans_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] filt,
        np.ndarray[DTYPE_int64_t, ndim=2] counts,
        np.ndarray[DTYPE_int_t, ndim=3] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        int distance_div,
        int distance_bins):
    cdef long long int i, j, fend1, fend2, index
    cdef long long int num_data = data.shape[0]
    cdef long long int num_features = all_indices.shape[1]
    with nogil:
        for i in range(num_data):
            fend1 = data[i, 0]
            if filt[fend1] == 0:
                continue
            fend2 = data[i, 1]
            if filt[fend2] == 0:
                continue
            index = 0
            for j in range(num_features):
                if all_indices[fend1, j, 0] < all_indices[fend2, j, 0]:
                    index += (all_indices[fend1, j, 1] + all_indices[fend2, j, 0]) * bin_divs[j]
                else:
                    index += (all_indices[fend2, j, 1] + all_indices[fend1, j, 0]) * bin_divs[j]
            if distance_div > 0:
                index += (distance_bins - 1) * distance_div
            counts[index, 0] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binning_bin_trans_expected(
        np.ndarray[DTYPE_int_t, ndim=1] filt,
        np.ndarray[DTYPE_int64_t, ndim=2] counts,
        np.ndarray[DTYPE_int_t, ndim=3] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        int distance_div,
        int distance_bins,
        int startfend1,
        int stopfend1,
        int startfend2,
        int stopfend2):
    cdef long long int i, j, fend1, fend2, index
    cdef long long int num_features = all_indices.shape[1]
    with nogil:
        for fend1 in range(startfend1, stopfend1):
            if filt[fend1] == 0:
                continue
            for fend2 in range(startfend2, stopfend2):
                if filt[fend2] == 0:
                    continue
                index = 0
                for j in range(num_features):
                    if all_indices[fend1, j, 0] < all_indices[fend2, j, 0]:
                        index += (all_indices[fend1, j, 1] + all_indices[fend2, j, 0]) * bin_divs[j]
                    else:
                        index += (all_indices[fend2, j, 1] + all_indices[fend1, j, 0]) * bin_divs[j]
                if distance_div > 0:
                    index += (distance_bins - 1) * distance_div
                counts[index, 1] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def remap_mrh_data(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2,
        int lowerlimit,
        int upperlimit,
        int offset,
        int invalid,
        int flip):
    cdef long long int i, j, k, map1, map2
    cdef long long int num_data = data.shape[0]
    cdef long long int num_fends = indices.shape[0] - 1
    with nogil:
        i = 0
        if mapping2 is None:   
            for j in range(num_fends - 1):
                if indices[j] == indices[j + 1]:
                    continue
                map1 = mapping[data[indices[j], 0] - offset]
                if map1 < 0:
                    continue
                for k in range(indices[j], indices[j + 1]):
                    if data[k, 1] < upperlimit:
                        map2 = mapping[data[k, 1] - offset]
                        if map2 < 0:
                            continue
                        data[i, 0] = map1
                        data[i, 1] = map2
                        data[i, 2] = data[k, 2]
                        i += 1
        else:   
            for j in range(num_fends):
                if indices[j] == indices[j + 1]:
                    continue
                map1 = mapping[data[indices[j], 0] - offset]
                if map1 == -1:
                    continue               
                for k in range(indices[j], indices[j + 1]):
                    if data[k, 1] < lowerlimit or data[k, 1] >= upperlimit:
                        continue
                    map2 = mapping2[data[k, 1] - lowerlimit]
                    if map2 == -1:
                        continue
                    data[i, 0] = map1
                    data[i, 1] = map2
                    data[i, 2] = data[k, 2]
                    i += 1
            if flip == 1:
                for j in range(i):
                    k = data[j, 0]
                    data[j, 0] = data[j, 1]
                    data[j, 1] = k
    return i


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_mrh_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=2] observed not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2):
    cdef long long int i, j, bin1, bin2
    cdef long long int num_fends = indices.shape[0] - 1
    with nogil:
        if mapping2 is None:
            for i in range(num_fends):
                bin1 = mapping[i]
                for j in range(indices[i], indices[i + 1]):
                    bin2 = mapping[data[j, 0]]
                    observed[bin1, bin2] += data[j, 1]
        else:
            for i in range(num_fends):
                bin1 = mapping[i]
                for j in range(indices[i], indices[i + 1]):
                    bin2 = mapping2[data[j, 0]]
                    observed[bin1, bin2] += data[j, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_mrh_cis_expected(
        np.ndarray[DTYPE_t, ndim=2] expected not None,
        np.ndarray[DTYPE_int_t, ndim=1] fend_nums not None,
        np.ndarray[DTYPE_int_t, ndim=1] binmapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices,
        np.ndarray[DTYPE_t, ndim=2] distance_parameters,
        double chrom_mean,
        int dt_int):
    cdef long long int i, j, k, l, bin1, bin2, num, num_parameters, map1, map2
    cdef double value, corr, distance
    cdef long long int n = expected.shape[0]
    cdef long long int valid_fends = binmapping.shape[0]
    cdef long long int num_fends = mapping.shape[0]
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    with nogil:
        if dt_int == 0:
            # if raw data is requested, only determine the number of valid fends per bin
            for i in range(n):
                for j in range(i, n):
                    expected[i, j] = (obs_indices[i + 1] - obs_indices[i]) * (obs_indices[j + 1] - obs_indices[j])
                    expected[j, i] += expected[i, j]
            for i in range(valid_fends - 1):
                bin1 = binmapping[i]
                num = fend_nums[i]
                if fend_nums[i + 1] - num == 1:
                    bin2 = binmapping[i + 1]
                    expected[bin1, bin2] -= 1
                    expected[bin2, bin1] -= 1
                if num % 2 == 0:
                    if fend_nums[i + 1] - num == 3:
                        bin2 = binmapping[i + 1]
                        expected[bin1, bin2] -= 1
                        expected[bin2, bin1] -= 1
                    elif i < valid_fends - 2 and fend_nums[i + 2] - num == 3:
                        bin2 = binmapping[i + 2]
                        expected[bin1, bin2] -= 1
                        expected[bin2, bin1] -= 1
                    elif i < valid_fends - 3 and fend_nums[i + 3] - num == 3:
                        bin2 = binmapping[i + 3]
                        expected[bin1, bin2] -= 1
                        expected[bin2, bin1] -= 1
        # check if can use shortcutting
        elif dt_int == 1 and fend_indices is None:
            for i in range(n):
                for j in range(i, n):
                    expected[i, j] = correction_sums[i] * correction_sums[j]
                expected[i, i] /= 2.0
            for i in range(num_fends):
                map1 = mapping[i]
                if map1 == -1:
                    continue
                bin1 = binmapping[map1]
                corr = corrections[map1]
                expected[bin1, bin1] -= corr * corr / 2.0
                if i < num_fends - 1:
                    map2 = mapping[i + 1]
                    if map2 > map1:
                        expected[bin1, binmapping[map2]] -= corr * corrections[map2]
                    if i < num_fends - 3 and fend_nums[map1] % 2 == 0:
                        map2 = mapping[i + 3]
                        if map2 > map1:
                            expected[bin1, binmapping[map2]] -= corr * corrections[map2]
        else:
            # if distance correction is requested, use chrom mean, else 1.0
            for i in range(valid_fends - 1):
                bin1 = binmapping[i]
                num = fend_nums[i]
                l = 0
                for j in range(i + 1, min(valid_fends, i + 4)):
                    if fend_nums[i + 1] - num == 1:
                        continue
                    if num % 2 == 0 and fend_nums[j] - num == 3:
                        continue
                    if dt_int < 2:
                        value = 1.0
                    else:                            
                        distance = log(<double>(mids[j] - mids[i]))
                        while distance > distance_parameters[l, 0]:
                            l += 1
                        value = exp(distance * distance_parameters[l, 1] + distance_parameters[l, 2] + chrom_mean)
                    if not corrections is None:
                        value *= corrections[i] * corrections[j]
                    if not fend_indices is None:
                        for k in range(num_parameters):
                            if fend_indices[i, k, 0] < fend_indices[j, k, 0]:
                                value *= binning_corrections[fend_indices[i, k, 1] + fend_indices[j, k, 0]]
                            else:
                                value *= binning_corrections[fend_indices[i, k, 0] + fend_indices[j, k, 1]]
                    expected[bin1, binmapping[j]] += value
                for j in range(i + 4, valid_fends):
                    if dt_int < 2:
                        value = 1.0
                    else:                            
                        distance = log(<double>(mids[j] - mids[i]))
                        while distance > distance_parameters[l, 0]:
                            l += 1
                        value = exp(distance * distance_parameters[l, 1] + distance_parameters[l, 2] + chrom_mean)
                    if not corrections is None:
                        value *= corrections[i] * corrections[j]
                    if not fend_indices is None:
                        for k in range(num_parameters):
                            if fend_indices[i, k, 0] < fend_indices[j, k, 0]:
                                value *= binning_corrections[fend_indices[i, k, 1] + fend_indices[j, k, 0]]
                            else:
                                value *= binning_corrections[fend_indices[i, k, 0] + fend_indices[j, k, 1]]
                    expected[bin1, binmapping[j]] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_mrh_trans_expected(
        np.ndarray[DTYPE_t, ndim=2] expected not None,
        np.ndarray[DTYPE_int_t, ndim=1] binmapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] binmapping2 not None,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices2 not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_t, ndim=1] corrections2,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        np.ndarray[DTYPE_t, ndim=1] correction_sums2,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices2,
        double chrom_mean,
        int dt_int):
    cdef long long int i, j, k, bin1, num_parameters
    cdef double value
    cdef long long int n = expected.shape[0]
    cdef long long int m = expected.shape[1]
    cdef long long int num_fends = binmapping.shape[0]
    cdef long long int num_fends2 = binmapping2.shape[0]
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    with nogil:
        if dt_int == 0:
            # if raw data is requested, only determine the number of valid fends per bin
            for i in range(n):
                for j in range(m):
                    expected[i, j] = (obs_indices[i + 1] - obs_indices[i]) * (obs_indices2[j + 1] - obs_indices2[j])
        else:
            # check if can use shortcutting
            if dt_int == 1 and fend_indices is None:
                for i in range(n):
                    for j in range(m):
                        expected[i, j] += correction_sums[i] * correction_sums2[j]
            else:
                # if distance correction is requested, use chrom mean, else 1.0
                for i in range(num_fends):
                    bin1 = binmapping[i]
                    for j in range(num_fends2):
                        if dt_int < 2:
                            value = 1.0
                        else:
                            value = chrom_mean
                        if not corrections is None:
                            value *= corrections[i] * corrections2[j]
                        if not fend_indices is None:
                            for k in range(num_parameters):
                                if fend_indices[i, k, 0] < fend_indices[j, k, 0]:
                                    value *= binning_corrections[fend_indices[i, k, 1] + fend_indices[j, k, 0]]
                                else:
                                    value *= binning_corrections[fend_indices[i, k, 0] + fend_indices[j, k, 1]]
                        expected[bin1, binmapping2[j]] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def make_trans_mrh_toplevel(
        np.ndarray[DTYPE_int_t, ndim=2] observed,
        np.ndarray[DTYPE_t, ndim=2] expected,
        np.ndarray[DTYPE_t, ndim=1] current_level_data,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices2,
        np.ndarray[DTYPE_int_t, ndim=1] bin_position,
        int minobservations):
    cdef long long int i, j, k, l, index
    cdef double value, log_2
    cdef int n_bins = obs_indices.shape[0] - 1
    cdef int m_bins = obs_indices2.shape[0] - 1
    cdef double nan = numpy.nan
    with nogil:
        log_2 = log(2.0)
        index = 0
        for i in range(n_bins):
            for j in range(m_bins):
                bin_position[index] = i * m_bins + j
                for k in range(obs_indices[i], obs_indices[i + 1]):
                    for l in range(obs_indices2[j], obs_indices2[j + 1]):
                        current_level_data[index] += observed[k, l]
                if current_level_data[index] >= minobservations:
                    # find expected values
                    value = 0.0
                    for k in range(obs_indices[i], obs_indices[i + 1]):
                        for l in range(obs_indices2[j], obs_indices2[j + 1]):
                            value += expected[k, l]
                    current_level_data[index] = log(current_level_data[index] / value) / log_2
                else:
                    current_level_data[index] = nan
                index += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def make_cis_mrh_toplevel(
        np.ndarray[DTYPE_int_t, ndim=2] observed,
        np.ndarray[DTYPE_t, ndim=2] expected,
        np.ndarray[DTYPE_t, ndim=1] current_level_data,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices,
        np.ndarray[DTYPE_int_t, ndim=1] bin_position,
        int minobservations):
    cdef long long int i, j, k, l, index
    cdef double value, log_2
    cdef long long int n_bins = obs_indices.shape[0] - 1
    cdef double nan = numpy.nan
    with nogil:
        index = 0
        log_2 = log(2.0)
        for i in range(n_bins):
            for j in range(i, n_bins):
                bin_position[index] = i * n_bins + j
                for k in range(obs_indices[i], obs_indices[i + 1]):
                    for l in range(obs_indices[j], obs_indices[j + 1]):
                        current_level_data[index] += observed[k, l]
                if current_level_data[index] >= minobservations:
                    # find expected values
                    value = 0.0
                    for k in range(obs_indices[i], obs_indices[i + 1]):
                        for l in range(obs_indices[j], obs_indices[j + 1]):
                            value += expected[k, l]
                    current_level_data[index] = log(current_level_data[index] / value) / log_2
                else:
                    current_level_data[index] = nan
                index += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def make_trans_mrh_midlevel(
        np.ndarray[DTYPE_int_t, ndim=2] observed,
        np.ndarray[DTYPE_t, ndim=2] expected,
        np.ndarray[DTYPE_t, ndim=1] current_level_data,
        np.ndarray[DTYPE_t, ndim=1] prev_level_data,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_indices,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_shapes,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices2,
        np.ndarray[DTYPE_int_t, ndim=1] prev_bin_position,
        np.ndarray[DTYPE_int_t, ndim=1] bin_position,
        int prev_m_bins,
        int m_bins,
        int minobservations,
        int pos):
    cdef long long int i, j, k, l, m, valid, shape, pos2, index, start1, stop1, start2, stop2, value1
    cdef double value2, log_2
    cdef long long int num_prev_data = prev_level_data.shape[0]
    cdef double nan = numpy.nan
    with nogil:
        pos2 = 0
        log_2 = log(2.0)
        for i in range(num_prev_data):
            if prev_level_data[i] == nan:
                prev_level_indices[i] = -1
            else:
                valid = 0
                # find position of previous resolution bin
                start1 = (prev_bin_position[i] / prev_m_bins) * 2
                stop1 = start1 + 2
                start2 = (prev_bin_position[i] % prev_m_bins) * 2
                stop2 = start2 + 2
                # find new resolution data for this set of bins
                index = 1
                shape = 0
                for j in range(start1, stop1):
                    for k in range(start2, stop2):
                        value1 = 0
                        value2 = 0.0
                        for l in range(obs_indices[j], obs_indices[j + 1]):
                            for m in range(obs_indices2[k], obs_indices2[k + 1]):
                                value1 += observed[l, m]
                                value2 += expected[l, m]
                        if value1 >= minobservations:
                            current_level_data[pos2] = log(value1 / value2) / log_2
                            bin_position[pos2] = j * m_bins + k
                            shape += index
                            valid += 1
                            pos2 += 1
                        index *= 2
                # if there are valid values at this new resolution, add the values to the new data set
                if valid > 0:
                    prev_level_indices[i] = pos2 - valid + pos
                    prev_level_shapes[i] = shape
                else:
                    prev_level_indices[i] = -1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def make_cis_mrh_midlevel(
        np.ndarray[DTYPE_int_t, ndim=2] observed,
        np.ndarray[DTYPE_t, ndim=2] expected,
        np.ndarray[DTYPE_t, ndim=1] current_level_data,
        np.ndarray[DTYPE_t, ndim=1] prev_level_data,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_indices,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_shapes,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices,
        np.ndarray[DTYPE_int_t, ndim=1] prev_bin_position,
        np.ndarray[DTYPE_int_t, ndim=1] bin_position,
        int prev_n_bins,
        int n_bins,
        int minobservations,
        int pos):
    cdef long long int i, j, k, l, m, valid, pos2, shape, index, start1, stop1, start2, stop2, value1
    cdef double value2, log_2
    cdef long long int num_prev_data = prev_level_data.shape[0]
    cdef double nan = numpy.nan
    with nogil:
        pos2 = 0
        log_2 = log(2.0)
        for i in range(num_prev_data):
            if prev_level_data[i] == nan:
                prev_level_indices[i] = -1
            else:
                # find position of previous resolution bin
                start1 = (prev_bin_position[i] / prev_n_bins) * 2
                stop1 = start1 + 2
                start2 = (prev_bin_position[i] % prev_n_bins) * 2
                stop2 = start2 + 2
                # find new resolution data for this set of bins
                valid = 0
                index = 1
                shape = 0
                for j in range(start1, stop1):
                    for k in range(start2, stop2):
                        # since it's possible to be below the diagonal, check and skip if necessary
                        if j > k:
                            index *= 2
                            continue
                        value1 = 0
                        value2 = 0.0
                        for l in range(obs_indices[j], obs_indices[j + 1]):
                            for m in range(obs_indices[k], obs_indices[k + 1]):
                                value1 += observed[l, m]
                                value2 += expected[l, m]
                        if value1 >= minobservations:
                            current_level_data[pos2] = log(value1 / value2) / log_2
                            bin_position[pos2] = j * n_bins + k
                            shape += index
                            valid += 1
                            pos2 += 1
                        index *= 2
                # if there are valid values at this new resolution, add the values to the new data set
                if valid > 0:
                    prev_level_indices[i] = pos2 - valid + pos
                    prev_level_shapes[i] = shape
                else:
                    prev_level_indices[i] = -1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def make_trans_mrh_lowerlevel(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] data_indices,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        np.ndarray[DTYPE_t, ndim=1] correction_sums2,
        np.ndarray[DTYPE_t, ndim=1] current_level_data,
        np.ndarray[DTYPE_t, ndim=1] prev_level_data,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_indices,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_shapes,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices2,
        np.ndarray[DTYPE_int_t, ndim=1] prev_bin_position,
        np.ndarray[DTYPE_int_t, ndim=1] bin_position,
        int prev_m_bins,
        int m_bins,
        int minobservations,
        int pos):
    cdef long long int i, j, k, l, m, valid, pos2, shape, index, index1, index2
    cdef long long int start1, stop1, start2, stop2, start, mid, stop, value1
    cdef double value2, log_2
    cdef long long int num_prev_data = prev_level_data.shape[0]
    cdef double nan = numpy.nan
    with nogil:
        log_2 = log(2.0)
        pos2 = 0
        for i in range(num_prev_data):
            if prev_level_data[i] == nan:
                prev_level_indices[i] = -1
            else:
                # find position of previous resolution bin
                start1 = (prev_bin_position[i] / prev_m_bins) * 2
                stop1 = start1 + 2
                start2 = (prev_bin_position[i] % prev_m_bins) * 2
                stop2 = start2 + 2
                # find new resolution data for this set of bins
                valid = 0
                index = 1
                shape = 0
                m = 0
                for j in range(start1, stop1):
                    for k in range(start2, stop2):
                        value1 = 0
                        mid = 0
                        for l in range(obs_indices[j], obs_indices[j + 1]):
                            start = data_indices[l]
                            stop = data_indices[l + 1]
                            index1 = obs_indices2[k]
                            index2 = obs_indices2[k + 1]
                            m = min(mid + start, stop)
                            while m > start and data[m - 1, 0] >= index1:
                                m -= 1
                            while m < stop and data[m, 0] < index1:
                                m += 1
                            mid = m - start
                            while m < stop and data[m, 0] < index2:
                                value1 += data[m, 1]
                                m += 1
                        if value1 >= minobservations:
                            # find expected value
                            if not correction_sums is None and not correction_sums2 is None:
                                value2 = correction_sums[j] * correction_sums2[k]
                            else:
                                value2 = 1.0
                            current_level_data[pos2] = log(value1 / value2) / log_2
                            bin_position[pos2] = j * m_bins + k
                            shape += index
                            valid += 1
                            pos2 += 1
                        index *= 2
                # if there are valid values at this new resolution, add the values to the new data set
                if valid > 0:
                    prev_level_indices[i] = pos2 - valid + pos
                    prev_level_shapes[i] = shape
                else:
                    prev_level_indices[i] = -1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def make_cis_mrh_lowerlevel(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] data_indices,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_t, ndim=1] correction_sums,
        np.ndarray[DTYPE_int_t, ndim=1] fend_nums,
        np.ndarray[DTYPE_t, ndim=1] current_level_data,
        np.ndarray[DTYPE_t, ndim=1] prev_level_data,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_indices,
        np.ndarray[DTYPE_int_t, ndim=1] prev_level_shapes,
        np.ndarray[DTYPE_int_t, ndim=1] obs_indices,
        np.ndarray[DTYPE_int_t, ndim=1] prev_bin_position,
        np.ndarray[DTYPE_int_t, ndim=1] bin_position,
        int prev_n_bins,
        int n_bins,
        int minobservations,
        int pos):
    cdef long long int i, j, k, l, m, valid, pos2, shape, index, index1, index2
    cdef long long int start1, stop1, start2, stop2, start, mid, stop, value1
    cdef double value2, log_2
    cdef long long int num_prev_data = prev_level_data.shape[0]
    cdef double nan = numpy.nan
    with nogil:
        pos2 = 0
        log_2 = log(2.0)
        for i in range(num_prev_data):
            if prev_level_data[i] == nan:
                prev_level_indices[i] = -1
            else:
                # find position of previous resolution bin
                start1 = (prev_bin_position[i] / prev_n_bins) * 2
                stop1 = start1 + 2
                start2 = (prev_bin_position[i] % prev_n_bins) * 2
                stop2 = start2 + 2
                # find new resolution data for this set of bins
                valid = 0
                index = 1
                shape = 0
                for j in range(start1, stop1):
                    for k in range(start2, stop2):
                        # since it's possible to be below the diagonal, check and skip if necessary
                        if j > k:
                            index *= 2
                            continue
                        value1 = 0
                        mid = 0
                        for l in range(obs_indices[j], obs_indices[j + 1]):
                            start = data_indices[l]
                            stop = data_indices[l + 1]
                            index1 = obs_indices[k]
                            index2 = obs_indices[k + 1]
                            m = min(mid + start, stop - 1)
                            while m > start and data[m, 0] > index1:
                                m -= 1
                            while m < stop and data[m, 0] < index1:
                                m += 1
                            mid = m - start
                            while m < stop and data[m, 0] < index2:
                                value1 += data[m, 1]
                                m += 1
                        if value1 >= minobservations:
                            # find expected value
                            if not correction_sums is None:
                                value2 = correction_sums[j] * correction_sums[k]
                                l = obs_indices[j + 1] - 1
                                m = obs_indices[k]
                                start = obs_indices[j]
                                stop = obs_indices[k + 1]
                                if fend_nums[l] + 1 == fend_nums[m]:
                                    value2 -= corrections[l] * corrections[m]
                                if fend_nums[l] % 2 == 0:
                                    if fend_nums[l] + 3 == fend_nums[m]:
                                        value2 -= corrections[l] * corrections[m]
                                    elif m + 1 < stop and fend_nums[l] + 3 == fend_nums[m + 1]:
                                        value2 -= corrections[l] * corrections[m + 1]
                                    elif m + 2 < stop and fend_nums[l] + 3 == fend_nums[m + 2]:
                                        value2 -= corrections[l] * corrections[m + 2]
                                l -= 1
                                if l >= start and fend_nums[l] % 2 == 0:
                                    if fend_nums[l] + 3 == fend_nums[m]:
                                        value2 -= corrections[l] * corrections[m]
                                    elif m + 1 < stop and fend_nums[l] + 3 == fend_nums[m + 1]:
                                        value2 -= corrections[l] * corrections[m + 1]
                                l -= 1
                                if l >= start and fend_nums[l] % 2 == 0:
                                    if fend_nums[l] + 3 == fend_nums[m]:
                                        value2 -= corrections[l] * corrections[m]
                            else:
                                value2 = 1.0
                            current_level_data[pos2] = log(value1 / value2) / log_2
                            bin_position[pos2] = j * n_bins + k
                            shape += index
                            valid += 1
                            pos2 += 1
                        index *= 2
                # if there are valid values at this new resolution, add the values to the new data set
                if valid > 0:
                    prev_level_indices[i] = pos2 - valid + pos
                    prev_level_shapes[i] = shape
                else:
                    prev_level_indices[i] = -1
    return None
