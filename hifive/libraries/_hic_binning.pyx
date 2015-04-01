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
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None):
    cdef long long int fend1, fend2, i, map1, map2
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        for fend1 in range(num_fends - 1):
            map1 = mapping[fend1]
            if map1 == -1:
                continue
            for i in range(indices[fend1], indices[fend1 + 1]):
                fend2 = data[i, 1]
                if fend2 >= num_fends:
                    continue
                map2 = mapping[fend2]
                if map2 == -1 or map1 == map2 or fend2 >= max_fend[fend1]:
                    continue
                signal[map1, map2 - map1 - 1, 0] += data[i, 2]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_upper_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None):
    cdef long long int fend1, fend2, i, index, map1, map2
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5))
    with nogil:
        for fend1 in range(num_fends - 1):
            map1 = mapping[fend1]
            if map1 == -1:
                continue
            index = map1 * (num_bins - 1) - map1 * (map1 + 1) / 2 - 1
            for i in range(indices[fend1], indices[fend1 + 1]):
                fend2 = data[i, 1]
                if fend2 >= num_fends:
                    continue
                map2 = mapping[fend2]
                if map2 == -1 or map2 == map1 or fend2 >= max_fend[fend1]:
                    continue
                signal[index + map2, 0] += data[i, 2]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_compact_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] correction_indices,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins,
        np.ndarray[DTYPE_int_t, ndim=2] fend_indices,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double chrom_mean,
        int startfend):
    cdef long long int fend1, fend2, j, k, map1, map2, bin1, bin2, index, num_parameters
    cdef double distance, value
    cdef long long int num_fends = mapping.shape[0]
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    else:
        num_parameters = 0
    with nogil:
        for fend1 in range(num_fends - 1):
            map1 = mapping[fend1]
            if map1 == -1:
                continue
            k = 0
            # find opposite strand adjacents, skipping same fragment and same strand adjacents
            fend2 = fend1 + 2
            map2 = mapping[fend2]
            if map2 >= 0 and map2 != map1 and fend2 < max_fend[fend1]:
                 # give starting expected value
                value = 1.0
                # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                if not corrections is None:
                    value *= corrections[fend1] * corrections[fend2]
                # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                if not binning_corrections is None:
                    for j in range(num_parameters):
                        bin1 = min(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        bin2 = max(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        index = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value *= binning_corrections[index]
                # if finding distance, enrichment, or expected, correct for distance
                if not parameters is None:
                    distance = log(mids[fend2] - mids[fend1])
                    while distance > parameters[k, 0]:
                        k += 1
                    value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                signal[map1, map2 - map1 - 1, 1] += value
            for fend2 in range((fend1 / 2) * 2 + 4, min(max_fend[fend1], num_fends)):
                map2 = mapping[fend2]
                if map2 == -1 or map2 == map1:
                    continue
                 # give starting expected value
                value = 1.0
                # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                if not corrections is None:
                    value *= corrections[fend1] * corrections[fend2]
                # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                if not binning_corrections is None:
                    for j in range(num_parameters):
                        bin1 = min(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        bin2 = max(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        index = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value *= binning_corrections[index]
                # if finding distance, enrichment, or expected, correct for distance
                if not parameters is None:
                    distance = log(mids[fend2] - mids[fend1])
                    while distance > parameters[k, 0]:
                        k += 1
                    value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                signal[map1, map2 - map1 - 1, 1] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_upper_expected(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] correction_indices,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins,
        np.ndarray[DTYPE_int_t, ndim=2] fend_indices,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        double chrom_mean,
        int startfend):
    cdef long long int fend1, fend2, j, k, index, map1, map2, bin1, bin2, index2, num_parameters
    cdef double distance, value
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5))
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    else:
        num_parameters = 0
    with nogil:
        for fend1 in range(num_fends - 1):
            map1 = mapping[fend1]
            if map1 == -1:
                continue
            k = 0
            index = map1 * (num_bins - 1) - map1 * (map1 + 1) / 2 - 1
            # find opposite strand adjacents, skipping same fragment and same strand adjacents
            fend2 = fend1 + 2
            map2 = mapping[fend2]
            if map2 >= 0 and map2 != map1 and fend2 < max_fend[fend1]:
                 # give starting expected value
                value = 1.0
                # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                if not corrections is None:
                    value *= corrections[fend1] * corrections[fend2]
                # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                if not binning_corrections is None:
                    for j in range(num_parameters):
                        bin1 = min(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        bin2 = max(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        index2 = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value *= binning_corrections[index2]
                # if finding distance, enrichment, or expected, correct for distance
                if not parameters is None:
                    distance = log(mids[fend2] - mids[fend1])
                    while distance > parameters[k, 0]:
                        k += 1
                    value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                signal[index + map2, 1] += value
            for fend2 in range((fend1 / 2) * 2 + 4, min(max_fend[fend1], num_fends)):
                map2 = mapping[fend2]
                if map2 == -1 or map2 == map1:
                    continue
                 # give starting expected value
                value = 1.0
                # if finding fend, enrichment, or expected, and using express or probability bias correction, correct for fend
                if not corrections is None:
                    value *= corrections[fend1] * corrections[fend2]
                # if finding fend, enrichment, or expected, and using binning bias correction, correct for fend
                if not binning_corrections is None:
                    for j in range(num_parameters):
                        bin1 = min(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        bin2 = max(fend_indices[fend1 + startfend, j], fend_indices[fend2 + startfend, j]) 
                        index2 = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value *= binning_corrections[index2]
                # if finding distance, enrichment, or expected, correct for distance
                if not parameters is None:
                    distance = log(mids[fend2] - mids[fend1])
                    while distance > parameters[k, 0]:
                        k += 1
                    value *= exp(distance * parameters[k, 1] + parameters[k, 2] + chrom_mean)
                signal[index + map2, 1] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_compact_to_compact(
        np.ndarray[DTYPE_t, ndim=3] binned,
        np.ndarray[DTYPE_t, ndim=3] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping):
    cdef long long int i, j, fend2
    cdef long long int num_bins = binned.shape[0]
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    with nogil:
        for i in range(num_fends - 1):
            if mapping[i] == -1:
                continue
            for j in range(min(num_fends - i - 1, max_fend)):
                fend2 = j + i + 1
                if mapping[fend2] == -1 or mapping[fend2] == mapping[i]:
                    continue
                else:
                    binned[mapping[i], mapping[fend2] - mapping[i] - 1, 0] += unbinned[i, j, 0]
                    binned[mapping[i], mapping[fend2] - mapping[i] - 1, 1] += unbinned[i, j, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_compact_to_upper(
        np.ndarray[DTYPE_t, ndim=2] binned,
        np.ndarray[DTYPE_t, ndim=3] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        int num_bins):
    cdef long long int i, j, index
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    with nogil:
        for i in range(num_fends - 1):
            if mapping[i] == -1:
                continue
            index = mapping[i] * num_bins - mapping[i] * (mapping[i] + 1) / 2 - mapping[i] - 1
            for j in range(i + 1, num_fends):
                if mapping[j] == -1 or mapping[j] == mapping[i]:
                    continue
                else:
                    binned[index + mapping[j], 0] += unbinned[i, j - i - 1, 0]
                    binned[index + mapping[j], 1] += unbinned[i, j - i - 1, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_upper_to_compact(
        np.ndarray[DTYPE_t, ndim=3] binned,
        np.ndarray[DTYPE_t, ndim=2] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping):
    cdef long long int i, j, index
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        for i in range(num_fends - 1):
            if mapping[i] == -1:
                continue
            index = i * num_fends - i * (i + 1) / 2 - i - 1
            for j in range(i + 1, num_fends):
                if mapping[j] == -1 or mapping[j] == mapping[i]:
                    continue
                else:
                    binned[mapping[i], mapping[j] - mapping[i] - 1, 0] += unbinned[index + j, 0]
                    binned[mapping[i], mapping[j] - mapping[i] - 1, 1] += unbinned[index + j, 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_upper_to_upper(
        np.ndarray[DTYPE_t, ndim=2] binned,
        np.ndarray[DTYPE_t, ndim=2] unbinned,
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        int num_bins):
    cdef long long int i, j, index, index2
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        for i in range(num_fends - 1):
            if mapping[i] == -1:
                continue
            index = i * num_fends - i * (i + 1) / 2 - i - 1
            index2 = mapping[i] * num_bins - mapping[i] * (mapping[i] + 1) / 2 - mapping[i] - 1
            for j in range(i + 1, num_fends):
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
        int removefailed):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist, index, index2
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int num_fends = ub_mids.shape[0]
    with nogil:
        for x in range(num_bins - 1):
            index = x * num_bins - x * (x + 1) / 2 - x - 1
            for y in range(x + 1, num_bins):
                # if bin already meets our criteria, skip
                if binned[index + y, 0] >= minobservations:
                    continue
                # otherwise, set bordering unbinned positions according to bounds
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
                        index2 = lX * num_fends - lX * (lX + 1) / 2 - lX - 1
                        for i in range(max(lX + 1, lY), uY + 1):
                            binned[index + y, 0] += unbinned[index2 + i, 0]
                            binned[index + y, 1] += unbinned[index2 + i, 1]
                    elif min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        index2 = uX * num_fends - uX * (uX + 1) / 2 - uX - 1
                        for i in range(max(uX + 1, lY), uY + 1):
                            binned[index + y, 0] += unbinned[index2 + i, 0]
                            binned[index + y, 1] += unbinned[index2 + i, 1]
                    elif min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(lX, min(uX + 1, lY)):
                            index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                            binned[index + y, 0] += unbinned[index2 + lY, 0]
                            binned[index + y, 1] += unbinned[index2 + lY, 1]
                    elif min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(lX, min(uX + 1, uY)):
                            index2 = i * num_fends - i * (i + 1) / 2 - i - 1
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
        int removefailed):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist, index
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int max_bin = binned.shape[1]
    cdef long long int num_fends = ub_mids.shape[0]
    with nogil:
        for x in range(num_bins - 1):
            for y in range(x + 1, min(x + max_bin, num_bins)):
                # if bin already meets our criteria, skip
                if binned[x, y - x - 1, 0] >= minobservations:
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
                while (maxsearch == 0 or min_dist < maxsearch) and binned[x, y - x - 1, 0] < minobservations:
                    # find min dist, update distance, and add new row or col of observations
                    if min_dist == lX_dist:
                        lX -= 1
                        if lX > 0:
                            lX_dist = b_mids[x] - ub_mids[lX - 1]
                        else:
                            lX_dist = 1000000000
                        index = lX * num_fends - lX * (lX + 1) / 2 - lX - 1
                        for i in range(max(lX + 1, lY), uY + 1):
                            binned[x, y - x - 1, 0] += unbinned[index + i, 0]
                            binned[x, y - x - 1, 1] += unbinned[index + i, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        index = uX * num_fends - uX * (uX + 1) / 2 - uX - 1
                        for i in range(max(uX + 1, lY), uY + 1):
                            binned[x, y - x - 1, 0] += unbinned[index + i, 0]
                            binned[x, y - x - 1, 1] += unbinned[index + i, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(lX, min(uX + 1, lY)):
                            index = i * num_fends - i * (i + 1) / 2 - i - 1
                            binned[x, y - x - 1, 0] += unbinned[index + lY, 0]
                            binned[x, y - x - 1, 1] += unbinned[index + lY, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(lX, min(uX + 1, uY)):
                            index = i * num_fends - i * (i + 1) / 2 - i - 1
                            binned[x, y - x - 1, 0] += unbinned[index + uY, 0]
                            binned[x, y - x - 1, 1] += unbinned[index + uY, 1]
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                if binned[x, y - x - 1, 0] < minobservations and removefailed == 1:
                    binned[x, y - x - 1, 0] = 0
                    binned[x, y - x - 1, 1] = 0
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
        int removefailed):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist, index
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int num_fends = ub_mids.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    with nogil:
        for x in range(num_bins - 1):
            index = x * num_bins - x * (x + 1) / 2 - x - 1
            for y in range(x + 1, num_bins):
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
                        for i in range(max(lX + 1, lY), min(uY + 1, lX + max_fend + 1)):
                            binned[index + y, 0] += unbinned[lX, i - lX - 1, 0]
                            binned[index + y, 1] += unbinned[lX, i - lX - 1, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        for i in range(max(uX + 1, lY), min(uY + 1, uX + max_fend + 1)):
                            binned[index + y, 0] += unbinned[uX, i - uX - 1, 0]
                            binned[index + y, 1] += unbinned[uX, i - uX - 1, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(max(lX, lY - max_fend - 1), min(uX + 1, lY)):
                            binned[index + y, 0] += unbinned[i, lY - i - 1, 0]
                            binned[index + y, 1] += unbinned[i, lY - i - 1, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(max(lX, uY - max_fend - 1), min(uX + 1, uY)):
                            binned[index + y, 0] += unbinned[i, uY - i - 1, 0]
                            binned[index + y, 1] += unbinned[i, uY - i - 1, 1]
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
        int removefailed):
    cdef long long int x, y, i, lX, lX_dist, uX, uX_dist, lY, lY_dist, uY, uY_dist, min_dist
    cdef long long int num_bins = bounds.shape[0]
    cdef long long int max_bin = binned.shape[1]
    cdef long long int num_fends = ub_mids.shape[0]
    cdef long long int max_fend = unbinned.shape[1]
    with nogil:
        for x in range(num_bins - 1):
            for y in range(x + 1, min(x + max_bin, num_bins)):
                # if bin already meets our criteria, skip
                if binned[x, y - x - 1, 0] >= minobservations:
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
                while (maxsearch == 0 or min_dist < maxsearch) and binned[x, y - x - 1, 0] < minobservations:
                    # find min dist, update distance, and add new row or col of observations
                    if min_dist == lX_dist:
                        lX -= 1
                        if lX > 0:
                            lX_dist = b_mids[x] - ub_mids[lX - 1]
                        else:
                            lX_dist = 1000000000
                        for i in range(max(lX + 1, lY), min(uY + 1, lX + max_fend + 1)):
                            binned[x, y - x - 1, 0] += unbinned[lX, i - lX - 1, 0]
                            binned[x, y - x - 1, 1] += unbinned[lX, i - lX - 1, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        for i in range(max(uX + 1, lY), min(uY + 1, uX + max_fend + 1)):
                            binned[x, y - x - 1, 0] += unbinned[uX, i - uX - 1, 0]
                            binned[x, y - x - 1, 1] += unbinned[uX, i - uX - 1, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(max(lX, lY - max_fend - 1), min(uX + 1, lY)):
                            binned[x, y - x - 1, 0] += unbinned[i, lY - i - 1, 0]
                            binned[x, y - x - 1, 1] += unbinned[i, lY - i - 1, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_fends - 1:
                            uY_dist = ub_mids[uY + 1] - b_mids[y]
                        else:
                            uY_dist = 1000000000
                        for i in range(max(lX, uY - max_fend - 1), min(uX + 1, uY)):
                            binned[x, y - x - 1, 0] += unbinned[i, uY - i - 1, 0]
                            binned[x, y - x - 1, 1] += unbinned[i, uY - i - 1, 1]
                    # determine min distance
                    min_dist = min(min(lX_dist, uX_dist), min(lY_dist, uY_dist))
                if binned[x, y - x - 1, 0] < minobservations and removefailed == 1:
                    binned[x, y - x - 1, 0] = 0
                    binned[x, y - x - 1, 1] = 0
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
                if fend2 >= num_fends2 or mapping2[fend2] == -1:
                    i += 1
                    continue
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
        np.ndarray[DTYPE_int_t, ndim=1] correction_indices,
        np.ndarray[DTYPE_int_t, ndim=1] binning_num_bins,
        np.ndarray[DTYPE_int_t, ndim=2] fend_indices,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double trans_mean,
        int startfend1,
        int startfend2):
    cdef long long int j, fend1, fend2, map1, map2, bin1, bin2, index, num_parameters
    cdef double value
    cdef long long int num_fends1 = mapping1.shape[0]
    cdef long long int num_fends2 = mapping2.shape[0]
    if not fend_indices is None:
        num_parameters = fend_indices.shape[1]
    else:
        num_parameters = 0
    with nogil:
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
                        bin1 = min(fend_indices[fend1 + startfend1, j], fend_indices[fend2 + startfend2, j]) 
                        bin2 = max(fend_indices[fend1 + startfend1, j], fend_indices[fend2 + startfend2, j]) 
                        index = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value *= binning_corrections[index]
                signal[map1, map2, 1] += value
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
        np.ndarray[DTYPE_int_t, ndim=2] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] distance_cutoffs,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        int distance_div,
        int distance_bins,
        int mindistance,
        int maxdistance):
    cdef long long int i, j, k, distance, fend1, fend2, index, bin1, bin2, prev_fend
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
                bin1 = min(all_indices[fend1, j], all_indices[fend2, j])
                bin2 = max(all_indices[fend1, j], all_indices[fend2, j])
                index += ((num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2) * bin_divs[j]
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
        np.ndarray[DTYPE_int_t, ndim=2] all_indices,
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
    cdef long long int i, j, k, distance, fend1, fend2, index, bin1, bin2
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
                    bin1 = min(all_indices[fend1, j], all_indices[fend2, j])
                    bin2 = max(all_indices[fend1, j], all_indices[fend2, j])
                    index += ((num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2) * bin_divs[j]
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
                    bin1 = min(all_indices[fend1, j], all_indices[fend2, j])
                    bin2 = max(all_indices[fend1, j], all_indices[fend2, j])
                    index += ((num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2) * bin_divs[j]
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
        np.ndarray[DTYPE_int_t, ndim=2] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        int distance_div,
        int distance_bins):
    cdef long long int i, j, fend1, fend2, index, bin1, bin2
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
                bin1 = min(all_indices[fend1, j], all_indices[fend2, j])
                bin2 = max(all_indices[fend1, j], all_indices[fend2, j])
                index += ((num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2) * bin_divs[j]
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
        np.ndarray[DTYPE_int_t, ndim=2] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        int distance_div,
        int distance_bins,
        int startfend1,
        int stopfend1,
        int startfend2,
        int stopfend2):
    cdef long long int i, j, fend1, fend2, index, bin1, bin2
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
                    bin1 = min(all_indices[fend1, j], all_indices[fend2, j])
                    bin2 = max(all_indices[fend1, j], all_indices[fend2, j])
                    index += ((num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2) * bin_divs[j]
                if distance_div > 0:
                    index += (distance_bins - 1) * distance_div
                counts[index, 1] += 1
    return None
