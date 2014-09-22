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
def find_distance_mapping(
        np.ndarray[DTYPE_int_t, ndim=2] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        int first_fend,
        int start,
        int binsize,
        int stepsize):
    cdef long long int start_fend, stop_fend, start_pos, stop_pos, i
    cdef long long int num_bins = mapping.shape[0]
    cdef long long int num_fends = mids.shape[0]
    with nogil:
        start_fend = 0
        stop_fend = 0
        start_pos = start
        stop_pos = start + binsize
        for i in range(num_bins):
            while start_fend < num_fends - 1 and mids[start_fend] < start_pos:
                start_fend += 1
            stop_fend = max(stop_fend, start_fend)
            while stop_fend < num_fends - 1 and mids[stop_fend] < stop_pos:
                stop_fend += 1
            mapping[i, 0] = start_fend + first_fend
            mapping[i, 1] = stop_fend + first_fend
            start_pos += stepsize
            stop_pos += stepsize
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fend_coverage(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] data_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_int_t, ndim=1] coverage not None,
        int mincoverage):
    cdef long long int i, fend1, valid
    cdef long long int num_fends = filter.shape[0]
    with nogil:
        valid = 0
        for fend1 in range(num_fends):
            if filter[fend1] == 0:
                continue
            i = data_indices[fend1]
            while i < data_indices[fend1 + 1] and data[i, 1] < max_fend[fend1]:
                if filter[data[i, 1]] == 1:
                    coverage[fend1] += 1
                    coverage[data[i, 1]] += 1
                i += 1
        for i in range(num_fends):
            if coverage[i] < mincoverage:
                filter[i] = 0
            else:
                valid += 1
    return valid


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def unbinned_signal_compact(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_means not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        int datatype,
        int startfend):
    cdef long long int fend1, fend2, i, j, k, distance, last_index
    cdef double frac
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int num_distance_bins = distance_mids.shape[0]
    with nogil:
        # if finding anything but expected, fill in actual signal
        if datatype < 4:
            for i in range(num_fends - 1):
                fend1 = mapping[i]
                if filter[fend1] == 0:
                    continue
                j = i + 1
                k = indices[fend1 - startfend]
                last_index = indices[fend1 - startfend + 1]
                while j < max_fend[i]:
                    fend2 = mapping[j]
                    if filter[fend2] > 0:
                        while k < last_index and data[k, 1] < fend2:
                            k += 1
                        if k < last_index and data[k, 1] == fend2:
                            signal[i, j - i - 1, 0] = data[k, 2]
                    j += 1
        # fill in expected signal
        for i in range(num_fends - 1):
            fend1 = mapping[i]
            if filter[fend1] == 0:
                continue
            j = i + 1
            fend2 = mapping[j]
            k = 0
            # find opposite strand adjacents, skipping same fragment and same strand adjacents
            while j < num_fends and fend2 < (fend1 / 2 + 2) * 2:
                if fend2 == fend1 + 2 and filter[fend2] == 1:
                     # give starting expected value
                    signal[i, j - i - 1, 1] = 1.0
                    # if finding fend, enrichment, or expected, correct for fend
                    if datatype > 0 and datatype != 2:
                        signal[i, j - i - 1, 1] *= corrections[fend1] * corrections[fend2]
                    # if finding distance, enrichment, or expected, correct for distance
                    if datatype > 1:
                        distance = mids[j] - mids[i]
                        while k < num_distance_bins and distance_mids[k] < distance:
                            k += 1
                        if k == 0:
                            signal[i, j - i - 1, 1] *= distance_means[0]
                        elif k < num_distance_bins:
                            frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                            signal[i, j - i - 1, 1] *= distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                        else:
                            signal[i, j - i - 1, 1] *= distance_means[num_distance_bins - 1]
                j += 1
                if j < num_fends:
                    fend2 = mapping[j]
            while j < max_fend[i]:
                fend2 = mapping[j]
                if filter[fend2] > 0:
                    # give starting expected value
                    signal[i, j - i - 1, 1] = 1.0
                    # if finding fend, enrichment, or expected, correct for fend
                    if datatype > 0 and datatype != 2:
                        signal[i, j - i - 1, 1] *= corrections[fend1] * corrections[fend2]
                    # if finding distance, enrichment, or expected, correct for distance
                    if datatype > 1:
                        distance = mids[j] - mids[i]
                        while k < num_distance_bins and distance_mids[k] < distance:
                            k += 1
                        if k == 0:
                            signal[i, j - i - 1, 1] *= distance_means[0]
                        elif k < num_distance_bins:
                            frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                            signal[i, j - i - 1, 1] *= distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                        else:
                            signal[i, j - i - 1, 1] *= distance_means[num_distance_bins - 1]
                    # if finding expected only, fill in filter values for observed signal
                    if datatype == 4:
                        signal[i, j - i - 1, 0] = 1.0
                j += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def unbinned_signal_upper(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_means not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        int datatype,
        int startfend):
    cdef long long int fend1, fend2, k, distance
    cdef long long int index, i, j
    cdef double frac
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int num_distance_bins = distance_mids.shape[0]
    cdef long long int signal_size = signal.shape[0]
    with nogil:
        # if finding anything but expected, fill in actual signal
        if datatype < 4:
            for i in range(num_fends - 1):
                fend1 = mapping[i]
                if filter[fend1] == 0:
                    continue
                index = i * num_fends - i * (i + 1) / 2 - i - 1
                j = i + 1
                for k in range(indices[fend1 - startfend], indices[fend1 - startfend + 1]):
                    fend2 = data[k, 1]
                    if filter[fend2] == 0:
                        continue
                    while j < num_fends and mapping[j] < fend2:
                        j += 1
                    if j < num_fends and mapping[j] == fend2:
                        signal[index + j, 0] = data[k, 2]
        # fill in expected signal
        for i in range(num_fends - 1):
            fend1 = mapping[i]
            if filter[fend1] == 0:
                continue
            index = i * num_fends - i * (i + 1) / 2 - i - 1
            j = i + 1
            fend2 = mapping[j]
            k = 0
            # find opposite strand adjacents, skipping same fragment and same strand adjacents
            while j < num_fends and fend2 < (fend1 / 2 + 2) * 2:
                if fend2 == fend1 + 2 and filter[fend2] == 1:
                     # give starting expected value
                    signal[index + j, 1] = 1.0
                    # if finding fend, enrichment, or expected, correct for fend
                    if datatype > 0 and datatype != 2:
                        signal[index + j, 1] *= corrections[fend1] * corrections[fend2]
                    # if finding distance, enrichment, or expected, correct for distance
                    if datatype > 1:
                        distance = mids[j] - mids[i]
                        while k < num_distance_bins and distance_mids[k] < distance:
                            k += 1
                        if k == 0:
                            signal[index + j, 1] *= distance_means[0]
                        elif k < num_distance_bins:
                            frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                            signal[index + j, 1] *= distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                        else:
                            signal[index + j, 1] *= distance_means[num_distance_bins - 1]
                j += 1
                if j < num_fends:
                    fend2 = mapping[j]
            while j < max_fend[i]:
                fend2 = mapping[j]
                if filter[fend2] > 0:
                    # give starting expected value
                    signal[index + j, 1] = 1.0
                    # if finding fend, enrichment, or expected, correct for fend
                    if datatype > 0 and datatype != 2:
                        signal[index + j, 1] *= corrections[fend1] * corrections[fend2]
                    # if finding distance, enrichment, or expected, correct for distance
                    if datatype > 1:
                        distance = mids[j] - mids[i]
                        while k < num_distance_bins and distance_mids[k] < distance:
                            k += 1
                        if k == 0:
                            signal[index + j, 1] *= distance_means[0]
                        elif k < num_distance_bins:
                            frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                            signal[index + j, 1] *= distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                        else:
                            signal[index + j, 1] *= distance_means[num_distance_bins - 1]
                    # if finding expected only, fill in filter values for observed signal
                    if datatype == 4:
                        signal[index + j, 0] = 1.0
                j += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binned_signal_compact(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_means not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        int datatype):
    cdef long long int i, j, k, l, distance
    cdef double frac, expected, observed
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int num_distance_bins = distance_mids.shape[0]
    with nogil:
        for i in range(num_fends - 1):
            if filter[i] == 0 or mapping[i] == -1:
                continue
            k = 0
            l = indices[i]
            # fill in signal
            for j in range(i + 2, max_fend[i]):
                if filter[j] == 0 or mapping[j] == -1 or mapping[i] == mapping[j] or (i / 2 == j / 2 and i != j + 2):
                    continue
                # if requested, find observed signal
                if datatype < 4:
                    while l < indices[i + 1] and data[l, 1] < j:
                        l += 1
                    if l < indices[i + 1] and data[l, 1] == j:
                        observed = data[l, 2]
                    else:
                        observed = 0.0
                else:
                    observed = 1.0
                # find expected signal
                expected = 1.0
                # if finding fend, enrichment, or expected, correct for fend
                if datatype > 0 and datatype != 2:
                    expected *= corrections[i] * corrections[j]
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    distance = mids[j] - mids[i]
                    while k < num_distance_bins and distance_mids[k] < distance:
                        k += 1
                    if k == 0:
                        expected *= distance_means[0]
                    elif k < num_distance_bins:
                        frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                        expected *= distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                    else:
                        expected *= distance_means[num_distance_bins - 1]
                signal[mapping[i], mapping[j] - mapping[i] - 1, 0] += observed
                signal[mapping[i], mapping[j] - mapping[i] - 1, 1] += expected
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binned_signal_upper(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_means not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        int datatype,
        int num_bins):
    cdef long long int i, j, k, l, distance, index
    cdef double frac, expected, observed
    cdef long long int num_fends = mapping.shape[0]
    cdef long long int num_distance_bins = distance_mids.shape[0]
    with nogil:
        for i in range(num_fends - 1):
            if filter[i] == 0 or mapping[i] == -1:
                continue
            index = mapping[i] * num_bins - mapping[i] * (mapping[i] + 1) / 2 - mapping[i] - 1
            k = 0
            l = indices[i]
            # fill in signal
            for j in range(i + 2, max_fend[i]):
                if filter[j] == 0 or mapping[j] == -1 or mapping[i] == mapping[j] or (i / 2 == j / 2 and i != j + 2):
                    continue
                # if requested, find observed signal
                if datatype < 4:
                    while l < indices[i + 1] and data[l, 1] < j:
                        l += 1
                    if l < indices[i + 1] and data[l, 1] == j:
                        observed = data[l, 2]
                    else:
                        observed = 0.0
                else:
                    observed = 1.0
                # find expected signal
                expected = 1.0
                # if finding fend, enrichment, or expected, correct for fend
                if datatype > 0 and datatype != 2:
                    expected *= corrections[i] * corrections[j]
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    distance = mids[j] - mids[i]
                    while k < num_distance_bins and distance_mids[k] < distance:
                        k += 1
                    if k == 0:
                        expected *= distance_means[0]
                    elif k < num_distance_bins:
                        frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                        expected *= distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                    else:
                        expected *= distance_means[num_distance_bins - 1]
                signal[index + mapping[j], 0] += observed
                signal[index + mapping[j], 1] += expected
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
                if binned[x, y - x - 1, 0] < minobservations and removefailed:
                    binned[x, y - x - 1, 0] = 0
                    binned[x, y - x - 1, 1] = 0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binned_signal_trans(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double trans_mean,
        int start1,
        int start2,
        int datatype):
    cdef long long int i, j, k, bin1, bin2, fend1, fend2
    cdef double expected, observed
    cdef long long int num_fends1 = mapping1.shape[0]
    cdef long long int num_fends2 = mapping2.shape[0]
    cdef long long int num_bins1 = data.shape[0]
    cdef long long int num_bins2 = data.shape[1]
    with nogil:
        for i in range(num_fends1):
            fend1 = i + start1
            if filter[fend1] == 0:
                continue
            bin1 = mapping1[i]
            k = indices[i]
            if bin1 == -1:
                continue
            # fill in signal
            for j in range(num_fends2):
                fend2 = j + start2
                if filter[fend2] == 0:
                    continue
                bin2 = mapping2[j]
                if bin2 == -1:
                    continue
                # if requested, find observed signal
                if datatype < 4:
                    while k < indices[i + 1] and data[k, 1] < fend2:
                        k += 1
                    if k < indices[i + 1] and data[k, 1] == fend2:
                        observed = data[k, 2]
                    else:
                        observed = 0.0
                else:
                    observed = 1.0
                # find expected signal
                expected = 1.0
                # if finding fend, enrichment, or expected, correct for fend
                if datatype > 0 and datatype != 2:
                    expected *= corrections[fend1] * corrections[fend2]
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    expected *= trans_mean
                signal[bin1, bin2, 0] += observed
                signal[bin1, bin2, 1] += expected
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
        int maxsearch):
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
                if binned[x, y, 0] < minobservations:
                    binned[x, y, 0] = 0
                    binned[x, y, 1] = 0
    return None
