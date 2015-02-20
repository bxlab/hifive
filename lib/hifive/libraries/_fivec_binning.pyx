# distutils: language = c++

"""These functions provide increased speed in handling the signal-binning
functions necessary for supporting fivec analysis using HiFive.
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
def find_fragment_coverage(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] data_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] starts not None,
        np.ndarray[DTYPE_int_t, ndim=1] stops not None,
        np.ndarray[DTYPE_int_t, ndim=1] coverage not None,
        int mincoverage):
    cdef long long int i, j, frag1, valid
    cdef long long int num_regions = starts.shape[0]
    cdef long long int num_frags = filter.shape[0]
    with nogil:
        valid = 0
        for i in range(num_regions):
            for frag1 in range(starts[i], stops[i] - 1):
                if filter[frag1] == 0:
                    continue
                for j in range(data_indices[frag1], data_indices[frag1 + 1]):
                    if filter[data[j, 1]] == 0:
                        continue
                    coverage[frag1] += 1
                    coverage[data[j, 1]] += 1
        for i in range(num_frags):
            if coverage[i] < mincoverage:
                filter[i] = 0
            elif filter[i] > 0:
                valid += 1
    return valid


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def unbinned_signal_compact(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double mu,
        double gamma,
        double sigma,
        int datatype):
    cdef long long int frag1, frag2, i, j, distance, num_data
    cdef long long int num_frags = mapping.shape[0]
    cdef long long int xdim = signal.shape[0]
    cdef long long int ydim = signal.shape[1]
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    with nogil:
        # if finding anything but expected, fill in actual signal
        if datatype < 4:
            for i in range(num_data):
                if data[i, 1] >= num_frags or mapping[data[i, 0]] == 0 or mapping[data[i, 1]] == 0:
                    continue
                if filter[data[i, 0]] == 0 or filter[data[i, 1]] == 0:
                    continue
                if mapping[data[i, 0]] > 0:
                    frag1 = mapping[data[i, 0]] - 1
                    frag2 = -mapping[data[i, 1]] - 1
                else:
                    frag1 = mapping[data[i, 1]] - 1
                    frag2 = -mapping[data[i, 0]] - 1
                signal[frag1, frag2, 0] = data[i, 2]
        # fill in expected signal
        frag1 = 0
        for i in range(xdim):
            while frag1 < num_frags - 1 and mapping[frag1] - 1 != i:
                frag1 += 1
            if filter[frag1] == 0:
                continue
            frag2 = 0
            for j in range(ydim):
                while frag2 < num_frags - 1 and -mapping[frag2] - 1 != j:
                    frag2 += 1
                if filter[frag2] == 0:
                    continue
                # give starting expected value
                if datatype > 1:
                    signal[i, j, 1] = exp(mu)
                else:
                    signal[i, j, 1] = 1.0
                # if finding fragment, enrichment, or expected, correct for fragment
                if datatype > 0 and datatype != 2:
                    signal[i, j, 1] *= exp(corrections[frag1] + corrections[frag2])
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    if frag1 < frag2:
                        distance = mids[frag2] - mids[frag1]
                    else:
                        distance = mids[frag1] - mids[frag2]
                    signal[i, j, 1] *= exp(-gamma * log(distance))
                # if finding expected only, fill in filter values for observed signal
                if datatype == 4:
                    signal[i, j, 0] = 1.0
                # if finding enrichment, scale both observed and expected values by sigma
                if datatype == 3:
                    signal[i, j, 0] = pow(signal[i, j, 0], 1.0 / sigma)
                    signal[i, j, 1] = pow(signal[i, j, 1], 1.0 / sigma)
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def unbinned_signal_upper(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] strands not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        double mu,
        double gamma,
        double sigma,
        int num_bins,
        int datatype):
    cdef long long int frag1, frag2, i, j, k, index, distance, num_data
    cdef long long int num_frags = mapping.shape[0]
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    with nogil:
        # if finding anything but expected, fill in actual signal
        if datatype < 4:
            for frag1 in range(num_frags - 1):
                if filter[frag1] == 0:
                    continue
                i = mapping[frag1]
                index = i * num_bins - i * (i + 1) / 2 - i - 1
                for k in range(indices[frag1], indices[frag1 + 1]):
                    frag2 = data[k, 1]
                    if frag2 >= num_frags or filter[frag2] == 0:
                        continue
                    j = mapping[frag2]
                    signal[index + j, 0] = data[k, 2]
        # fill in expected signal
        for frag1 in range(num_frags - 1):
            if filter[frag1] == 0:
                continue
            i = mapping[frag1]
            index = i * num_bins - i * (i + 1) / 2 - i - 1
            for frag2 in range(frag1 + 1, num_frags):
                if filter[frag2] == 0 or strands[frag1] == strands[frag2]:
                    continue
                j = mapping[frag2]
                # give starting expected value
                if datatype > 1:
                    signal[index + j, 1] = exp(mu)
                else:
                    signal[index + j, 1] = 1.0
                # if finding fragment, enrichment, or expected, correct for fragment
                if datatype > 0 and datatype != 2:
                    signal[index + j, 1] *= exp(corrections[frag1] + corrections[frag2])
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    distance = mids[frag2] - mids[frag1]
                    signal[index + j, 1] *= exp(-gamma * log(distance))
                # if finding expected only, fill in filter values for observed signal
                if datatype == 4:
                    signal[index + j, 0] = 1.0
                # if finding enrichment, scale both observed and expected values by sigma
                if datatype == 3:
                    signal[index + j, 0] = pow(signal[index + j, 0], 1.0 / sigma)
                    signal[index + j, 1] = pow(signal[index + j, 1], 1.0 / sigma)
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binned_signal_upper(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] strands not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        double mu,
        double gamma,
        double sigma,
        int datatype,
        int num_bins):
    cdef long long int frag1, frag2, i, j, k, index, distance, num_data
    cdef double observed, expected
    cdef long long int num_frags = mapping.shape[0]
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    with nogil:
        # if finding anything but expected, fill in actual signal
        if datatype < 4:
            for frag1 in range(num_frags - 1):
                if filter[frag1] == 0:
                    continue
                i = mapping[frag1]
                index = i * num_bins - i * (i + 1) / 2 - i - 1
                for k in range(indices[frag1], indices[frag1 + 1]):
                    frag2 = data[k, 1]
                    if frag2 >= num_frags or filter[frag2] == 0 or strands[frag1] == strands[frag2]:
                        continue
                    j = mapping[frag2]
                    if i == j:
                        continue
                    observed = data[k, 2]
                    if datatype == 3:
                        observed = pow(observed, 1.0 / sigma)
                    signal[index + j, 0] += observed
        # fill in expected signal
        for frag1 in range(num_frags - 1):
            if filter[frag1] == 0:
                continue
            i = mapping[frag1]
            index = i * num_bins - i * (i + 1) / 2 - i - 1
            for frag2 in range(frag1 + 1, num_frags):
                j = mapping[frag2]
                if filter[frag2] == 0 or strands[frag1] == strands[frag2] or i == j:
                    continue
                # give starting expected value
                if datatype > 1:
                    expected = exp(mu)
                else:
                    expected = 1.0
                # if finding fragment, enrichment, or expected, correct for fragment
                if datatype > 0 and datatype != 2:
                    expected *= exp(corrections[frag1] + corrections[frag2])
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    distance = mids[frag2] - mids[frag1]
                    expected *= exp(-gamma * log(distance))
                # if finding expected only, fill in filter values for observed signal
                if datatype == 4:
                    signal[index + j, 0] += 1.0
                # if finding enrichment, scale both observed and expected values by sigma
                if datatype == 3:
                    expected = pow(expected, 1.0 / sigma)
                signal[index + j, 1] += expected
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bin_upper_to_upper(
        np.ndarray[DTYPE_t, ndim=2] binned not None,
        np.ndarray[DTYPE_t, ndim=2] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
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
                        index2 = lX * num_fends - lX * (lX + 1) / 2 - lX - 1
                        for i in range(max(lX + 1, lY), uY + 1):
                            binned[index + y, 0] += unbinned[index2 + i, 0]
                            binned[index + y, 1] += unbinned[index2 + i, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_fends - 1:
                            uX_dist = ub_mids[uX + 1] - b_mids[x]
                        else:
                            uX_dist = 1000000000
                        index2 = uX * num_fends - uX * (uX + 1) / 2 - uX - 1
                        for i in range(max(uX + 1, lY), uY + 1):
                            binned[index + y, 0] += unbinned[index2 + i, 0]
                            binned[index + y, 1] += unbinned[index2 + i, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids[y] - ub_mids[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(lX, min(uX + 1, lY)):
                            index2 = i * num_fends - i * (i + 1) / 2 - i - 1
                            binned[index + y, 0] += unbinned[index2 + lY, 0]
                            binned[index + y, 1] += unbinned[index2 + lY, 1]
                    if min_dist == uY_dist:
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
def unbinned_signal_trans_full(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] strands not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double trans_mu,
        double sigma,
        int startfrag1,
        int startfrag2,
        int datatype):
    cdef long long int frag1, frag2, i, j
    cdef long long int num_frags1 = mapping1.shape[0]
    cdef long long int num_frags2 = mapping2.shape[0]
    cdef long long int num_data = data.shape[0]
    with nogil:
        # if finding anything but expected, fill in actual signal
        for i in range(num_data):
            frag1 = data[i, 0] - startfrag1
            frag2 = data[i, 1] - startfrag2
            if frag2 >= num_frags2 or frag2 < 0 or filter[frag1 + startfrag1] == 0 or filter[frag2 + startfrag2] == 0:
                continue
            signal[mapping1[frag1], mapping2[frag2], 0] = data[i, 2]
        # fill in expected signal
        for frag1 in range(num_frags1):
            if filter[frag1 + startfrag1] == 0:
                continue
            i = mapping1[frag1]
            for frag2 in range(num_frags2):
                if filter[frag2 + startfrag2] == 0 or strands[frag1 + startfrag1] == strands[frag2 + startfrag2]:
                    continue
                j = mapping2[frag2]
                # give starting expected value
                signal[i, j, 1] = 1.0
                # if finding fragment, enrichment, or expected, correct for fragment
                if datatype > 0 and datatype != 2:
                    signal[i, j, 1] *= exp(corrections[frag1 + startfrag1] + corrections[frag2 + startfrag2])
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    signal[i, j, 1] *= trans_mu
                # if finding expected only, fill in filter values for observed signal
                if datatype == 4:
                    signal[i, j, 0] = 1.0
                # if finding enrichment, scale both observed and expected values by sigma
                if datatype == 3:
                    signal[i, j, 0] = pow(signal[i, j, 0], 1.0 / sigma)
                    signal[i, j, 1] = pow(signal[i, j, 1], 1.0 / sigma)
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def unbinned_signal_trans_compact(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] strands not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=3] signal1 not None,
        np.ndarray[DTYPE_t, ndim=3] signal2 not None,
        double trans_mu,
        double sigma,
        int startfrag1,
        int startfrag2,
        int datatype):
    cdef long long int frag1, frag2, i, j
    cdef long long int num_frags1 = mapping1.shape[0]
    cdef long long int num_frags2 = mapping2.shape[0]
    cdef long long int num_data = data.shape[0]
    with nogil:
        # if finding anything but expected, fill in actual signal
        for i in range(num_data):
            frag1 = data[i, 0] - startfrag1
            frag2 = data[i, 1] - startfrag2
            if frag2 >= num_frags2 or frag2 < 0 or filter[frag1 + startfrag1] == 0 or filter[frag2 + startfrag2] == 0:
                continue
            if mapping1[frag1] > 0:
                signal1[mapping1[frag1] - 1, -mapping2[frag2] - 1, 0] = data[i, 2]
            else:
                signal2[mapping2[frag2] - 1, -mapping1[frag1] - 1, 0] = data[i, 2]
        # fill in expected signal
        for frag1 in range(num_frags1):
            if filter[frag1 + startfrag1] == 0:
                continue
            if mapping1[frag1] > 0:
                i = mapping1[frag1] - 1
            else:
                i = -mapping1[frag1] - 1
            for frag2 in range(num_frags2):
                if filter[frag2 + startfrag2] == 0 or strands[frag1 + startfrag1] == strands[frag2 + startfrag2]:
                    continue
                if mapping2[frag2] > 0:
                    j = mapping2[frag2] - 1
                    # give starting expected value
                    signal2[j, i, 1] = 1.0
                    # if finding fragment, enrichment, or expected, correct for fragment
                    if datatype > 0 and datatype != 2:
                        signal2[j, i, 1] *= exp(corrections[frag1 + startfrag1] + corrections[frag2 + startfrag2])
                    # if finding distance, enrichment, or expected, correct for distance
                    if datatype > 1:
                        signal2[j, i, 1] *= trans_mu
                    # if finding expected only, fill in filter values for observed signal
                    if datatype == 4:
                        signal2[j, i, 0] = 1.0
                    # if finding enrichment, scale both observed and expected values by sigma
                    if datatype == 3:
                        signal2[j, i, 0] = pow(signal2[j, i, 0], 1.0 / sigma)
                        signal2[j, i, 1] = pow(signal2[j, i, 1], 1.0 / sigma)
                else:
                    j = -mapping2[frag2] - 1
                    # give starting expected value
                    signal1[i, j, 1] = 1.0
                    # if finding fragment, enrichment, or expected, correct for fragment
                    if datatype > 0 and datatype != 2:
                        signal1[i, j, 1] *= exp(corrections[frag1 + startfrag1] + corrections[frag2 + startfrag2])
                    # if finding distance, enrichment, or expected, correct for distance
                    if datatype > 1:
                        signal1[i, j, 1] *= trans_mu
                    # if finding expected only, fill in filter values for observed signal
                    if datatype == 4:
                        signal1[i, j, 0] = 1.0
                    # if finding enrichment, scale both observed and expected values by sigma
                    if datatype == 3:
                        signal1[i, j, 0] = pow(signal1[i, j, 0], 1.0 / sigma)
                        signal1[i, j, 1] = pow(signal1[i, j, 1], 1.0 / sigma)
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binned_signal_trans(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] strands not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double trans_mu,
        double sigma,
        int startfrag1,
        int startfrag2,
        int datatype):
    cdef long long int frag1, frag2, i, j
    cdef double observed, expected
    cdef long long int num_frags1 = mapping1.shape[0]
    cdef long long int num_frags2 = mapping2.shape[0]
    cdef long long int num_data = data.shape[0]
    with nogil:
        # if finding anything but expected, fill in actual signal
        for i in range(num_data):
            frag1 = data[i, 0] - startfrag1
            frag2 = data[i, 1] - startfrag2
            if frag2 >= num_frags2 or frag2 < 0 or filter[frag1 + startfrag1] == 0 or filter[frag2 + startfrag2] == 0:
                continue
            observed = data[i, 2]
            if datatype == 3:
                observed = pow(observed, 1.0 / sigma)
            signal[mapping1[frag1], mapping2[frag2], 0] += sigma
        # fill in expected signal
        for frag1 in range(num_frags1):
            if filter[frag1 + startfrag1] == 0:
                continue
            i = mapping1[frag1]
            for frag2 in range(num_frags2):
                if filter[frag2 + startfrag2] == 0 or strands[frag1 + startfrag1] == strands[frag2 + startfrag2]:
                    continue
                j = mapping2[frag2]
                # give starting expected value
                expected = 1.0
                # if finding fragment, enrichment, or expected, correct for fragment
                if datatype > 0 and datatype != 2:
                    expected *= exp(corrections[frag1 + startfrag1] + corrections[frag2 + startfrag2])
                # if finding distance, enrichment, or expected, correct for distance
                if datatype > 1:
                    expected *= trans_mu
                # if finding expected only, fill in filter values for observed signal
                if datatype == 4:
                    expected = 1.0
                # if finding enrichment, scale both observed and expected values by sigma
                if datatype == 3:
                    expected = pow(expected, 1.0 / sigma)
                signal[i, j, 1] += expected
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin_trans(
        np.ndarray[DTYPE_t, ndim=3] unbinned not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids2 not None,
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
    cdef long long int num_frags1 = mids1.shape[0]
    cdef long long int num_frags2 = mids2.shape[0]
    with nogil:
        for x in range(num_bins1):
            for y in range(num_bins2):
                # if bin already meets our criteria, skip
                if binned[x, y, 0] >= minobservations:
                    continue
                # otherwise, set boarding unbinned positions according to bounds
                lX = bounds1[x, 0]
                uX = bounds1[x, 1] - 1
                lY = bounds2[y, 0]
                uY = bounds2[y, 1] - 1
                # find distance in each direction
                if lX > 0:
                    lX_dist = b_mids1[x] - mids1[lX - 1]
                else:
                    lX_dist = 1000000000
                if uX < num_frags1 - 1:
                    uX_dist = mids1[uX + 1] - b_mids1[x]
                else:
                    uX_dist = 1000000000
                if lY > 0:
                    lY_dist = b_mids2[y] - mids2[lY - 1]
                else:
                    lY_dist = 1000000000
                if uY < num_frags2 - 1:
                    uY_dist = mids2[uY + 1] - b_mids2[y]
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
                            lX_dist = b_mids1[x] - mids1[lX - 1]
                        else:
                            lX_dist = 1000000000
                        for i in range(lY, uY + 1):
                            binned[x, y, 0] += unbinned[lX, i, 0]
                            binned[x, y, 1] += unbinned[lX, i, 1]
                    if min_dist == uX_dist:
                        uX += 1
                        if uX < num_frags1 - 1:
                            uX_dist = mids1[uX + 1] - b_mids1[x]
                        else:
                            uX_dist = 1000000000
                        for i in range(lY, uY + 1):
                            binned[x, y, 0] += unbinned[uX, i, 0]
                            binned[x, y, 1] += unbinned[uX, i, 1]
                    if min_dist == lY_dist:
                        lY -= 1
                        if lY > 0:
                            lY_dist = b_mids2[y] - mids2[lY - 1]
                        else:
                            lY_dist = 1000000000
                        for i in range(lX, uX + 1):
                            binned[x, y, 0] += unbinned[i, lY, 0]
                            binned[x, y, 1] += unbinned[i, lY, 1]
                    if min_dist == uY_dist:
                        uY += 1
                        if uY < num_frags2 - 1:
                            uY_dist = mids2[uY + 1] - b_mids2[y]
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
