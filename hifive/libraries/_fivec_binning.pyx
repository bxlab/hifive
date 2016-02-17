# distutils: language = c++

"""These functions provide increased speed in handling the signal-binning functions necessary for supporting fivec analysis using HiFive."""

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
def find_cis_compact_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None):
    cdef long long int frag1, frag2, i, map1, map2
    cdef long long int num_frags = mapping.shape[0]
    with nogil:
        for frag1 in range(num_frags - 1):
            map1 = mapping[frag1]
            if map1 == 0:
                continue
            if map1 < 0:
                map1 = -1 - map1
            else:
                map1 = map1 - 1
            for i in range(indices[frag1], indices[frag1 + 1]):
                frag2 = data[i, 1]
                if frag2 >= num_frags:
                    continue
                map2 = mapping[frag2]
                if map2 == 0:
                    continue
                if map2 < 0:
                    signal[map1, -1 - map2, 0] += data[i, 2]
                else:
                    signal[map2 - 1, map1, 0] += data[i, 2]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_cis_upper_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None):
    cdef long long int frag1, frag2, i, index, map1, map2
    cdef long long int num_frags = mapping.shape[0]
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5))
    with nogil:
        for frag1 in range(num_frags - 1):
            map1 = mapping[frag1]
            if map1 == -1:
                continue
            index = map1 * (num_bins - 1) - map1 * (map1 + 1) / 2 - 1
            for i in range(indices[frag1], indices[frag1 + 1]):
                frag2 = data[i, 1]
                if frag2 >= num_frags:
                    continue
                map2 = mapping[frag2]
                if map2 == -1 or map2 == map1:
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
        np.ndarray[DTYPE_int_t, ndim=2] frag_indices,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_int_t, ndim=1] strands,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double gamma,
        double region_mean,
        int startfrag):
    cdef long long int frag1, frag2, map1, map2, strand1, bin1, bin2, index, num_parameters
    cdef double distance, value
    cdef long long int num_frags = mapping.shape[0]
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5))
    if not frag_indices is None:
        num_parameters = frag_indices.shape[1]
    else:
        num_parameters = 0
    with nogil:
        for frag1 in range(num_frags - 1):
            map1 = mapping[frag1]
            if map1 == 0:
                continue
            if map1 < 0:
                map1 = -1 - map1
            else:
                map1 = map1 - 1
            strand1 = strands[frag1]
            for frag2 in range(frag1 + 1, num_frags):
                map2 = mapping[frag2]
                if strands[frag2] == strand1 or map2 == 0:
                    continue
                 # give starting expected value
                value = 0.0
                # if finding frag, enrichment, or expected, and using express or probability bias correction, correct for frag
                if not corrections is None:
                    value += corrections[frag1] + corrections[frag2]
                # if finding frag, enrichment, or expected, and using binning bias correction, correct for frag
                if not binning_corrections is None:
                    for j in range(num_parameters):
                        bin1 = min(frag_indices[frag1 + startfrag, j], frag_indices[frag2 + startfrag, j]) 
                        bin2 = max(frag_indices[frag1 + startfrag, j], frag_indices[frag2 + startfrag, j]) 
                        index = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value += binning_corrections[index]
                # if finding distance, enrichment, or expected, correct for distance
                if gamma != 0.0:
                    distance = log(<double>(mids[frag2] - mids[frag1]))
                    value += region_mean - distance * gamma
                if strand1 == 0:
                    signal[map1, -1 - map2, 1] += exp(value)
                else:
                    signal[map2 - 1, map1, 1] += exp(value)
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
        np.ndarray[DTYPE_int_t, ndim=2] frag_indices,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] strands not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        double gamma,
        double region_mean,
        int startfrag):
    cdef long long int frag1, frag2, index, map1, map2, strand1, bin1, bin2, index2, num_parameters
    cdef double distance, value
    cdef long long int num_frags = mapping.shape[0]
    cdef long long int num_bins = int(0.5 + pow(0.25 + 2 * signal.shape[0], 0.5))
    if not frag_indices is None:
        num_parameters = frag_indices.shape[1]
    else:
        num_parameters = 0
    with nogil:
        for frag1 in range(num_frags - 1):
            map1 = mapping[frag1]
            if map1 == -1:
                continue
            strand1 = strands[frag1]
            index = map1 * (num_bins - 1) - map1 * (map1 + 1) / 2 - 1
            for frag2 in range(frag1 + 1, num_frags):
                map2 = mapping[frag2]
                if strands[frag2] == strand1 or map2 == -1 or map2 == map1:
                    continue
                 # give starting expected value
                value = 0.0
                # if finding frag, enrichment, or expected, and using express or probability bias correction, correct for frag
                if not corrections is None:
                    value += corrections[frag1] + corrections[frag2]
                # if finding frag, enrichment, or expected, and using binning bias correction, correct for frag
                if not binning_corrections is None:
                    for j in range(num_parameters):
                        bin1 = min(frag_indices[frag1 + startfrag, j], frag_indices[frag2 + startfrag, j]) 
                        bin2 = max(frag_indices[frag1 + startfrag, j], frag_indices[frag2 + startfrag, j]) 
                        index2 = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value += binning_corrections[index2]
                # if finding distance, enrichment, or expected, correct for distance
                if gamma != 0.0:
                    distance = log(<double>(mids[frag2] - mids[frag1]))
                    value += region_mean - distance * gamma
                signal[index + map2, 1] += exp(value)
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
def find_trans_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int64_t, ndim=1] indices,
        np.ndarray[DTYPE_int_t, ndim=1] mapping1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] mapping2 not None,
        np.ndarray[DTYPE_t, ndim=3] signal not None):
    cdef long long int frag1, frag2, i, map1, map2
    cdef long long int num_frags1 = mapping1.shape[0]
    cdef long long int num_frags2 = mapping2.shape[0]
    cdef long long int num_data = data.shape[0]
    with nogil:
        for frag1 in range(num_frags1):
            map1 = mapping1[frag1]
            if map1 == -1:
                continue
            for i in range(indices[frag1], indices[frag1 + 1]):
                frag2 = data[i, 1]
                if frag2 < 0 or frag2 >= num_frags2:
                    continue
                map2 = mapping2[frag2]
                if map2 == -1:
                    continue
                signal[map1, map2, 0] += data[i, 2]
        """
        for i in range(num_data):
            frag1 = data[i, 0]
            frag2 = data[i, 1]
            if frag2 < 0 or frag2 >= num_frags2:
                continue
            map1 = mapping1[frag1]
            map2 = mapping2[frag2]
            if map1 == -1 or map2 == -1:
                continue
            signal[map1, map2, 0] += data[i, 2]
        """
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
        np.ndarray[DTYPE_int_t, ndim=2] frag_indices,
        np.ndarray[DTYPE_int_t, ndim=1] strands1,
        np.ndarray[DTYPE_int_t, ndim=1] strands2,
        np.ndarray[DTYPE_t, ndim=3] signal not None,
        double trans_mean,
        int startfrag1,
        int startfrag2):
    cdef long long int frag1, frag2, map1, map2, strand1, bin1, bin2, index, num_parameters
    cdef double distance, value
    cdef long long int num_frags1 = mapping1.shape[0]
    cdef long long int num_frags2 = mapping2.shape[0]
    if not frag_indices is None:
        num_parameters = frag_indices.shape[1]
    else:
        num_parameters = 0
    with nogil:
        for frag1 in range(num_frags1):
            map1 = mapping1[frag1]
            if map1 == -1:
                continue
            strand1 = strands1[frag1]
            for frag2 in range(num_frags2):
                map2 = mapping2[frag2]
                if strands2[frag2] == strand1 or map2 == -1:
                    continue
                 # give starting expected value
                value = trans_mean
                # if finding frag, enrichment, or expected, and using express or probability bias correction, correct for frag
                if not corrections1 is None:
                    value += corrections1[frag1] + corrections2[frag2]
                # if finding frag, enrichment, or expected, and using binning bias correction, correct for frag
                if not binning_corrections is None:
                    for j in range(num_parameters):
                        bin1 = min(frag_indices[frag1 + startfrag1, j], frag_indices[frag2 + startfrag2, j]) 
                        bin2 = max(frag_indices[frag1 + startfrag1, j], frag_indices[frag2 + startfrag2, j]) 
                        index = bin1 * (binning_num_bins[j] - 1) - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                        value += binning_corrections[index]
                signal[map1, map2, 1] += exp(value)
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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binning_bin_observed(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=2] trans_data,
        np.ndarray[DTYPE_int_t, ndim=1] regions,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_int64_t, ndim=1] counts,
        np.ndarray[DTYPE_64_t, ndim=2] sums,
        np.ndarray[DTYPE_int_t, ndim=2] all_indices,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=1] bin_divs,
        np.ndarray[DTYPE_t, ndim=1] region_means,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        int mindistance,
        int maxdistance,
        float gamma,
        float trans_mean):
    cdef long long int i, k, distance, frag1, frag2, index, bin1, bin2, prev_frag, num_data, num_trans
    cdef double log_dist, signal, log_count
    cdef long long int num_features = all_indices.shape[1]
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    if not trans_data is None:
        num_trans = trans_data.shape[0]
    else:
        num_trans = 0
    with nogil:
        prev_frag = -1
        for i in range(num_data):
            frag1 = data[i, 0]
            if frag1 != prev_frag:
                prev_frag = frag1
                k = 0
            frag2 = data[i, 1]
            distance = mids[frag2] - mids[frag1]
            if distance < mindistance or distance > maxdistance:
                continue
            index = 0
            for j in range(num_features):
                bin1 = min(all_indices[frag1, j], all_indices[frag2, j])
                bin2 = max(all_indices[frag1, j], all_indices[frag2, j])
                index += ((num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2) * bin_divs[j]
            counts[index] += 1
            log_dist = log(<double>distance)
            signal = region_means[regions[frag1]] - log_dist * gamma
            log_count = log(<double>data[i, 2])
            if not corrections is None:
                log_count -= corrections[frag1] - corrections[frag2]
            sums[index, 0] += log_count - signal
            sums[index, 1] += pow(log_count - signal, 2.0)
        for i in range(num_trans):
            frag1 = trans_data[i, 0]
            frag2 = trans_data[i, 1]
            index = 0
            for j in range(num_features):
                bin1 = min(all_indices[frag1, j], all_indices[frag2, j])
                bin2 = max(all_indices[frag1, j], all_indices[frag2, j])
                index += ((num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2) * bin_divs[j]
            counts[index] += 1
            log_count = log(<double>trans_data[i, 2])
            if not corrections is None:
                log_count -= corrections[frag1] - corrections[frag2]
            sums[index, 0] += log_count - trans_mean
            sums[index, 1] += pow(log_count - trans_mean, 2.0)
    return None
