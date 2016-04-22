# distutils: language = c++

"""These functions provide increased speed in handling the distance-
dependent signal functions necessary for supporting hic analysis using
HiFive.
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
def find_distance_bin_sums(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] rev_mapping not None,
        np.ndarray[DTYPE_t, ndim=1] cutoffs not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_64_t, ndim=1] counts not None,
        np.ndarray[DTYPE_int_t, ndim=2] indices not None,
        np.ndarray[DTYPE_int64_t, ndim=2] bin_size not None,
        np.ndarray[DTYPE_64_t, ndim=1] count_sum not None,
        np.ndarray[DTYPE_64_t, ndim=1] logdistance_sum not None,
        int start,
        int stop,
        int binned):
    cdef long long int i, j, temp, fend1, fend2, previous_fend
    cdef double log_dist
    cdef long long int num_data = indices.shape[0]
    cdef long long int num_fends = rev_mapping.shape[0]
    with nogil:
        previous_fend = -1
        j = 0
        for i in range(num_data):
            fend1 = mapping[indices[i, 0]]
            fend2 = mapping[indices[i, 1]]
            if fend1 == -1 or fend2 == -1:
                continue
            if fend1 != previous_fend:
                j = 0
                previous_fend = fend1
            log_dist = log(<double>(max(1, mids[fend2] - mids[fend1])))
            while log_dist > cutoffs[j]:
                j += 1
            count_sum[j] += counts[i]
            bin_size[j, 0] += 1
        for fend1 in range(start, stop):
            j = 0
            if binned == 0:
                for fend2 in range(fend1 + 1, min(fend1 + 4, num_fends)):
                    if rev_mapping[fend1] % 2 == 0:
                        temp = rev_mapping[fend2] - rev_mapping[fend1]
                        if temp == 1 or temp == 3:
                            continue
                    else:
                        if rev_mapping[fend2] - rev_mapping[fend1] == 1:
                            continue
                    log_dist = log(<double>(mids[fend2] - mids[fend1]))
                    while log_dist > cutoffs[j]:
                        j += 1
                    bin_size[j, 1] += 1
                    logdistance_sum[j] += log_dist
                for fend2 in range(min(fend1 + 4, num_fends), num_fends):
                    log_dist = log(<double>(mids[fend2] - mids[fend1]))
                    while log_dist > cutoffs[j]:
                        j += 1
                    bin_size[j, 1] += 1
                    logdistance_sum[j] += log_dist
            else:
                for fend2 in range(fend1, num_fends):
                    log_dist = log(<double>(max(1, mids[fend2] - mids[fend1])))
                    while log_dist > cutoffs[j]:
                        j += 1
                    bin_size[j, 1] += 1
                    logdistance_sum[j] += log_dist
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_binary_distance_bin_sums(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] rev_mapping not None,
        np.ndarray[DTYPE_t, ndim=1] cutoffs not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=2] indices not None,
        np.ndarray[DTYPE_int64_t, ndim=2] counts not None,
        np.ndarray[DTYPE_64_t, ndim=1] logdistance_sum not None,
        int start,
        int stop,
        int binned):
    cdef long long int i, j, temp, fend1, fend2, previous_fend
    cdef double log_dist
    cdef long long int num_data = indices.shape[0]
    cdef long long int num_fends = rev_mapping.shape[0]
    with nogil:
        previous_fend = -1
        j = 0
        for i in range(num_data):
            fend1 = mapping[indices[i, 0]]
            fend2 = mapping[indices[i, 1]]
            if fend1 == -1 or fend2 == -1:
                continue
            if fend1 != previous_fend:
                j = 0
                previous_fend = fend1
            log_dist = log(<double>(max(1, mids[fend2] - mids[fend1])))
            while log_dist > cutoffs[j]:
                j += 1
            counts[j, 0] += 1
        for fend1 in range(start, stop):
            j = 0
            if binned == 0:
                for fend2 in range(fend1 + 1, min(fend1 + 4, num_fends)):
                    if rev_mapping[fend1] % 2 == 0:
                        temp = rev_mapping[fend2] - rev_mapping[fend1]
                        if temp == 1 or temp == 3:
                            continue
                    else:
                        if rev_mapping[fend2] - rev_mapping[fend1] == 1:
                            continue
                    log_dist = log(<double>(mids[fend2] - mids[fend1]))
                    while log_dist > cutoffs[j]:
                        j += 1
                    counts[j, 1] += 1
                    logdistance_sum[j] += log_dist
                for fend2 in range(min(fend1 + 4, num_fends), num_fends):
                    log_dist = log(<double>(mids[fend2] - mids[fend1]))
                    while log_dist > cutoffs[j]:
                        j += 1
                    counts[j, 1] += 1
                    logdistance_sum[j] += log_dist
            else:
                for fend2 in range(fend1, num_fends):
                    log_dist = log(<double>(max(1, mids[fend2] - mids[fend1])))
                    while log_dist > cutoffs[j]:
                        j += 1
                    counts[j, 1] += 1
                    logdistance_sum[j] += log_dist
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_remapped_distance_means(
        np.ndarray[DTYPE_int_t, ndim=1] indices0,
        np.ndarray[DTYPE_int_t, ndim=1] indices1,
        np.ndarray[DTYPE_int_t, ndim=1] mids,
        np.ndarray[DTYPE_t, ndim=1] means,
        np.ndarray[DTYPE_t, ndim=2] parameters,
        float chrom_mean):
    cdef long long int i, j, index0, previous_index
    cdef double frac, distance
    cdef long long int num_pairs = indices0.shape[0]
    with nogil:
        previous_index = -1
        j = 0
        for i in range(num_pairs):
            index0 = indices0[i]
            distance = log(<double>(max(1, mids[indices1[i]] - mids[index0])))
            if index0 != previous_index:
                previous_index = index0
                j = 0
            while distance > parameters[j, 0]:
                j += 1
            means[i] = exp(parameters[j, 1] * distance + parameters[j, 2] + chrom_mean)
    return None
