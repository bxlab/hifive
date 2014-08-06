# distutils: language = c++
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

"""These functions provide support for handling BI calculation and analysis."""

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
def find_bi(
        np.ndarray[DTYPE_t, ndim=3] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] BI_mids not None,
        np.ndarray[DTYPE_t, ndim=1] BI not None,
        np.ndarray[DTYPE_t, ndim=2] temp not None,
        int width,
        int mincount):
    cdef int i, j, k, count, sum1, sum2, upbound, downbound
    cdef double diff
    cdef int num_bins = data.shape[0]
    cdef int max_bin = data.shape[1]
    with nogil:
        for i in range(1, num_bins):
            # find width around boundary point
            upbound = i - 1
            while upbound > 0 and BI_mids[i - 1] - mids[upbound - 1] < width:
                upbound -= 1
            downbound = i + 1
            while downbound < num_bins and mids[downbound] - BI_mids[i - 1] < width:
                downbound += 1
            # zero out temp array
            for j in range(2 * max_bin):
                temp[j, 0] = 0
                temp[j, 1] = 0
                temp[j, 2] = 0
                temp[j, 3] = 0
            # fill temp array with current interaction set
            # find upstream interactions
            for j in range(1, min(max_bin, i)):
                for k in range(max(upbound, i - j), i):
                    temp[max_bin - j, 0] += data[i - j - 1, k - i + j, 0]
                    temp[max_bin - j, 1] += data[i - j - 1, k - i + j, 1]
                for k in range(i, min(downbound, max_bin + i - j)):
                    temp[max_bin - j, 2] += data[i - j - 1, k - i + j, 0]
                    temp[max_bin - j, 3] += data[i - j - 1, k - i + j, 1]
            # find downstream interactions
            for j in range(1, min(max_bin, num_bins - i - 1)):
                for k in range(max(upbound, i + j - max_bin), i):
                    temp[max_bin + j - 1, 0] += data[k, i + j - k - 1, 0]
                    temp[max_bin + j - 1, 1] += data[k, i + j - k - 1, 1]
                for k in range(i, min(downbound, i + j)):
                    temp[max_bin + j - 1, 2] += data[k, i + j - k - 1, 0]
                    temp[max_bin + j - 1, 3] += data[k, i + j - k - 1, 1]
            # find set counts
            count = 0
            for j in range(max_bin * 2):
                if temp[j, 0] > 0 and temp[j, 2] > 0:
                    count += 1
                    temp[j, 0] = log(temp[j, 0] / temp[j, 1])
                    temp[j, 2] = log(temp[j, 2] / temp[j, 3])
                    temp[j, 1] = 1
                else:
                    temp[j, 1] = 0
            if count < mincount:
                BI[i - 1] = Inf
                continue
            # find difference
            for j in range(max_bin * 2):
                if temp[j, 1] == 1:
                    diff = temp[j, 0] - temp[j, 2]
                    if diff < 0:
                        BI[i - 1] += diff
                    else:
                        BI[i - 1] -= diff
            BI[i - 1] /= count
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_bi_height(
        np.ndarray[DTYPE_t, ndim=3] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] BI_mids not None,
        np.ndarray[DTYPE_t, ndim=1] BI not None,
        np.ndarray[DTYPE_t, ndim=2] temp not None,
        int width,
        int window,
        int height,
        int mincount):
    cdef int i, j, k, count, sum1, sum2, upbound, downbound, index
    cdef double diff
    cdef int num_bins = data.shape[0]
    cdef int max_bin = data.shape[1]
    cdef int temp_bins = temp.shape[0]
    cdef int half_temp = temp_bins / 2
    with nogil:
        for i in range(1, num_bins):
            # find width around boundary point
            upbound = i - 1
            while upbound > 0 and BI_mids[i - 1] - mids[upbound - 1] < width:
                upbound -= 1
            downbound = i + 1
            while downbound < num_bins and mids[downbound] - BI_mids[i - 1] < width:
                downbound += 1
            # zero out temp array
            for j in range(temp_bins):
                temp[j, 0] = 0
                temp[j, 1] = 0
                temp[j, 2] = 0
                temp[j, 3] = 0
            # fill temp array with current interaction set
            # find upstream interactions
            for j in range(1, min(max_bin, i)):
                if BI_mids[i - 1] - mids[i - j - 1] > window:
                    continue
                index = (BI_mids[i - 1] - mids[i - j - 1]) / height
                for k in range(max(upbound, i - j), i):
                    temp[half_temp - index - 1, 0] += data[i - j - 1, k - i + j, 0]
                    temp[half_temp - index - 1, 1] += data[i - j - 1, k - i + j, 1]
                for k in range(i, min(downbound, max_bin + i - j)):
                    temp[half_temp - index - 1, 2] += data[i - j - 1, k - i + j, 0]
                    temp[half_temp - index - 1, 3] += data[i - j - 1, k - i + j, 1]
            # find downstream interactions
            for j in range(1, min(max_bin, num_bins - i - 1)):
                if mids[i + j] - BI_mids[i - 1] > window:
                    continue
                index = (mids[i + j] - BI_mids[i - 1]) / height
                for k in range(max(upbound, i + j - max_bin), i):
                    temp[half_temp + index, 0] += data[k, i + j - k - 1, 0]
                    temp[half_temp + index, 1] += data[k, i + j - k - 1, 1]
                for k in range(i, min(downbound, i + j)):
                    temp[half_temp + index, 2] += data[k, i + j - k - 1, 0]
                    temp[half_temp + index, 3] += data[k, i + j - k - 1, 1]
            # find set counts
            count = 0
            for j in range(temp_bins):
                if temp[j, 0] > 0 and temp[j, 2] > 0:
                    count += 1
                    temp[j, 0] = log(temp[j, 0] / temp[j, 1])
                    temp[j, 2] = log(temp[j, 2] / temp[j, 3])
                    temp[j, 1] = 1
                else:
                    temp[j, 1] = 0
            if count < mincount:
                BI[i - 1] = Inf
                continue
            # find difference
            for j in range(temp_bins):
                if temp[j, 1] == 1:
                    diff = temp[j, 0] - temp[j, 2]
                    if diff < 0:
                        BI[i - 1] += diff
                    else:
                        BI[i - 1] -= diff
            BI[i - 1] /= count
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_bi_old(
        np.ndarray[DTYPE_t, ndim=3] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] BI_mids not None,
        np.ndarray[DTYPE_t, ndim=1] BI not None,
        np.ndarray[DTYPE_t, ndim=2] temp not None,
        int width,
        int mincount):
    cdef int i, j, k, count, sum1, sum2, upbound, downbound
    cdef double std1, std2, mean1, mean2, temp1
    cdef int num_bins = data.shape[0]
    cdef int max_bin = data.shape[1]
    with nogil:
        for i in range(1, num_bins):
            # find width around boundary point
            upbound = i - 1
            while upbound > 0 and BI_mids[i - 1] - mids[upbound - 1] < width:
                upbound -= 1
            downbound = i + 1
            while downbound < num_bins and mids[downbound] - BI_mids[i - 1] < width:
                downbound += 1
            # zero out temp array
            for j in range(2 * max_bin):
                temp[j, 0] = 0
                temp[j, 1] = 0
                temp[j, 2] = 0
                temp[j, 3] = 0
            # fill temp array with current interaction set
            # find upstream interactions
            for j in range(1, min(max_bin, i)):
                for k in range(max(upbound, i - j), i):
                    temp[max_bin - j, 0] += data[i - j - 1, k - i + j, 0]
                    temp[max_bin - j, 1] += data[i - j - 1, k - i + j, 1]
                for k in range(i, min(downbound, max_bin + i - j)):
                    temp[max_bin - j, 2] += data[i - j - 1, k - i + j, 0]
                    temp[max_bin - j, 3] += data[i - j - 1, k - i + j, 1]
            # find downstream interactions
            for j in range(1, min(max_bin, num_bins - i - 1)):
                for k in range(max(upbound, i + j - max_bin), i):
                    temp[max_bin + j - 1, 0] += data[k, i + j - k - 1, 0]
                    temp[max_bin + j - 1, 1] += data[k, i + j - k - 1, 1]
                for k in range(i, min(downbound, i + j)):
                    temp[max_bin + j - 1, 2] += data[k, i + j - k - 1, 0]
                    temp[max_bin + j - 1, 3] += data[k, i + j - k - 1, 1]
            # find set means
            mean1 = 0.0
            mean2 = 0.0
            count = 0
            for j in range(max_bin * 2):
                if temp[j, 0] > 0 and temp[j, 2] > 0:
                    count += 1
                    temp[j, 0] = log(temp[j, 0] / temp[j, 1])
                    temp[j, 2] = log(temp[j, 2] / temp[j, 3])
                    temp[j, 1] = 1
                    mean1 += temp[j, 0]
                    mean2 += temp[j, 2]
                else:
                    temp[j, 1] = 0
            if count < mincount:
                BI[i - 1] = Inf
                continue
            mean1 /= count
            mean2 /= count
            # find set stds
            std1 = 0.0
            std2 = 0.0
            for j in range(max_bin * 2):
                if temp[j, 1] == 1:
                    temp1 = temp[j, 0] - mean1
                    std1 += temp1 * temp1
                    temp1 = temp[j, 2] - mean2
                    std2 += temp1 * temp1
            std1 /= count - 1
            std2 /= count - 1
            # find correlation
            for j in range(max_bin * 2):
                if temp[j, 1] == 1:
                    BI[i - 1] += (temp[j, 0] - mean1) / std1 * (temp[j, 2] - mean2) / std2
            BI[i - 1] /= count
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_bi_profile(
        np.ndarray[DTYPE_int_t, ndim=1] starts not None,
        np.ndarray[DTYPE_int_t, ndim=1] stops not None,
        np.ndarray[DTYPE_t, ndim=1] scores not None,
        np.ndarray[DTYPE_int_t, ndim=1] chr_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] coordinates not None,
        np.ndarray[DTYPE_t, ndim=2] signal not None,
        int width,
        int chrint):
    cdef int i, j, k, l, start, stop, stop_index
    cdef double frac
    cdef int num_bins = signal.shape[0]
    cdef int num_coords = coordinates.shape[0]
    cdef double width_f = width
    with nogil:
        j = chr_indices[chrint]
        stop_index = chr_indices[chrint + 1]
        for i in range(num_coords):
            start = coordinates[i] - (width * num_bins) / 2
            while j < stop_index and stops[j] < start:
                j += 1
            l = j
            for k in range(num_bins):
                while l < stop_index and stops[l] < start + width:
                    frac = (min(start + width, stops[l]) - max(start, starts[l])) / width_f
                    signal[k, 0] += scores[l] * frac
                    signal[k, 1] += frac
                    l += 1
                if l < stop_index and starts[l] < start + width:
                    frac = (min(start + width, stops[l]) - max(start, starts[l])) / width_f
                    signal[k, 0] += scores[l] * frac
                    signal[k, 1] += frac
                start += width
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def gaussian_smoothing(
        np.ndarray[DTYPE_int_t, ndim=1] X not None,
        np.ndarray[DTYPE_t, ndim=1] Y not None,
        np.ndarray[DTYPE_t, ndim=1] smoothed not None,
        int width):
    cdef int i, j, start, stop
    cdef double total_weight, weight
    cdef int dim = X.shape[0]
    cdef double cutoff = width * 2.5
    with nogil:
        start = 0
        stop = 0
        for i in range(dim):
            total_weight = 0.0
            smoothed[i] = 0.0
            while start < i and X[i] - X[start] > cutoff:
                start += 1
            while stop < dim and X[stop] - X[i] < cutoff:
                stop += 1
            for j in range(start, stop):
                weight = exp(-pow(X[i] - X[j], 2.0) / (2.0 * pow(width, 2.0)))
                total_weight += weight
                smoothed[i] += weight * Y[j]
            smoothed[i] /= total_weight
    return None


