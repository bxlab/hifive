# distutils: language = c++

"""These functions provide increased speed in handling the distance-
dependent signal functions necessary for supporting fivec analysis using
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

cdef extern from "_normal.hpp":
    double stdnorm_cdf( double z ) nogil
    double stdnorm_pdf( double z ) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_max_frag(
        np.ndarray[DTYPE_int_t, ndim=1] max_fragment not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] starts not None,
        np.ndarray[DTYPE_int_t, ndim=1] stops not None,
        int maxdistance):
    cdef int i, j, k
    cdef int num_regions = starts.shape[0]
    with nogil:
        for i in range(num_regions):
            k = starts[i] + 1
            for j in range(starts[i], stops[i]):
                if maxdistance == 0:
                    max_fragment[j] = stops[i]
                else:
                    while k < stops[i] and mids[k] - mids[j] < maxdistance:
                        k += 1
                    max_fragment[j] = k
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fragment_interactions(
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fragment not None):
    cdef int i, j, k
    cdef int num_frags = filter.shape[0]
    with nogil:
        for i in range(num_frags - 1):
            if filter[i] == 0:
                continue
            for j in range(indices[i], indices[i + 1]):
                k = data[j, 1]
                if filter[k] > 0 and k < max_fragment[i]:
                    interactions[i] += 1
                    interactions[k] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_gradients(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts_n not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts_p not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=1] distance_signal not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] gradients not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fragment not None,
        double sigma,
        int find_cost):
    cdef int frag1, frag2, i
    cdef double value0, value, cost, expected, z_p, z_n, pdf_p, pdf_n, cdf_p, cdf_n
    cdef int num_frags = filter.shape[0]
    cdef double sigma_2 = pow(sigma, 2.0)
    with nogil:
        cost = 0.0
        for frag1 in range(num_frags-1):
            if filter[frag1] == 0:
                continue
            for i in range(indices[frag1], indices[frag1 + 1]):
                frag2 = data[i, 1]
                if filter[frag2] == 0 or frag2 >= max_fragment[frag1]:
                    continue
                expected = max(0.01, corrections[frag1] + corrections[frag2] + distance_signal[i])
                z_p = (log_counts_p[i] - expected) / sigma
                cdf_p = stdnorm_cdf(z_p)
                pdf_p = stdnorm_pdf(z_p)
                z_n = (log_counts_n[i] - expected) / sigma
                cdf_n = stdnorm_cdf(z_n)
                pdf_n = stdnorm_pdf(z_n)
                value0 = (cdf_p - cdf_n) * sigma
                if value0 != 0:
                    value = (pdf_p - pdf_n) / value0
                    gradients[frag1] += value
                    gradients[frag2] += value
                    if find_cost:
                        cost -= log(cdf_p - cdf_n)
                else:
                    value = (expected - log_counts[i]) / sigma_2
                    gradients[frag1] += value
                    gradients[frag2] += value
                    if find_cost:
                        value0 = pow(log_counts[i] - expected, 2.0) / sigma_2
                        cost += value0 * 0.5
    return cost


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def update_corrections(
    np.ndarray[DTYPE_int_t, ndim=1] filter not None,
    np.ndarray[DTYPE_t, ndim=1] corrections not None,
    np.ndarray[DTYPE_t, ndim=1] gradients not None,
    np.ndarray[DTYPE_int_t, ndim=1] interactions not None,
    double learning_rate ):
    cdef int i
    cdef int num_frags = filter.shape[0]
    with nogil:
        for i in range(num_frags):
            if filter[i] == 0:
                continue
            corrections[i] = max(-10.0, min(10.0, corrections[i] - learning_rate * gradients[i] / interactions[i]))
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fragment_means(
        np.ndarray[DTYPE_t, ndim=1] distance_means,
        np.ndarray[DTYPE_int_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_64_t, ndim=1] fragment_means not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None):
    cdef int i, frag1, frag2
    cdef double cost, temp
    cdef int num_frags = interactions.shape[0]
    cdef int num_data = data.shape[0]
    with nogil:
        cost = 0.0
        for i in range(num_frags):
            fragment_means[i] = 0.0
        for i in range(num_data):
            frag1 = data[i, 0]
            frag2 = data[i, 1]
            if distance_means is None:
                temp = log_counts[i] - corrections[frag1] - corrections[frag2]
            else:
                temp = log_counts[i] - distance_means[i] - corrections[frag1] - corrections[frag2]
            fragment_means[frag1] += temp
            fragment_means[frag2] += temp
        for i in range(num_frags):
            if interactions[i] == 0:
                continue
            temp = fragment_means[i] / interactions[i]
            corrections[i] += temp * 0.5
            cost += pow(temp, 2.0)
        cost = pow(cost, 0.5)
    return cost
