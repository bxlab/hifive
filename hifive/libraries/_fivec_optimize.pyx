# distutils: language = c++

"""These functions provide increased speed in handling the correction
optimization functions necessary for supporting fivec analysis using HiFive.
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
    double stdnorm_pdf(double x) nogil
    double stdnorm_cdf(double x) nogil


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_binning_correction_adjustment(
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=2] indices not None,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] correction_indices,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=2] frag_indices):
    cdef long long int i, frag1, frag2, bin1, bin2, index
    cdef long long int num_frags = corrections.shape[0]
    cdef long long int num_parameters = frag_indices.shape[1]
    with nogil:
        for i in range(num_frags):
            frag1 = indices[i, 0]
            frag2 = indices[i, 1]
            for j in range(num_parameters):
                bin1 = min(frag_indices[frag1, j], frag_indices[frag2, j])
                bin2 = max(frag_indices[frag1, j], frag_indices[frag2, j])
                index = (num_bins[j] - 1) * bin1 - bin1 * (bin1 - 1) / 2 + bin2 + correction_indices[j]
                corrections[i] *= binning_corrections[index]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_gradients(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts_n not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts_p not None,
        np.ndarray[DTYPE_t, ndim=1] distance_signal not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] gradients not None,
        double sigma):
    cdef long long int frag1, frag2, h, i
    cdef double value0, value, expected, z_p, z_n, pdf_p, pdf_n, cdf_p, cdf_n
    cdef long long int num_data = data.shape[0]
    cdef double sigma_2 = pow(sigma, 2.0)
    with nogil:
        for i in range(num_data):
            frag1 = data[i, 0]
            frag2 = data[i, 1]
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
            else:
                value = (expected - log_counts[i]) / sigma_2
                gradients[frag1] += value
                gradients[frag2] += value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_prob_cost(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts_n not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts not None,
        np.ndarray[DTYPE_t, ndim=1] log_counts_p not None,
        np.ndarray[DTYPE_t, ndim=1] distance_signal not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        double sigma):
    cdef long long int frag1, frag2, i
    cdef double value, cost, expected, z_p, z_n, cdf_p, cdf_n
    cdef long long int num_data = data.shape[0]
    cdef double sigma_2 = pow(sigma, 2.0)
    with nogil:
        cost = 0.0
        for i in range(num_data):
            frag1 = data[i, 0]
            frag2 = data[i, 1]
            expected = max(0.01, corrections[frag1] + corrections[frag2] + distance_signal[i])
            z_p = (log_counts_p[i] - expected) / sigma
            cdf_p = stdnorm_cdf(z_p)
            z_n = (log_counts_n[i] - expected) / sigma
            cdf_n = stdnorm_cdf(z_n)
            value = (cdf_p - cdf_n) * sigma
            if value != 0:
                cost -= log(cdf_p - cdf_n)
            else:
                cost += 0.5 * pow(log_counts[i] - expected, 2.0) / sigma_2
    return cost


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def update_corrections(
    np.ndarray[DTYPE_int_t, ndim=1] filter not None,
    np.ndarray[DTYPE_t, ndim=1] corrections not None,
    np.ndarray[DTYPE_t, ndim=1] new_corrections not None,
    np.ndarray[DTYPE_t, ndim=1] gradients not None,
    double learning_rate ):
    cdef long long int i
    cdef long long int num_frags = filter.shape[0]
    with nogil:
        for i in range(num_frags):
            if filter[i] == 0:
                continue
            new_corrections[i] = max(-10.0, min(10.0, corrections[i] - learning_rate * gradients[i]))
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_log_fragment_means(
        np.ndarray[DTYPE_t, ndim=1] distance_means,
        np.ndarray[DTYPE_t, ndim=1] trans_means,
        np.ndarray[DTYPE_int_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_64_t, ndim=1] fragment_means not None,
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=2] trans_data,
        np.ndarray[DTYPE_64_t, ndim=1] log_counts,
        np.ndarray[DTYPE_64_t, ndim=1] trans_log_counts,
        np.ndarray[DTYPE_t, ndim=1] corrections not None):
    cdef long long int i, frag1, frag2, num_data, num_trans
    cdef double cost, temp
    cdef long long int num_frags = interactions.shape[0]
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    if not trans_data is None:
        num_trans = trans_data.shape[0]
    else:
        num_trans = 0
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
        for i in range(num_trans):
            frag1 = trans_data[i, 0]
            frag2 = trans_data[i, 1]
            if trans_means is None:
                temp = trans_log_counts[i] - corrections[frag1] - corrections[frag2]
            else:
                temp = trans_log_counts[i] - trans_means[i] - corrections[frag1] - corrections[frag2]
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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fragment_means(
        np.ndarray[DTYPE_t, ndim=1] distance_means,
        np.ndarray[DTYPE_t, ndim=1] trans_means,
        np.ndarray[DTYPE_int_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_64_t, ndim=1] fragment_means not None,
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=2] trans_data,
        np.ndarray[DTYPE_64_t, ndim=1] log_counts,
        np.ndarray[DTYPE_64_t, ndim=1] log_trans_counts,
        np.ndarray[DTYPE_t, ndim=1] corrections not None):
    cdef long long int i, frag1, frag2, num_data, num_trans
    cdef double cost, temp
    cdef long long int num_frags = interactions.shape[0]
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    if not trans_data is None:
        num_trans = trans_data.shape[0]
    else:
        num_trans = 0
    with nogil:
        cost = 0.0
        for i in range(num_frags):
            fragment_means[i] = 0.0
        for i in range(num_data):
            frag1 = data[i, 0]
            frag2 = data[i, 1]
            if distance_means is None:
                temp = exp(log_counts[i] - corrections[frag1] - corrections[frag2])
            else:
                temp = exp(log_counts[i] - distance_means[i] - corrections[frag1] - corrections[frag2])
            fragment_means[frag1] += temp
            fragment_means[frag2] += temp
        for i in range(num_trans):
            frag1 = trans_data[i, 0]
            frag2 = trans_data[i, 1]
            if trans_means is None:
                temp = exp(log_trans_counts[i] - corrections[frag1] - corrections[frag2])
            else:
                temp = exp(log_trans_counts[i] - trans_means[i] - corrections[frag1] - corrections[frag2])
            fragment_means[frag1] += temp
            fragment_means[frag2] += temp
        for i in range(num_frags):
            if interactions[i] == 0:
                continue
            temp = fragment_means[i] / interactions[i]
            corrections[i] += log(temp) * 0.5
            cost += pow(temp - 1.0, 2.0)
        cost = pow(cost, 0.5)
    return cost


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_v(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=2] trans_data,
        np.ndarray[DTYPE_64_t, ndim=1] counts,
        np.ndarray[DTYPE_64_t, ndim=1] trans_counts,
        np.ndarray[DTYPE_64_t, ndim=2] corrections,
        np.ndarray[DTYPE_64_t, ndim=2] v):
    cdef long long int i, frag1, frag2, num_data, num_trans_data
    cdef double correction
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    if not trans_data is None:
        num_trans = trans_data.shape[0]
    else:
        num_trans = 0
    with nogil:
        for i in range(num_data):
            frag1 = data[i, 0]
            frag2 = data[i, 1]
            correction = corrections[frag1, 0] * corrections[frag2, 0] * counts[i]
            v[frag1, 0] += correction
            v[frag2, 0] += correction
        for i in range(num_trans):
            frag1 = trans_data[i, 0]
            frag2 = trans_data[i, 1]
            correction = corrections[frag1, 0] * corrections[frag2, 0] * trans_counts[i]
            v[frag1, 0] += correction
            v[frag2, 0] += correction
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_w(
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=2] trans_data,
        np.ndarray[DTYPE_64_t, ndim=1] counts,
        np.ndarray[DTYPE_64_t, ndim=1] trans_counts,
        np.ndarray[DTYPE_64_t, ndim=2] corrections,
        np.ndarray[DTYPE_64_t, ndim=2] p,
        np.ndarray[DTYPE_64_t, ndim=2] w):
    cdef long long int i, frag1, frag2, num_data, num_trans_data
    cdef double correction
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    if not trans_data is None:
        num_trans = trans_data.shape[0]
    else:
        num_trans = 0
    with nogil:
        for i in range(num_data):
            frag1 = data[i, 0]
            frag2 = data[i, 1]
            correction = corrections[frag1, 0] * corrections[frag2, 0] * counts[i]
            w[frag1, 0] += correction * p[frag2, 0]
            w[frag2, 0] += correction * p[frag1, 0]
        for i in range(num_trans):
            frag1 = trans_data[i, 0]
            frag2 = trans_data[i, 1]
            correction = corrections[frag1, 0] * corrections[frag2, 0] * trans_counts[i]
            w[frag1, 0] += correction * p[frag2, 0]
            w[frag2, 0] += correction * p[frag1, 0]
    return None
