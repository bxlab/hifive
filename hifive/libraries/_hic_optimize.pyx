# distutils: language = c++

"""These functions provide increased speed in handling functions dealing with correction
parameter optimization necessary for supporting hic analysis using HiFive.
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
def find_binning_correction_adjustment(
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices1 not None,
        np.ndarray[DTYPE_t, ndim=1] binning_corrections,
        np.ndarray[DTYPE_int_t, ndim=1] num_bins,
        np.ndarray[DTYPE_int_t, ndim=3] fend_indices):
    cdef long long int i, j, fend1, fend2, index
    cdef long long int num_fends = corrections.shape[0]
    cdef long long int num_parameters = fend_indices.shape[1]
    with nogil:
        for i in range(num_fends):
            fend1 = indices0[i]
            fend2 = indices1[i]
            for j in range(num_parameters):
                if fend_indices[fend1, j, 0] < fend_indices[fend2, j, 0]:
                    corrections[i] *= binning_corrections[fend_indices[fend1, j, 1] + fend_indices[fend2, j, 0]]
                else:
                    corrections[i] *= binning_corrections[fend_indices[fend2, j, 1] + fend_indices[fend1, j, 0]]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_binom_gradients(
        np.ndarray[DTYPE_int_t, ndim=1] counts,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices1 not None,
        np.ndarray[DTYPE_t, ndim=1] nonzero_means,
        np.ndarray[DTYPE_t, ndim=1] zero_means not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] inv_corrections not None,
        np.ndarray[DTYPE_64_t, ndim=1] gradients not None):
    cdef long long int i, index0, index1
    cdef double value, distance_mean, correction0, correction1
    cdef long long int num_zero_pairs = zero_indices0.shape[0]
    cdef long long int num_nonzero_pairs = nonzero_indices0.shape[0]
    with nogil:
        for i in range(num_nonzero_pairs):
            index0 = nonzero_indices0[i]
            index1 = nonzero_indices1[i]
            gradients[index0] -= inv_corrections[index0]
            if index1 != index0:
                gradients[index1] -= inv_corrections[index1]
        for i in range(num_zero_pairs):
            index0 = zero_indices0[i]
            index1 = zero_indices1[i]
            distance_mean = zero_means[i]
            correction0 = corrections[index0]
            correction1 = corrections[index1]
            value = 1.0 / (1.0 - distance_mean * corrections[index0] * corrections[index1])
            gradients[index0] += (distance_mean * corrections[index1]) * value
            if index1 != index0:
                gradients[index1] += (distance_mean * corrections[index0]) * value
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_binom_cost(
        np.ndarray[DTYPE_int_t, ndim=1] counts,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices1 not None,
        np.ndarray[DTYPE_t, ndim=1] nonzero_means not None,
        np.ndarray[DTYPE_t, ndim=1] zero_means not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] log_corrections not None):
    cdef long long int i
    cdef double cost
    cdef long long int num_zero_pairs = zero_indices0.shape[0]
    cdef long long int num_nonzero_pairs = nonzero_indices0.shape[0]
    with nogil:
        cost = 0.0
        for i in range(num_nonzero_pairs):
            cost -= nonzero_means[i] + log_corrections[nonzero_indices0[i]] + log_corrections[nonzero_indices1[i]]
        for i in range(num_zero_pairs):
            cost -= log(max(0.0000001, 1.0 - zero_means[i] * corrections[zero_indices0[i]] * corrections[zero_indices1[i]]))
    return cost


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_poisson_gradients(
        np.ndarray[DTYPE_int_t, ndim=1] counts not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices1 not None,
        np.ndarray[DTYPE_t, ndim=1] nonzero_means not None,
        np.ndarray[DTYPE_t, ndim=1] zero_means not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] inv_corrections not None,
        np.ndarray[DTYPE_64_t, ndim=1] gradients not None):
    cdef long long int i, index0, index1
    cdef long long int num_zero_pairs = zero_indices0.shape[0]
    cdef long long int num_nonzero_pairs = nonzero_indices0.shape[0]
    with nogil:
        for i in range(num_nonzero_pairs):
            index0 = nonzero_indices0[i]
            index1 = nonzero_indices1[i]
            gradients[index0] += nonzero_means[i] * corrections[index0] - counts[i] * inv_corrections[index0]
            if index1 != index0:
                gradients[index1] += nonzero_means[i] * corrections[index1] - counts[i] * inv_corrections[index1]
        for i in range(num_zero_pairs):
            index0 = zero_indices0[i]
            index1 = zero_indices1[i]
            gradients[index0] += zero_means[i] * corrections[index0]
            if index1 != index0:
                gradients[index1] += zero_means[i] * corrections[index1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_poisson_cost(
        np.ndarray[DTYPE_int_t, ndim=1] counts not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices1 not None,
        np.ndarray[DTYPE_t, ndim=1] nonzero_means not None,
        np.ndarray[DTYPE_t, ndim=1] zero_means not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] log_corrections not None):
    cdef long long int i
    cdef double cost
    cdef long long int num_zero_pairs = zero_indices0.shape[0]
    cdef long long int num_nonzero_pairs = nonzero_indices0.shape[0]
    with nogil:
        cost = 0.0
        for i in range(num_nonzero_pairs):
            cost += corrections[nonzero_indices0[i]] * corrections[nonzero_indices1[i]] * nonzero_means[i] - counts[i] * ( log(nonzero_means[i]) + log_corrections[nonzero_indices0[i]] + log_corrections[nonzero_indices1[i]] )
        for i in range(num_zero_pairs):
            cost += corrections[zero_indices0[i]] * corrections[zero_indices1[i]] * zero_means[i]
    return cost


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fend_means(
        np.ndarray[DTYPE_t, ndim=1] distance_means,
        np.ndarray[DTYPE_t, ndim=1] trans_means,
        np.ndarray[DTYPE_64_t, ndim=1] fend_means not None,
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=2] trans_data,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        double mu,
        double trans_mu,
        int binary):
    cdef long long int i, fend1, fend2, num_trans_data, num_data
    cdef double temp
    cdef long long int num_fends = fend_means.shape[0]
    if not trans_data is None:
        num_trans_data = trans_data.shape[0]
    else:
        num_trans_data = 0
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    with nogil:
        for i in range(num_fends):
            fend_means[i] = 0.0
        if binary == 0:
            if distance_means is None:
                for i in range(num_data):
                    fend1 = data[i, 0]
                    fend2 = data[i, 1]
                    temp = data[i, 2] / (mu * corrections[fend1] * corrections[fend2])
                    fend_means[fend1] += temp
                    fend_means[fend2] += temp
            else:
                for i in range(num_data):
                    fend1 = data[i, 0]
                    fend2 = data[i, 1]
                    temp = data[i, 2] / (distance_means[i] * corrections[fend1] * corrections[fend2])
                    fend_means[fend1] += temp
                    fend_means[fend2] += temp
            for i in range(num_trans_data):
                fend1 = trans_data[i, 0]
                fend2 = trans_data[i, 1]
                temp = trans_data[i, 2] / (trans_mu * corrections[fend1] * corrections[fend2])
                if not trans_means is None:
                    temp /= trans_means[i]
                fend_means[fend1] += temp
                fend_means[fend2] += temp
        else:
            if distance_means is None:
                for i in range(num_data):
                    fend1 = data[i, 0]
                    fend2 = data[i, 1]
                    temp = 1.0 / (mu * corrections[fend1] * corrections[fend2])
                    fend_means[fend1] += temp
                    fend_means[fend2] += temp
            else:
                for i in range(num_data):
                    fend1 = data[i, 0]
                    fend2 = data[i, 1]
                    temp = 1.0 / (distance_means[i] * corrections[fend1] * corrections[fend2])
                    fend_means[fend1] += temp
                    fend_means[fend2] += temp
            for i in range(num_trans_data):
                fend1 = trans_data[i, 0]
                fend2 = trans_data[i, 1]
                temp = 1.0 / (trans_mu * corrections[fend1] * corrections[fend2])
                if not trans_means is None:
                    temp /= trans_means[i]
                fend_means[fend1] += temp
                fend_means[fend2] += temp
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def update_express_corrections(
        np.ndarray[DTYPE_int_t, ndim=1] filt,
        np.ndarray[DTYPE_int64_t, ndim=1] interactions,
        np.ndarray[DTYPE_64_t, ndim=1] fend_means,
        np.ndarray[DTYPE_t, ndim=1] corrections,
        np.ndarray[DTYPE_64_t, ndim=1] change):
    cdef long long int i
    cdef double cost, temp
    cdef long long int num_fends = filt.shape[0]
    with nogil:
        cost = 0.0
        change[0] = 0
        for i in range(num_fends):
            if filt[i] == 0:
                continue
            temp = fend_means[i] / interactions[i]
            cost += pow(temp - 1.0, 2.0)
            temp = pow(temp, 0.5)
            if temp > 1.0:
                change[0] = max(change[0], temp - 1.0)
            else:
                change[0] = max(change[0], 1.0 - temp)
            corrections[i] *= temp
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
    cdef long long int i, fend1, fend2, num_data, num_trans_data
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
            fend1 = data[i, 0]
            fend2 = data[i, 1]
            correction = corrections[fend1, 0] * corrections[fend2, 0] * counts[i]
            v[fend1, 0] += correction
            v[fend2, 0] += correction
        for i in range(num_trans):
            fend1 = trans_data[i, 0]
            fend2 = trans_data[i, 1]
            correction = corrections[fend1, 0] * corrections[fend2, 0] * trans_counts[i]
            v[fend1, 0] += correction
            v[fend2, 0] += correction
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
    cdef long long int i, fend1, fend2, num_data, num_trans_data
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
            fend1 = data[i, 0]
            fend2 = data[i, 1]
            correction = corrections[fend1, 0] * corrections[fend2, 0] * counts[i]
            w[fend1, 0] += correction * p[fend2, 0]
            w[fend2, 0] += correction * p[fend1, 0]
        for i in range(num_trans):
            fend1 = trans_data[i, 0]
            fend2 = trans_data[i, 1]
            correction = corrections[fend1, 0] * corrections[fend2, 0] * trans_counts[i]
            w[fend1, 0] += correction * p[fend2, 0]
            w[fend2, 0] += correction * p[fend1, 0]
    return None
