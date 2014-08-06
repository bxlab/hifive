# distutils: language = c++
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

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
def find_possible_interactions(
        np.ndarray[DTYPE_int_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] chr_indices not None,
        int maxdistance):
    cdef int h, i, j, chrom_sum
    cdef int num_chromosomes = chr_indices.shape[0] - 1
    with nogil:
        # cycle through each chromosome
        for h in range(num_chromosomes):
            if maxdistance > 0:
                # cycle through each fend in chromosome, skipping if filtered
                for i in range(chr_indices[h], chr_indices[h + 1] - 2):
                    if filter[i] == 0:
                        continue
                    # check same strand adjacent fragment fends
                    j = i + 2
                    if filter[j] == 1 and (maxdistance == 0 or mids[j] - mids[i] <= maxdistance):
                        interactions[i] += 1
                        interactions[j] += 1
                    # check non-adjacent fragment fends
                    j = (i / 2 + 2) * 2
                    while j < chr_indices[h + 1] and (maxdistance == 0 or mids[j] - mids[i] <= maxdistance):
                        if filter[j] == 1:
                            interactions[i] += 1
                            interactions[j] += 1
                        j += 1
            else:
                chrom_sum = 0
                for i in range(chr_indices[h], chr_indices[h + 1]):
                    chrom_sum += filter[i]
                for i in range(chr_indices[h], chr_indices[h + 1]):
                    if filter[i] == 0:
                        continue
                    interactions[i] += chrom_sum
                    # remove same fragment fend if needed
                    if i % 2 == 0 and filter[i + 1] != 0:
                        interactions[i] -= 1
                        interactions[i + 1] -= 1
                    # remove adjacent fragment same strand fends if needed
                    if i / 2 < (chr_indices[h + 1] - 1) / 2:
                        j = i + 3 - (i % 2) * 2
                        if filter[j] != 0:
                            interactions[i] -= 1
                            interactions[j] -= 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fend_distance_bin_values(
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] data_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] distance_bins not None,
        np.ndarray[DTYPE_64_t, ndim=1] bin_sums not None,
        np.ndarray[DTYPE_int64_t, ndim=1] bin_counts not None,
        int fend,
        int stopfend,
        int corrected):
    cdef int i, j, distance, fend2
    cdef int num_bins = distance_bins.shape[0]
    with nogil:
        j = 0
        if corrected == 0:
            for i in range(data_indices[fend], data_indices[fend + 1]):
                fend2 = data[i, 1]
                if filter[fend2] == 0 or fend2 >= stopfend:
                    continue
                distance = mids[fend2] - mids[fend]
                while j < num_bins - 1 and distance > distance_bins[j]:
                    j += 1
                bin_sums[j] += data[i, 2]
        else:
            for i in range(data_indices[fend], data_indices[fend + 1]):
                fend2 = data[i, 1]
                if filter[fend2] == 0 or fend2 >= stopfend:
                    continue
                distance = mids[fend2] - mids[fend]
                while j < num_bins - 1 and distance > distance_bins[j]:
                    j += 1
                bin_sums[j] += data[i, 2] / (corrections[fend] * corrections[fend2])
        fend2 = fend + 2
        j = 0
        if fend2 < stopfend and filter[fend2] == 1:
            distance = mids[fend2] - mids[fend]
            while j < num_bins - 1 and distance > distance_bins[j]:
                j += 1
            bin_counts[j] += 1
        for fend2 in range((fend / 2 + 2) * 2, stopfend):
            if filter[fend2] == 0:
                continue
            distance = mids[fend2] - mids[fend]
            while j < num_bins - 1 and distance > distance_bins[j]:
                j += 1
            bin_counts[j] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_mindistance_interactions(
        np.ndarray[DTYPE_int64_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_int_t, ndim=1] chr_indices not None,
        np.ndarray[DTYPE_int_t, ndim=2] min_fend not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        int useread_int):
    cdef int i, j, k, chr_total, temp, all_total
    cdef int num_chroms = chr_indices.shape[0] - 1
    cdef int num_fends = filter.shape[0]
    with nogil:
        all_total = 0
        chr_total = 0
        if useread_int != 0:
            for i in range(num_fends):
                all_total += filter[i]
        for i in range(num_chroms):
            chr_total = 0
            for j in range(chr_indices[i], chr_indices[i + 1]):
                chr_total += filter[j]
            if useread_int > 0:
                for j in range(chr_indices[i], chr_indices[i + 1]):
                    if filter[j] > 0:
                        interactions[j] += all_total - chr_total
            if useread_int < 2:
                for j in range(chr_indices[i], chr_indices[i + 1]):
                    if filter[j] == 0:
                        continue
                    temp = chr_total
                    for k in range(min_fend[j, 0], min_fend[j, 1]):
                        temp -= filter[k]
                    if j % 2 == 0:
                        if filter[j + 1] == 1 and min_fend[j, 1] <= j + 1:
                            temp -= 1
                        if j - 1 >= chr_indices[i] and filter[j - 1] == 1 and min_fend[j, 0] > j - 1:
                            temp -= 1
                        if j + 3 < chr_indices[i + 1] and filter[j + 3] == 1 and min_fend[j, 1] <= j + 3:
                            temp -= 1
                    else:
                        if filter[j - 1] == 1 and min_fend[j, 0] > j - 1:
                            temp -= 1
                        if j - 3 >= chr_indices[i] and filter[j - 3] == 1 and min_fend[j, 0] > j - 3:
                            temp -= 1
                        if j + 1 < chr_indices[i + 1] and filter[j + 1] == 1 and min_fend[j, 1] <= j + 1:
                            temp -= 1
                    interactions[j] += temp
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_max_fend(
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] chromosomes not None,
        np.ndarray[DTYPE_int_t, ndim=1] chr_indices not None,
        int start,
        int maxdistance):
    cdef int h, i, j, chr_max
    cdef int num_fends = max_fend.shape[0]
    cdef int max_num_fends = mids.shape[0]
    with nogil:
        j = 1
        for i in range(num_fends):
            chr_max = min(max_num_fends, chr_indices[chromosomes[i] + 1] - start)
            if maxdistance == 0:
                max_fend[i] = chr_max
            else:
                while j < chr_max and mids[j] - mids[i] < maxdistance:
                    j += 1
                max_fend[i] = j
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_interaction_distance_means(
        np.ndarray[DTYPE_t, ndim=2] means not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_means not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None):
    cdef int i, j, k, distance
    cdef double frac
    cdef int num_bins = distance_mids.shape[0]
    cdef int num_fends = means.shape[0]
    with nogil:
        for i in range(num_fends):
            if filter[i] == 0:
                continue
            j = i + 2
            k = 0
            if j < max_fend[i] and filter[j] == 1:
                distance = mids[j] - mids[i]
                while k < num_bins and distance > distance_mids[k]:
                    k += 1
                if k == 0:
                    means[i, j - i - 1] = distance_means[0]
                elif k < num_bins:
                    frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                    means[i, j - i - 1] = distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                else:
                    means[i, j - i - 1] = distance_means[num_bins - 1]
            j = (i / 2 + 2) * 2
            while j < max_fend[i]:
                if filter[j] == 1:
                    distance = mids[j] - mids[i]
                    while k < num_bins and distance > distance_mids[k]:
                        k += 1
                    if k == 0:
                        means[i, j - i - 1] = distance_means[0]
                    elif k < num_bins:
                        frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                        means[i, j - i - 1] = distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                    else:
                        means[i, j - i - 1] = distance_means[num_bins - 1]
                j += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_data_distance_means(
        np.ndarray[DTYPE_t, ndim=1] means not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] data_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_means not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mids not None,
        int maxdistance):
    cdef int i, j, k, distance, fend2
    cdef double frac
    cdef int num_bins = distance_mids.shape[0]
    cdef int num_fends = filter.shape[0]
    with nogil:
        for i in range(num_fends - 1):
            if filter[i] == 0:
                continue
            k = 0
            for j in range(data_indices[i], data_indices[i + 1]):
                fend2 = data[j, 1]
                if filter[fend2] == 0:
                    continue
                distance = mids[fend2] - mids[i]
                if maxdistance == 0 or distance <= maxdistance:
                    while k < num_bins and distance > distance_mids[k]:
                        k += 1
                    if k == 0:
                        means[j] = distance_means[0]
                    elif k < num_bins:
                        frac = (distance - distance_mids[k - 1]) / (distance_mids[k] - distance_mids[k - 1])
                        means[j] = distance_means[k] * frac + distance_means[k - 1] * (1.0 - frac)
                    else:
                        means[j] = distance_means[num_bins - 1]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_gradients(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=2] interaction_means not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] gradients not None,
        np.ndarray[DTYPE_int_t, ndim=1] max_fend not None,
        int findcost):
    cdef int fend1, fend2, i, max_index, fend1_max
    cdef double value, cost, distance_mean
    cdef int num_fends = max_fend.shape[0]
    with nogil:
        cost = 0.0
        for fend1 in range(num_fends):
            if filter[fend1] == 0:
                continue
            i = indices[fend1]
            max_index = indices[fend1 + 1]
            fend1_max = max_fend[fend1]
            fend2 = fend1 + 2
            if fend2 < fend1_max and filter[fend2] == 1 and interaction_means[fend1, fend2 - fend1 - 1] > 0:
                distance_mean = interaction_means[fend1, fend2 - fend1 - 1]
                while i < max_index and data[i, 1] < fend2:
                    i += 1
                if i < max_index and data[i, 1] == fend2:
                    gradients[fend1] += distance_mean * corrections[fend2] - data[i, 2] / corrections[fend1]
                    gradients[fend2] += distance_mean * corrections[fend1] - data[i, 2] / corrections[fend2]
                    if findcost == 1:
                        value = distance_mean * corrections[fend1] * corrections[fend2]
                        cost += value - data[i, 2] * log(value)
                else:
                    gradients[fend1] += distance_mean * corrections[fend2]
                    gradients[fend2] += distance_mean * corrections[fend1]
                    if findcost == 1:
                        cost += distance_mean * corrections[fend1] * corrections[fend2]
            fend2 = (fend1 / 2 + 2) * 2
            while fend2 < fend1_max:
                if filter[fend2] == 1 and interaction_means[fend1, fend2 - fend1 - 1] > 0:
                    distance_mean = interaction_means[fend1, fend2 - fend1 - 1]
                    while i < max_index and data[i, 1] < fend2:
                        i += 1
                    if i < max_index and data[i, 1] == fend2:
                        gradients[fend1] += distance_mean * corrections[fend2] - data[i, 2] / corrections[fend1]
                        gradients[fend2] += distance_mean * corrections[fend1] - data[i, 2] / corrections[fend2]
                        if findcost == 1:
                            value = distance_mean * corrections[fend1] * corrections[fend2]
                            cost += value - data[i, 2] * log(value)
                    else:
                        gradients[fend1] += distance_mean * corrections[fend2]
                        gradients[fend2] += distance_mean * corrections[fend1]
                        if findcost == 1:
                            cost += distance_mean * corrections[fend1] * corrections[fend2]
                fend2 += 1
    return cost


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def update_corrections(
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_t, ndim=1] gradients not None,
        np.ndarray[DTYPE_int_t, ndim=1] interactions not None,
        double learningrate):
    cdef int i
    cdef int num_fends = corrections.shape[0]
    with nogil:
        for i in range(num_fends):
            if interactions[i] > 0:
                corrections[i] = max(0.01, corrections[i] - learningrate * gradients[i] / interactions[i])
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fend_means(
        np.ndarray[DTYPE_t, ndim=1] distance_means,
        np.ndarray[DTYPE_int64_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_64_t, ndim=1] fend_means not None,
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=2] trans_data,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        double mu,
        double trans_mu):
    cdef int i, fend1, fend2, num_trans_data, num_data
    cdef double cost, temp
    cdef int num_fends = filter.shape[0]
    if not trans_data is None:
        num_trans_data = trans_data.shape[0]
    else:
        num_trans_data = 0
    if not data is None:
        num_data = data.shape[0]
    else:
        num_data = 0
    with nogil:
        cost = 0.0
        for i in range(num_fends):
            fend_means[i] = 0.0
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
            fend_means[fend1] += temp
            fend_means[fend2] += temp
        for i in range(num_fends):
            if filter[i] == 0:
                continue
            temp = fend_means[i] / interactions[i]
            corrections[i] *= pow(temp, 0.5)
            cost += pow(temp - 1.0, 2.0)
        cost = pow(cost, 0.5)
    return cost
