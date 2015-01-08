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
        int useread_int,
        int binned):
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
                    if binned == 0:
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
def find_min_fend(
        np.ndarray[DTYPE_int_t, ndim=1] min_fend not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] chromosomes not None,
        np.ndarray[DTYPE_int_t, ndim=1] chr_indices not None,
        int mindistance):
    cdef int h, i, j, chr_max
    cdef int num_fends = min_fend.shape[0]
    cdef int max_num_fends = mids.shape[0]
    with nogil:
        j = 2
        for i in range(num_fends):
            chr_max = min(max_num_fends, chr_indices[chromosomes[i] + 1])
            if mindistance == 0:
                min_fend[i] = min(i + 2, chr_max)
            else:
                while j < chr_max and mids[j] - mids[i] < mindistance:
                    j += 1
                min_fend[i] = j
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_interaction_distance_means(
        np.ndarray[DTYPE_t, ndim=2] means not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_means not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mean_logs not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mids not None,
        np.ndarray[DTYPE_t, ndim=1] distance_mid_logs not None,
        np.ndarray[DTYPE_t, ndim=1] distance_log_spaces not None,
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
                    frac = (log(distance) - distance_mid_logs[k - 1]) / distance_log_spaces[k - 1]
                    means[i, j - i - 1] = exp(distance_mean_logs[k] * frac + distance_mean_logs[k - 1] * (1.0 - frac))
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
                        frac = (log(distance) - distance_mid_logs[k - 1]) / distance_log_spaces[k - 1]
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
        np.ndarray[DTYPE_int_t, ndim=1] chroms not None,
        np.ndarray[DTYPE_t, ndim=2] parameters not None,
        np.ndarray[DTYPE_t, ndim=1] chrom_means not None):
    cdef int i, j, k, fend2
    cdef double distance, chrom_mean
    cdef int num_fends = filter.shape[0]
    with nogil:
        for i in range(num_fends - 1):
            if filter[i] == 0:
                continue
            chrom_mean = chrom_means[chroms[i]]
            k = 0
            for j in range(data_indices[i], data_indices[i + 1]):
                fend2 = data[j, 1]
                if filter[fend2] == 0:
                    continue
                distance = log(mids[fend2] - mids[i])
                while distance > parameters[k, 0]:
                    k += 1
                means[j] = exp(parameters[k, 1] * distance + parameters[k, 2] + chrom_mean)
    return None


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


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fend_ranges(
        np.ndarray[DTYPE_int_t, ndim=1] mapping not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int64_t, ndim=2] ranges not None,
        int mindistance,
        int maxdistance,
        long long int start,
        long long int stop,
        int binned):
    cdef long long int i, j, pos, total
    cdef long long int total_fends = mapping.shape[0]
    with nogil:
        for i in range(start, stop):
            pos = mids[i]
            j = i + 1
            while j < total_fends and mids[j] - pos < mindistance:
                j += 1
            ranges[i, 1] = j
            j = total_fends
            while j > i and mids[j - 1] - pos > maxdistance:
                j -= 1
            ranges[i, 2] = j
            total = ranges[i, 2] - ranges[i, 1]
            if binned == 0:
                if ranges[i, 1] < total_fends and mapping[ranges[i, 1]] - mapping[i] == 1:
                    total -= 1
                if mapping[i] % 2 == 0:
                    if ranges[i, 1] < total_fends and mapping[ranges[i, 1]] - mapping[i] == 3:
                        total -= 1
                    if ranges[i, 1] + 1 < total_fends and mapping[ranges[i, 1] + 1] - mapping[i] == 3:
                        total -= 1
                    if ranges[i, 1] + 2 < total_fends and mapping[ranges[i, 1] + 2] - mapping[i] == 3:
                        total -= 1
            ranges[i, 0] = total
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_node_indices(
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        np.ndarray[DTYPE_int_t, ndim=1] indices0,
        np.ndarray[DTYPE_int_t, ndim=1] indices1,
        np.ndarray[DTYPE_int64_t, ndim=2] ranges,
        long long int start,
        long long int stop,
        int binned):
    cdef long long int total, i, j, pos
    cdef long long int num_fends = mapping.shape[0]
    with nogil:
        pos = 0
        total = 0
        i = 0
        while total + ranges[i, 0] < start:
            total += ranges[i, 0]
            i += 1
        while total < stop:
            j = ranges[i, 1]
            if binned == 0:
                while j < ranges[i, 2] and mapping[j] < mapping[i] + 4:
                    if mapping[j] - mapping[i] == 1:
                        j += 1
                    elif mapping[j] - mapping[i] == 3 and mapping[i] % 2 == 0:
                        j += 1
                    elif total < start:
                        j += 1
                        total += 1
                    elif total < stop:
                        indices0[pos] = i
                        indices1[pos] = j
                        j += 1
                        total += 1
                        pos += 1
                    else:
                        j += 1
            while j < ranges[i, 2]:
                if total < start:
                    j += 1
                    total += 1
                elif total < stop:
                    indices0[pos] = i
                    indices1[pos] = j
                    j += 1
                    total += 1
                    pos += 1
                else:
                    j += 1 
            i += 1
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
            distance = log(mids[indices1[i]] - mids[index0])
            if index0 != previous_index:
                previous_index = index0
                j = 0
            while distance > parameters[j, 0]:
                j += 1
            means[i] = exp(parameters[j, 1] * distance + parameters[j, 2] + chrom_mean)
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_gradients(
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices1 not None,
        np.ndarray[DTYPE_t, ndim=1] nonzero_means not None,
        np.ndarray[DTYPE_t, ndim=1] zero_means not None,
        np.ndarray[DTYPE_int_t, ndim=1] counts not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None,
        np.ndarray[DTYPE_64_t, ndim=1] gradients not None):
    cdef long long int i, index0, index1, count
    cdef double value, distance_mean, correction0, correction1
    cdef long long int num_zero_pairs = zero_indices0.shape[0]
    cdef long long int num_nonzero_pairs = nonzero_indices0.shape[0]
    with nogil:
        cost = 0.0
        for i in range(num_nonzero_pairs):
            index0 = nonzero_indices0[i]
            index1 = nonzero_indices1[i]
            count = counts[i]
            distance_mean = nonzero_means[i]
            correction0 = corrections[index0]
            correction1 = corrections[index1]
            gradients[index0] += distance_mean * correction1 - count / correction0
            gradients[index1] += distance_mean * correction0 - count / correction1
        for i in range(num_zero_pairs):
            index0 = zero_indices0[i]
            index1 = zero_indices1[i]
            distance_mean = zero_means[i]
            correction0 = corrections[index0]
            correction1 = corrections[index1]
            gradients[index0] += distance_mean * correction1
            gradients[index1] += distance_mean * correction0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calculate_cost(
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] zero_indices1 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] nonzero_indices1 not None,
        np.ndarray[DTYPE_t, ndim=1] nonzero_means not None,
        np.ndarray[DTYPE_t, ndim=1] zero_means not None,
        np.ndarray[DTYPE_int_t, ndim=1] counts not None,
        np.ndarray[DTYPE_t, ndim=1] corrections not None):
    cdef long long int i
    cdef double value, cost
    cdef long long int num_zero_pairs = zero_indices0.shape[0]
    cdef long long int num_nonzero_pairs = nonzero_indices0.shape[0]
    with nogil:
        cost = 0.0
        for i in range(num_nonzero_pairs):
            value = nonzero_means[i] * corrections[nonzero_indices0[i]] * corrections[nonzero_indices1[i]]
            cost += value - counts[i] * log(value)
        for i in range(num_zero_pairs):
            cost += zero_means[i] * corrections[zero_indices0[i]] * corrections[zero_indices1[i]]
    return cost


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
        np.ndarray[DTYPE_64_t, ndim=2] count_sum not None,
        np.ndarray[DTYPE_64_t, ndim=2] logdistance_sum not None,
        int start,
        int stop,
        int index,
        int binned):
    cdef int i, j, temp, fend1, fend2, previous_fend
    cdef double log_dist
    cdef int num_data = indices.shape[0]
    cdef int num_fends = rev_mapping.shape[0]
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
            log_dist = log(mids[fend2] - mids[fend1])
            while log_dist > cutoffs[j]:
                j += 1
            count_sum[index, j] += counts[i]
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
                    log_dist = log(mids[fend2] - mids[fend1])
                    while log_dist > cutoffs[j]:
                        j += 1
                    bin_size[index, j] += 1
                    logdistance_sum[index, j] += log_dist
            else:
                for fend2 in range(fend1 + 1, min(fend1 + 4, num_fends)):
                    log_dist = log(mids[fend2] - mids[fend1])
                    while log_dist > cutoffs[j]:
                        j += 1
                    bin_size[index, j] += 1
                    logdistance_sum[index, j] += log_dist
            for fend2 in range(fend1 + 4, num_fends):
                log_dist = log(mids[fend2] - mids[fend1])
                while log_dist > cutoffs[j]:
                    j += 1
                bin_size[index, j] += 1
                logdistance_sum[index, j] += log_dist
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_chromosome_sums(
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_t, ndim=2] parameters not None,
        np.ndarray[DTYPE_64_t, ndim=1] sums not None,
        int h):
    cdef int i, j, fend1, fend2, previous_index
    cdef double log_dist
    cdef int num_data = data.shape[0]
    with nogil:
        previous_index = -1
        j = 0
        for i in range(num_data):
            fend1 = data[i, 0]
            if fend1 != previous_index:
                j = 0
                previous_index = fend1
            fend2 = data[i, 1]
            log_dist = log(mids[fend2] - mids[fend1])
            while log_dist > parameters[j, 0]:
                j += 1
            sums[h] += data[i, 2] / exp(parameters[j, 1] * log_dist + parameters[j, 2])
    return None
