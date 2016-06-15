# distutils: language = c++

"""These functions provide increased speed in handling functions dealing with finding valid
interactions and interaction ranges necessary for supporting hic analysis using HiFive.
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
def find_max_bin(
        np.ndarray[DTYPE_int_t, ndim=2] binbounds not None,
        int maxdistance):
    cdef int i, j, maxbin, start
    cdef int num_bins = binbounds.shape[0]
    with nogil:
        maxbin = 0
        for i in range(num_bins):
            start = binbounds[i, 1] - 1
            j = i + 1
            while j < num_bins and binbounds[j, 0] - start < maxdistance:
                j += 1
            maxbin = max(maxbin, j - i - 1)
    return maxbin

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
    cdef long long int h, i, j, chr_max
    cdef long long int num_fends = max_fend.shape[0]
    cdef long long int max_num_fends = mids.shape[0]
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
        int mindistance,
        int binned):
    cdef long long int h, i, j, chr_max
    cdef long long int num_fends = min_fend.shape[0]
    cdef long long int max_num_fends = mids.shape[0]
    with nogil:
        if binned == 1:
            j = 0
        else:
            j = 2
        for i in range(num_fends):
            chr_max = min(max_num_fends, chr_indices[chromosomes[i] + 1])
            if mindistance == 0:
                if binned == 1:
                    min_fend[i] = min(i, chr_max)
                else:
                    min_fend[i] = min(i + 2, chr_max)
            else:
                while j < chr_max and mids[j] - mids[i] < mindistance:
                    j += 1
                min_fend[i] = j
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_fend_coverage(
        np.ndarray[DTYPE_int_t, ndim=2] data not None,
        np.ndarray[DTYPE_int64_t, ndim=1] data_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] min_fend not None,
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
            while i < data_indices[fend1 + 1] and data[i, 1] < min_fend[fend1]:
                i += 1
            while i < data_indices[fend1 + 1] and data[i, 1] < max_fend[fend1]:
                if filter[data[i, 1]] == 1:
                    coverage[fend1] += 1
                    if fend1 != data[i, 1]:
                        coverage[data[i, 1]] += 1
                i += 1
    return None


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
        int startfend,
        int binned):
    cdef long long int i, j, pos, total
    cdef long long int total_fends = mapping.shape[0]
    with nogil:
        for i in range(start, stop):
            pos = mids[i]
            j = i + 1 - binned
            while j < total_fends and mids[j] - pos < mindistance:
                j += 1
            ranges[i, 1] = j
            j = total_fends
            while j > ranges[i, 1] and mids[j - 1] - pos > maxdistance:
                j -= 1
            ranges[i, 2] = j
            total = ranges[i, 2] - ranges[i, 1]
            if binned == 0:
                if ranges[i, 1] < total_fends and mapping[ranges[i, 1]] - mapping[i] == 1:
                    total -= 1
                if (mapping[i] + startfend) % 2 == 0:
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
def find_nonzeros_in_range(
        np.ndarray[DTYPE_int64_t, ndim=2] ranges not None,
        np.ndarray[DTYPE_int_t, ndim=2] data not None):
    cdef long long int i, count
    cdef long long int num_data = data.shape[0]
    with nogil:
        count = 0
        for i in range(num_data):
            if data[i, 0] == -1:
                continue
            if data[i, 1] >= ranges[data[i, 0], 1] and data[i, 1] < ranges[data[i, 0], 2]:
                count += 1
    return count


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_nonzero_node_indices(
        np.ndarray[DTYPE_int64_t, ndim=2] ranges,
        np.ndarray[DTYPE_int_t, ndim=1] indices0,
        np.ndarray[DTYPE_int_t, ndim=1] indices1,
        np.ndarray[DTYPE_int_t, ndim=1] counts,
        np.ndarray[DTYPE_int_t, ndim=2] data,
        np.ndarray[DTYPE_int_t, ndim=1] interactions):
    cdef long long int i, j
    cdef long long int num_data = data.shape[0]
    with nogil:
        i = 0
        for j in range(num_data):
            if data[j, 0] == -1:
                continue
            if data[j, 1] >= ranges[data[j, 0], 1] and data[j, 1] < ranges[data[j, 0], 2]:
                indices0[i] = data[j, 0]
                indices1[i] = data[j, 1]
                interactions[data[j, 0]] += 1
                if data[j, 0] != data[j, 1]:
                    interactions[data[j, 1]] += 1
                if not counts is None:
                    counts[i] = data[j, 2]
                i += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_zero_node_indices(
        np.ndarray[DTYPE_int_t, ndim=1] mapping,
        np.ndarray[DTYPE_int64_t, ndim=2] ranges,
        np.ndarray[DTYPE_int_t, ndim=1] nzindices0,
        np.ndarray[DTYPE_int_t, ndim=1] nzindices1,
        np.ndarray[DTYPE_int_t, ndim=1] zindices0,
        np.ndarray[DTYPE_int_t, ndim=1] zindices1,
        np.ndarray[DTYPE_int_t, ndim=1] interactions,
        long long int start,
        long long int stop,
        int startfend,
        int binned):
    cdef long long int i, j, m1, m2, nzpos, zpos
    cdef long long int num_nz = nzindices0.shape[0]
    with nogil:
        nzpos = 0
        zpos = 0
        for i in range(start, stop):
            while nzpos < num_nz and nzindices0[nzpos] < i:
                nzpos += 1
            for j in range(ranges[i, 1], ranges[i, 2]):
                while nzpos < num_nz and nzindices0[nzpos] == i and nzindices1[nzpos] < j:
                    nzpos += 1
                m1 = mapping[i]
                m2 = mapping[j]
                if binned == 0 and m2 - m1 < 4 and (m2 - m1 == 1 or ((m1 + startfend) % 2 == 0 and m2 - m1 == 3)):
                    continue
                elif nzpos == num_nz or nzindices1[nzpos] != j or nzindices0[nzpos] != i:
                    zindices0[zpos] = i
                    zindices1[zpos] = j
                    interactions[i] += 1
                    if i != j:
                        interactions[j] += 1
                    zpos += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_distancebound_possible_interactions(
        np.ndarray[DTYPE_int64_t, ndim=1] interactions not None,
        np.ndarray[DTYPE_int_t, ndim=1] chr_indices not None,
        np.ndarray[DTYPE_int_t, ndim=1] mids not None,
        np.ndarray[DTYPE_int_t, ndim=1] filter not None,
        np.ndarray[DTYPE_int_t, ndim=1] chrints not None,
        int useread_int,
        int mindistance,
        int maxdistance,
        int binned):
    cdef long long int h, i, j, k, chr_total, all_total, total, stop
    cdef long long int num_chroms = chrints.shape[0]
    cdef long long int num_fends = filter.shape[0]
    with nogil:
        all_total = 0
        chr_total = 0
        # if trans interactions are included, find total valid fends
        if useread_int != 0:
            for i in range(num_fends):
                all_total += filter[i]
        for h in range(num_chroms):
            i = chrints[h]
            # for each chrom, find valid number of fends
            chr_total = 0
            for j in range(chr_indices[i], chr_indices[i + 1]):
                chr_total += filter[j]
            # if using trans, for each valid fend add on all possible trans interactions tot total
            if useread_int > 0:
                for j in range(chr_indices[i], chr_indices[i + 1]):
                    if filter[j] > 0:
                        interactions[j] += all_total - chr_total
            # if using cis, find number of interactions within distance range and
            # not adjacent fends or opposite strand adjacent fragment fends
            if useread_int < 2:
                for j in range(chr_indices[i], chr_indices[i + 1]):
                    if filter[j] == 0:
                        continue
                    if binned == 0:
                        total = chr_total - 1
                        if j > chr_indices[i] + 4:
                            # remove upstream interactions outside of maxdistance
                            k = chr_indices[i]
                            stop = (j / 2) * 2 - 2
                            while k < stop and mids[j] - mids[k] >= maxdistance:
                                total -= filter[k]
                                k += 1
                            stop = k - 1
                            # remove upstream interactions inside of mindistance
                            k = (j / 2) * 2 - 3
                            while k > stop and mids[j] - mids[k] < mindistance:
                                total -= filter[k]
                                k -= 1
                        if j < chr_indices[i + 1] - 4:
                            # remove downstream interactions outside of maxdistance
                            k = chr_indices[i + 1] - 1
                            stop = (j / 2) * 2 + 3
                            while k > stop and mids[k] - mids[j] >= maxdistance:
                                total -= filter[k]
                                k -= 1
                            stop = k + 1
                            # remove downstream interactions inside of mindistance
                            k = (j / 2) * 2 + 4
                            while k < stop and mids[k] - mids[j] < mindistance:
                                total -= filter[k]
                                k += 1
                        if j % 2 == 0:
                            # remove same fragment fend
                            total -= filter[j + 1]
                            # remove previous fragment fends if necessary
                            if j > chr_indices[i]:
                                total -= filter[j - 1]
                                if mids[j] - mids[j - 2] < mindistance:
                                    total -= filter[j - 2]
                            # remove next fragment fends if necessary
                            if j < chr_indices[i + 1] - 2:
                                total -= filter[j + 3]
                                if mids[j + 2] - mids[j] < mindistance:
                                    total -= filter[j + 2]
                        else:
                            # remove same fragment fend
                            total -= filter[j - 1]
                            # remove previous fragment fends if necessary
                            if j > chr_indices[i] + 2:
                                total -= filter[j - 3]
                                if mids[j] - mids[j - 2] < mindistance:
                                    total -= filter[j - 2]
                            # remove next fragment fends if necessary
                            if j < chr_indices[i + 1] - 2:
                                total -= filter[j + 1]
                                if mids[j + 2] - mids[j] < mindistance:
                                    total -= filter[j + 2]
                    else:
                        total = chr_total
                        if mindistance > 0:
                            total -= 1
                        # remove upstream interactions outside of maxdistance
                        k = chr_indices[i]
                        while k < j and mids[j] - mids[k] >= maxdistance:
                            total -= filter[k]
                            k += 1
                        stop = k - 1
                        # remove upstream interactions inside of mindistance
                        k = j - 1
                        while k > stop and mids[j] - mids[k] < mindistance:
                            total -= filter[k]
                            k -= 1
                        # remove downstream interactions outside of maxdistance
                        k = chr_indices[i + 1] - 1
                        while k > j and mids[k] - mids[j] >= maxdistance:
                            total -= filter[k]
                            k -= 1
                        stop = k + 1
                        # remove downstream interactions inside of mindistance
                        k = j + 1
                        while k < stop and mids[k] - mids[j] < mindistance:
                            total -= filter[k]
                            k += 1
                    interactions[j] += total
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def sum_weighted_indices(
        np.ndarray[DTYPE_int_t, ndim=1] indices0,
        np.ndarray[DTYPE_int_t, ndim=1] indices1,
        np.ndarray[DTYPE_t, ndim=1] weights,
        np.ndarray[DTYPE_int_t, ndim=1] int_weights,
        np.ndarray[DTYPE_64_t, ndim=1] sums):
    cdef long long int i
    cdef long long int num_indices = indices0.shape[0]
    with nogil:
        if not weights is None:
            for i in range(num_indices):
                sums[indices0[i]] += weights[i]
                sums[indices1[i]] += weights[i]
        elif not int_weights is None:
            for i in range(num_indices):
                sums[indices0[i]] += int_weights[i]
                sums[indices1[i]] += int_weights[i]
        else:
            for i in range(num_indices):
                sums[indices0[i]] += 1
                sums[indices1[i]] += 1
    return None




