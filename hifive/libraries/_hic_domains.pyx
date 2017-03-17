# distutils: language = c++

"""These functions provide increased speed in handling HiC domain finding using HiFive.
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
    double log2(double x) nogil
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
def find_BIs(
        np.ndarray[DTYPE_t, ndim=3] data,
        np.ndarray[DTYPE_t, ndim=3] BIs,
        int minbin,
        int maxbin,
        int width):
    cdef long long int i, j, k, l, n, bi_count
    cdef double a0, a1, b0, b1, temp, bi_sum, bi2_sum
    cdef int num_bins = data.shape[0]
    cdef int max_bins = data.shape[1]
    with nogil:
        bi_sum = 0.0
        bi2_sum = 0.0
        bi_count = 0
        for i in range(1, num_bins - 1):
            n = 0
            a0 = 0
            a1 = 0
            b0 = 0
            b1 = 0
            for j in range(1, min(num_bins - i - 1, maxbin)):
                l = min(i, width)
                for k in range(l):
                    b0 += data[i - 1 - k, j + k, 0]
                    b1 += data[i - 1 - k, j + k, 1]
                l = min(j, width)
                for k in range(l):
                    a0 += data[i + k, j - k - 1, 0]
                    a1 += data[i + k, j - k - 1, 1]
                if a0 >= 10 and b0 >= 10:
                    #BIs[i, j + 1, 0] = log(a0 * b1 / (a1 * b0))
                    n += 1
                    if n > 1:
                        BIs[i, j + 1, 0] = (log(a0 * b1 / (a1 * b0)) + (j - 1) * BIs[i, j, 0]) / j
                    else:
                        BIs[i, j + 1, 0] = log(a0 * b1 / (a1 * b0)) / j
                elif n > 0:
                    BIs[i, j + 1, 0] = BIs[i, j, 0] * (j - 1.0) / j
        for i in range(2, num_bins):
            n = 0
            a0 = 0.0
            a1 = 0.0
            b0 = 0.0
            b1 = 0.0
            for j in range(1, min(i, maxbin)):
                l = min(j, width)
                for k in range(l):
                    a0 += data[i - j - 1, j - k - 1, 0]
                    a1 += data[i - j - 1, j - k - 1, 1]
                l = min(num_bins - i, width)
                for k in range(l):
                    b0 += data[i - j - 1, k + j, 0]
                    b1 += data[i - j - 1, k + j, 1]
                if a0 >= 10 and b0 >= 10:
                    #BIs[i - j - 1, j + 1, 1] = log(a0 * b1 / (a1 * b0))
                    n += 1
                    if n > 1:
                        BIs[i - j - 1, j + 1, 1] = (log(a0 * b1 / (a1 * b0)) + (j - 1) * BIs[i - j, j, 1]) / j
                    else:
                        BIs[i - j - 1, j + 1, 1] = log(a0 * b1 / (a1 * b0)) / j
                elif n > 0:
                    BIs[i - j - 1, j + 1, 1] = BIs[i - j, j, 1] * (j - 1.0) / j
        for i in range(num_bins):
            for j in range(minbin):
                BIs[i, j, 0] = 0
                BIs[i, j, 1] = 0
            for j in range(minbin, maxbin + 1):
                #if BIs[i, j, 0] > -Inf and BIs[i, j, 1] > -Inf:
                if BIs[i, j, 0] > 0 and BIs[i, j, 1] > 0:
                    BIs[i, j, 0] += BIs[i, j, 1]
                    BIs[i, j, 1] = 1
                    bi_sum += BIs[i, j, 0]
                    bi2_sum += BIs[i, j, 0] * BIs[i, j, 0]
                    bi_count += 1
                else:
                    BIs[i, j, 0] = 0
                    BIs[i, j, 1] = 0
        bi2_sum = pow(bi2_sum / (bi_count - 1), 0.5)
        for i in range(num_bins):
            for j in range(minbin, maxbin + 1):
                if BIs[i, j, 1] == 1:
                    if BIs[i, j, 0] > bi2_sum:
                        BIs[i, j, 0] /= bi2_sum
                    else:
                        BIs[i, j, 0] = 0
                        BIs[i, j, 1] = 0
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_BI_path(
        np.ndarray[DTYPE_t, ndim=3] bscores,
        np.ndarray[DTYPE_int_t, ndim=1] path,
        np.ndarray[DTYPE_64_t, ndim=1] scores,
        int minbin,
        int maxbin):
    cdef int i, j, k
    cdef double score
    cdef int num_bins = path.shape[0]
    with nogil:
        path[0] = 1
        for i in range(1, num_bins):
            scores[i] = scores[i - 1]
            path[i] = 1
            for j in range(minbin, min(maxbin + 1, i)):
                if bscores[i - j + 1, j - 1, 1] == 1:
                    score = scores[i - j] + pow(j, 0.5) * bscores[i - j + 1, j - 1, 0]
                    if score > scores[i]:
                        scores[i] = score
                        path[i] = j
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_arrowhead_transformation(
        np.ndarray[DTYPE_t, ndim=3] data,
        np.ndarray[DTYPE_t, ndim=2] score,
        int maxbin):
    cdef long long int i, j
    cdef double a, b
    cdef long long int num_bins = data.shape[0]
    cdef long long int max_bins = data.shape[1]
    with nogil:
        for i in range(1, num_bins - 1):
            for j in range(1, min(min(i, num_bins - i - 1), maxbin)):
                if data[i, j - 1, 1] > 0 and data[i - j, j - 1, 1] > 0:
                    a = data[i, j - 1, 0]
                    b = data[i - j, j - 1, 0]
                    if a > 0 or b > 0:
                        score[i, j - 1] = (b - a) / (a + b)
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_arrowhead_scores(
        np.ndarray[DTYPE_t, ndim=2] score,
        np.ndarray[DTYPE_t, ndim=3] sums,
        np.ndarray[DTYPE_t, ndim=3] signs,
        np.ndarray[DTYPE_t, ndim=4] variances,
        np.ndarray[DTYPE_t, ndim=2] domain_score,
        int minbin):
    cdef long long int h, i, j, k, l, dcount
    cdef double value, max_sum, max_sign, max_variance, v0, v1, e0, e1, s0, s1, dmean
    cdef long long int num_bins = score.shape[0]
    cdef long long int max_bins = score.shape[1]
    with nogil:
        max_sum = 0
        max_sign = 0
        max_variance = 0
        for i in range(minbin, num_bins - minbin):
            j = minbin
            for k in range(j / 2):
                for l in range(k, j - k):
                    value = score[i + k, l]
                    if value != 0:
                        sums[i, j - 1, 0] += value
                        if value > 0:
                            signs[i, j - 1, 0] += 1
                        else:
                            signs[i, j - 1, 0] -= 1
                        variances[i, j - 1, 0, 0] += 1
                        variances[i, j - 1, 0, 1] += value * value
            for k in range(j - j / 2, j + 1):
                for l in range(j - k, k):
                    value = score[i + k, l]
                    if value != 0:
                        sums[i, j - 1, 1] += value
                        if value > 0:
                            signs[i, j - 1, 1] += 1
                        else:
                            signs[i, j - 1, 1] -= 1
                        variances[i, j - 1, 1, 0] += 1
                        variances[i, j - 1, 1, 1] += value * value
        for j in range(minbin + 1, max_bins):
            for i in range(num_bins - j - 1):
                sums[i, j - 1, 0] = sums[i, j - 2, 0]
                signs[i, j - 1, 0] = signs[i, j - 2, 0]
                variances[i, j - 1, 0, 0] = variances[i, j - 2, 0, 0]
                variances[i, j - 1, 0, 1] = variances[i, j - 2, 0, 1]
                for k in range(j / 2):
                    l = j - k - 1
                    value = score[i + k, l]
                    if value != 0:
                        sums[i, j - 1, 0] += value
                        if value > 0:
                            signs[i, j - 1, 0] += 1
                        else:
                            signs[i, j - 1, 0] -= 1
                        variances[i, j - 1, 0, 0] += 1
                        variances[i, j - 1, 0, 1] += value * value
                sums[i, j - 1, 1] = sums[i + 1, j - 2, 1]
                signs[i, j - 1, 1] = signs[i + 1, j - 2, 1]
                variances[i, j - 1, 1, 0] = variances[i + 1, j - 2, 1, 0]
                variances[i, j - 1, 1, 1] = variances[i + 1, j - 2, 1, 1]
                for k in range(j - j / 2, j + 1):
                    l = k - 1
                    value = score[i + k, l]
                    if value != 0:
                        sums[i, j - 1, 1] += value
                        if value > 0:
                            signs[i, j - 1, 1] += 1
                        else:
                            signs[i, j - 1, 1] -= 1
                        variances[i, j - 1, 1, 0] += 1
                        variances[i, j - 1, 1, 1] += value * value
        for i in range(minbin, num_bins - minbin):
            for j in range(minbin - 1, min(min(i, num_bins - i - 1), max_bins) - 1):
                if variances[i, j, 0, 0] > 1 and variances[i, j, 1, 0] > 1:
                    e0 = sums[i, j, 0] / variances[i, j, 0, 0]
                    e1 = sums[i, j, 1] / variances[i, j, 1, 0]
                    v0 = variances[i, j, 0, 1] / variances[i, j, 0, 0] - e0 * e0
                    v1 = variances[i, j, 1, 1] / variances[i, j, 1, 0] - e1 * e1
                    s0 = signs[i, j, 0] / variances[i, j, 0, 0]
                    s1 = signs[i, j, 1] / variances[i, j, 1, 0]
                    #if s0 > -0.5 or s1 < 0.5:
                    #    continue
                    sums[i, j, 0] = max(0, e1) - min(0, e0)
                    #sums[i, j, 0] = pow(max(0, e1) * max(0, -e0), 0.5)
                    if sums[i, j, 0] > max_sum:
                        max_sum = sums[i, j, 0]
                    #sums[i, j, 1] = s1 - s0
                    sums[i, j, 1] = max(0, s1) * max(0, -s0)
                    if sums[i, j, 1] > max_sign:
                        max_sign = sums[i, j, 1]
                    variances[i, j, 0, 0] = v0 + v1
                    #variances[i, j, 0, 0] = pow(v0 * v1, 0.5)
                    if variances[i, j, 0, 0] > max_variance:
                        max_variance = variances[i, j, 0, 0]
                    domain_score[i, j] = 1
        for i in range(minbin, num_bins - minbin):
            for j in range(minbin - 1, min(min(i, num_bins - i - 1), max_bins) - 1):
                if domain_score[i, j] == 1:
                    e0 = sums[i, j, 0] / max_sum
                    e1 = sums[i, j, 1] / max_sign
                    v0 = variances[i, j, 0, 0] / max_variance
                    #if v0 < 0.2:
                    domain_score[i, j] = (sums[i, j, 1] / max_sign) * ((sums[i, j, 0] / max_sum) + (variances[i, j, 0, 0] / max_variance))
                    dmean += domain_score[i, j]
                    dcount += 1
        dmean /= dcount * 2
        for i in range(minbin, num_bins - minbin):
            for j in range(minbin - 1, min(min(i, num_bins - i - 1), max_bins) - 1):
                if domain_score[i, j] < dmean:
                    domain_score[i, j] = 0 
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_arrowhead_path(
        np.ndarray[DTYPE_t, ndim=2] dscores,
        np.ndarray[DTYPE_int_t, ndim=1] path,
        np.ndarray[DTYPE_64_t, ndim=1] scores,
        int minbin,
        int maxbin):
    cdef int i, j, k
    cdef double score
    cdef int num_bins = path.shape[0]
    with nogil:
        path[0] = 1
        for i in range(1, num_bins):
            scores[i] = scores[i - 1]
            path[i] = 1
            for j in range(minbin, min(maxbin, i)):
                score = scores[i - j] + pow(j, 0.5) * dscores[i - j, j - 1]
                if score > scores[i]:
                    scores[i] = score
                    path[i] = j
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_band_score(
        np.ndarray[DTYPE_t, ndim=3] hm,
        np.ndarray[DTYPE_t, ndim=2] scores,
        int minband,
        int maxband,
        int band):
    cdef int i, j, k
    cdef double observed, expected
    cdef int num_bins = hm.shape[0]
    cdef int hminband = minband / 2
    cdef int hmaxband = maxband / 2
    cdef int width = hmaxband - hminband
    with nogil:
        """
        observed = 0.0
        expected = 0.0
        for j in range(width):
            for k in range(maxband - width, maxband):
                observed += hm[j, k - j - 1, 0]
                expected += hm[j, k - j - 1, 1]
        if observed > 0.0:#expected > 0.0:
            scores[hmaxband, band] = log(observed / expected)
        for i in range(hmaxband + 1, num_bins - hmaxband + 1):
            for j in range(i - hmaxband, i - hminband):
                observed -= hm[j - 1, i + hminband - j - 1, 0]
                observed += hm[j, i + hmaxband - j - 2, 0]
                expected -= hm[j - 1, i + hminband - j - 1, 1]
                expected += hm[j, i + hmaxband - j - 2, 1]
            for j in range(i + hminband, i + hmaxband - 1):
                observed -= hm[i - hmaxband - 1, j - i + hmaxband, 0]
                expected -= hm[i - hmaxband - 1, j - i + hmaxband, 1]
                observed += hm[i - hminband - 1, j - i + hminband, 0]
                expected += hm[i - hminband - 1, j - i + hminband, 1]
            if observed > 0.0:#expected > 0.0:
                scores[i, band] = log(observed / expected)
        """
        for i in range(hmaxband, num_bins - hmaxband + 1):
            observed = 0.0
            expected = 0.0
            for j in range(i - hmaxband, i - hminband):
                for k in range(i + hminband, i + hmaxband):
                    observed += hm[j, k - j - 1, 0]
                    expected += hm[j, k - j - 1, 1]
            if observed > 0.0:#expected > 0.0:
                scores[i, band] = log(observed / expected)
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def dynamically_bin(
        np.ndarray[DTYPE_int_t, ndim=2] unbinned_c not None,
        np.ndarray[DTYPE_t, ndim=2] unbinned_e not None,
        np.ndarray[DTYPE_int_t, ndim=2] binned_c not None,
        np.ndarray[DTYPE_t, ndim=2] binned_e not None,
        np.ndarray[DTYPE_int_t, ndim=1] valid not None,
        np.ndarray[DTYPE_int_t, ndim=1] X_coords not None,
        np.ndarray[DTYPE_int_t, ndim=1] Y_coords not None,
        int minobservations):
    cdef long long int X, Y, lX, elX, uX, euX, eX, lY, elY, uY, euY, eY, i, j, k, iteration
    cdef long long int num_bins = binned_c.shape[0]
    cdef long long int num_coords = X_coords.shape[0]
    with nogil:
        for i in range(num_coords):
            X = X_coords[i]
            Y = Y_coords[i]
            if valid[X] == 0 or valid[Y] == 0:
                continue
            if unbinned_c[X, Y] >= minobservations:
                binned_c[X, Y] = unbinned_c[X, Y]
                binned_c[Y, X] = unbinned_c[X, Y]
                binned_e[X, Y] = unbinned_e[X, Y]
                binned_e[Y, X] = unbinned_e[X, Y]
                lX = X - 1
                uX = X + 1
                lY = Y - 1
                uY = Y + 1
                continue
            elif (i > 0 and Y == Y_coords[i - 1] + 1 and X == X_coords[i - 1] and
                  binned_c[X, Y - 1] > 0 and lX - uX > 3):
                lX += 1
                uX -= 1
                lY += 2
            else:
                lX = X - 1
                uX = X + 1
                lY = Y - 1
                uY = Y + 1
            elX = max(lX, 0)
            euX = min(uX, num_bins - 1)
            elY = max(lY, 0)
            euY = min(uY, num_bins - 1)
            for j in range(max(lX + 1, 0), euX):
                for k in range(max(j, lY + 1), euY):
                    binned_c[X, Y] += unbinned_c[j, k]
                    binned_e[X, Y] += unbinned_e[j, k]
            iteration  = 1
            while binned_c[X, Y] < minobservations:
                if lX >= 0 and valid[lX] == 1:
                    for j in range(elY, euY + 1):
                        binned_c[X, Y] += unbinned_c[lX, j]
                        binned_e[X, Y] += unbinned_e[lX, j]
                if uY < num_bins and valid[uY] == 1:
                    for j in range(elX + 1, euX + 1):
                        binned_c[X, Y] += unbinned_c[j, uY]
                        binned_e[X, Y] += unbinned_e[j, uY]
                if uX < uY and uX < num_bins and valid[uX] == 1:
                    eY = max(uX, lY)
                    for j in range(eY, euY):
                        binned_c[X, Y] += unbinned_c[uX, j]
                        binned_e[X, Y] += unbinned_e[uX, j]
                if lY > lX and lY >= 0 and valid[lY] == 1:
                    eX = min(uX - 1, lY)
                    for j in range(elX + 1, eX + 1):
                        binned_c[X, Y] += unbinned_c[j, lY]
                        binned_e[X, Y] += unbinned_e[j, lY]
                lX -= 1
                elX = max(lX, 0)
                uX += 1
                euX = min(uX, num_bins - 1)
                lY -= 1
                elY = max(lY, 0)
                uY += 1
                euY = min(uY, num_bins - 1)
                iteration += 1
            binned_c[Y, X] = binned_c[X, Y]
            binned_e[Y, X] = binned_e[X, Y]
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_distance_parameters(
        np.ndarray[DTYPE_int_t, ndim=1] counts not None,
        np.ndarray[DTYPE_t, ndim=1] expected not None,
        np.ndarray[DTYPE_int_t, ndim=1] valid not None,
        np.ndarray[DTYPE_int_t, ndim=1] states not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices1 not None,
        np.ndarray[DTYPE_int_t, ndim=2] dcounts not None,
        np.ndarray[DTYPE_64_t, ndim=2] dexpected not None,
        np.ndarray[DTYPE_int_t, ndim=2] binsizes not None):
    cdef long long int i, j, x, y, start, dist, index, bin1, bin2
    cdef long long int num_bins = indices0.shape[0]
    with nogil:
        for i in range(num_bins):
            bin1 = indices0[i]
            if valid[bin1] <= 0:
                continue
            bin2 = indices1[i]
            if valid[bin2] <= 0:
                continue
            if expected[i] == 0.0:
                continue
            if states[bin1] == states[bin2]:
                if states[bin1] == 1:
                    index = 1
                else:
                    index = 0
            else:
                index = 2
            dist = bin2 - bin1
            dcounts[dist, index] += counts[i]
            dexpected[dist, index] += expected[i]
            binsizes[dist, index] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def optimize_compartment_bounds(
        np.ndarray[DTYPE_int_t, ndim=1] counts not None,
        np.ndarray[DTYPE_t, ndim=1] expected not None,
        np.ndarray[DTYPE_int_t, ndim=1] valid not None,
        np.ndarray[DTYPE_int_t, ndim=1] states not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices0 not None,
        np.ndarray[DTYPE_int_t, ndim=1] indices1 not None,
        np.ndarray[DTYPE_64_t, ndim=1] paramsA not None,
        np.ndarray[DTYPE_64_t, ndim=1] paramsB not None,
        np.ndarray[DTYPE_64_t, ndim=1] paramsAB not None,
        np.ndarray[DTYPE_64_t, ndim=2] scores not None):
    cdef long long int i, X, Y, dist
    cdef double muA, muB, lambdaA, lambdaB, temp1, temp2
    cdef long long int num_indices = indices0.shape[0]
    with nogil:
        for i in range(num_indices):
            X = indices0[i]
            if valid[X] == 0:
                continue
            Y = indices1[i]
            if valid[Y] == 0:
                continue
            dist = Y - X
            if states[Y] == 1:
                muA = paramsAB[dist]
                muB = paramsB[dist]
            else:
                muA = paramsA[dist]
                muB = paramsAB[dist]
            lambdaA = muA * expected[i]
            lambdaB = muB * expected[i]
            temp1 = counts[i] * log(lambdaA) - lambdaA
            temp2 = counts[i] * log(lambdaB) - lambdaB
            scores[X, 1] += temp1
            scores[X, 0] += temp2
            if states[X] == states[Y]:
                scores[Y, 1] += temp1
                scores[Y, 0] += temp2
            else:
                if states[X] == 1:
                    muA = paramsAB[dist]
                    muB = paramsB[dist]
                else:
                    muA = paramsA[dist]
                    muB = paramsAB[dist]
                lambdaA = muA * expected[i]
                lambdaB = muB * expected[i]
                scores[Y, 1] += counts[i] * log(lambdaA) - lambdaA
                scores[Y, 0] += counts[i] * log(lambdaB) - lambdaB
    return None
















