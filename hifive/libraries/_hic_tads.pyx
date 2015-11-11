# distutils: language = c++

"""These functions provide increased speed in handling TAD finding using HiFive.
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
def find_betadeltas(
        np.ndarray[DTYPE_t, ndim=3] data,
        np.ndarray[DTYPE_t, ndim=2] betas,
        np.ndarray[DTYPE_t, ndim=2] deltas,
        np.ndarray[DTYPE_t, ndim=2] fits,
        np.ndarray[DTYPE_t, ndim=2] errors,
        int maxbin):
    cdef long long int i, j, k, l, n
    cdef double x, y, x2, xy, temp, temp1, beta, delta, fit, error
    cdef int num_bins = data.shape[0]
    cdef int minsize = 5
    with nogil:
        for i in range(num_bins - minsize + 1):
            l = minsize
            n = 0
            x = 0.0
            x2 = 0.0
            y = 0.0
            xy = 0.0
            error = 0.0
            fit = 0.0
            delta = 0.0
            beta = 0.0
            for j in range(i, i + l):
                for k in range(j + 1, i + l + 1):
                    if data[j, k, 1] > 0.0:
                        n += 1
                        temp = k - j
                        x += temp
                        x2 += temp * temp
                        temp1 = data[j, k, 0] / data[j, k, 1]
                        y += temp1
                        xy += temp * temp1
                        error += pow(data[j, k, 0] - data[j, k, 1], 2)
            if n >= 3:
                delta = (xy - x * y / n) / (x2 - x * x / n)
                beta = (y - delta * x) / n
                for j in range(i, i + l -  1):
                    for k in range(j + 1, i + l):
                        if data[j, k, 1] > 0.0:
                            fit += pow(data[j, k, 1] * (beta + delta * (k - j)) - data[j, k, 0], 2)
            betas[i, i + l] = beta
            deltas[i, i + l] = delta
            if delta > 0.0:
                fits[i, i + l] = fit - error
            else:
                fits[i, i + minsize] = Inf
            errors[i, i + l] = error
            for l in range(minsize + 1, min(maxbin, num_bins - i + 1)):
                k = i + l
                for j in range(i, i + l):
                    if data[j, k, 1] > 0.0:
                        n += 1
                        temp = k - j
                        x += temp
                        x2 += temp * temp
                        temp1 = data[j, k, 0] / data[j, k, 1]
                        y += temp1
                        xy += temp * temp1
                        error += pow(data[j, k, 0] - data[j, k, 1], 2)
                if n >= 3:
                    delta = (xy - x * y / n) / (x2 - x * x / n)
                    beta = (y - delta * x) / n
                    fit = 0.0
                    for j in range(i, i + l -  1):
                        for k in range(j + 1, i + l):
                            if data[j, k, 1] > 0.0:
                                fit += pow(data[j, k, 1] * (beta + delta * (k - j)) - data[j, k, 0], 2)
                betas[i, i + l] = beta
                deltas[i, i + l] = delta
                if delta > 0.0:
                    fits[i, i + l] = fit - error
                else:
                    fits[i, i + minsize] = Inf
    return None

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
def find_initial_TAD_std_params(
        np.ndarray[DTYPE_t, ndim=3] data,
        np.ndarray[DTYPE_t, ndim=3] BI_score,
        np.ndarray[DTYPE_t, ndim=2] scores,
        np.ndarray[DTYPE_t, ndim=3] params,
        int maxbin,
        int minsize,
        double gamma):
    cdef long long int i, j, k, l, m, n
    cdef double x, x2, temp
    cdef int num_bins = data.shape[0]
    with nogil:
        for i in range(num_bins - minsize + 1):
            l = minsize
            m = 0
            n = 0
            x = 0.0
            x2 = 0.0
            for j in range(i, i + l - 1):
                for k in range(i + l - j - 1):
                    if data[j, k, 1] > 0.0:
                        n += 1
                        temp = data[j, k, 0]
                        x += temp
                        x2 += temp * temp
            params[i, l, 0] = n
            params[i, l, 1] = x
            params[i, l, 2] = x2
            if n >= 3 and BI_score[i, l, 1] > 0:
                #scores[i, m] = BI_score[i, m] / pow(x2 / n - (x / n) * (x / n), 0.5)
                scores[i, l] = pow(x2 / n - (x / n) * (x / n), 0.5) - BI_score[i, l, 0] * gamma
            else:
                scores[i, l] = Inf
            for l in range(minsize + 1, min(maxbin + 1, num_bins - i + 1)):
                k = i + l - 1
                m = l - minsize
                for j in range(i, i + l - 1):
                    if data[j, k - j - 1, 1] > 0.0:
                        n += 1
                        temp = data[j, k - j - 1, 0]
                        x += temp
                        x2 += temp * temp
                params[i, l, 0] = n
                params[i, l, 1] = x
                params[i, l, 2] = x2
                if n >= 3 and BI_score[i, l, 1] > 0:
                    #scores[i, m] = BI_score[i, m] / pow(x2 / n - (x / n) * (x / n), 0.5)
                    scores[i, l] = pow(x2 / n - (x / n) * (x / n), 0.5) - BI_score[i, l, 0] * gamma
                else:
                    scores[i, l] = Inf
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_TAD_path(
        np.ndarray[DTYPE_t, ndim=2] scores,
        np.ndarray[DTYPE_int_t, ndim=2] paths,
        np.ndarray[DTYPE_t, ndim=3] path_scores,
        np.ndarray[DTYPE_int_t, ndim=1] final_path,
        int minbin,
        int maxbin):
    cdef int i, j, k, index
    cdef int num_bins = scores.shape[0]
    cdef int zeroindex = maxbin - minbin + 1
    cdef double score, score1, temp_score
    with nogil:
        path_scores[num_bins, zeroindex, 0] = 0.0
        path_scores[num_bins, zeroindex, 1] = 0.0
        for i in range(1, num_bins + 1):
            index = num_bins - i
            score = path_scores[index + 1, zeroindex, 0] #/ max(1, path_scores[index + 1, zeroindex, 1])
            path_scores[index, zeroindex, 0] = path_scores[index + 1, zeroindex, 0]
            path_scores[index, zeroindex, 1] = path_scores[index + 1, zeroindex, 1]
            paths[index, zeroindex] = zeroindex
            for j in range(maxbin - minbin + 1):
                temp_score = path_scores[index + 1, j, 0] #/ max(1, path_scores[index + 1, j, 1])
                if temp_score < score:
                    score = temp_score
                    path_scores[index, zeroindex, 0] = path_scores[index + 1, j, 0]
                    path_scores[index, zeroindex, 1] = path_scores[index + 1, j, 1]
                    paths[index, zeroindex] = j
                #path_scores[index, zeroindex, 0] = Inf
                if j + minbin <= i:
                    score1 = (path_scores[index + j + minbin, zeroindex, 0] + scores[index, j]) #/ max( 1,
                    #    path_scores[index + j + minbin, zeroindex, 1] + j + minbin)
                    path_scores[index, j, 0] = path_scores[index + j + minbin, zeroindex, 0] + scores[index, j] #* j
                    path_scores[index, j, 1] = path_scores[index + j + minbin, zeroindex, 1] + j + minbin
                    paths[index, j] = zeroindex
                    for k in range(maxbin - minbin + 1):
                        temp_score = (path_scores[index + j + minbin, k, 0] + scores[index, j]) #/ max( 1,
                        #    path_scores[index + j + minbin, k, 1] + j + minbin)
                        if temp_score < score1:
                            score1 = temp_score
                            path_scores[index, j, 0] = path_scores[index + j + minbin, k, 0] + scores[index, j] #* j
                            path_scores[index, j, 1] = path_scores[index + j + minbin, k, 1] + j + minbin
                            paths[index, j] = k
        score = Inf
        index = -1
        while score == Inf:
            index += 1
            k = zeroindex
            for i in range(maxbin - minbin + 1):
                #if path_scores[0, i, 0] / max(1, path_scores[0, i, 1]) < path_scores[0, k, 0] / max(1,
                #    path_scores[0, k, 1]):
                if path_scores[index, i, 0] < path_scores[index, k, 0]:
                    k = i
            score = path_scores[index, k, 0]
        while index < num_bins:
            if k != zeroindex:
                final_path[index] = minbin + k
                index += k + minbin
            else:
                index += 1
            if index < num_bins:
                k = paths[index, k]
    return None



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_TAD_subparts(
        np.ndarray[DTYPE_t, ndim=3] scores,
        np.ndarray[DTYPE_t, ndim=4] sub_params,
        np.ndarray[DTYPE_t, ndim=2] BIs,
        np.ndarray[DTYPE_t, ndim=3] params,
        int minbins,
        int TADsize,
        int treesize,
        double gamma):
    cdef int i, j, k, l, n
    cdef double temp, left, right
    cdef int num_bins = scores.shape[0]
    cdef int maxbins = scores.shape[1] - minbins
    cdef int maxtree = scores.shape[2]
    with nogil:
        n = TADsize - minbins
        if treesize == 1:
            for i in range(num_bins - TADsize + 1):
                if params[i, n, 0] >= 3:
                    temp = params[i, n, 1] / params[i, n, 0]
                    scores[i, n, 0] = (pow(BIs[i, 0] * BIs[i + TADsize - 1, 1], gamma) -
                        ((params[i, n, 2] / params[i, n, 0]) - temp * temp))
                    sub_params[i, n, 0, 0] = params[i, n, 0]
                    sub_params[i, n, 0, 1] = params[i, n, 1]
                    sub_params[i, n, 0, 2] = params[i, n, 2]
                else:
                    scores[i, n, 0] = -Inf
        """
        else:
            # for each possible TAD starting position
            for i in range(num_bins - TADsize + 1):
                # for each possible TAD subdivision
                best_subdiv_score = -Inf
                best_subdiv = -1
                best_subcombo = -1
                for j in range(i + 1, i + TADsize):
                    # for each combination of values from left and right subdivisions summing to treesize
                    best_combo_score = -Inf
                    best_combo = -1
                    best_eff_size = -1
                    best_eff_x = -1
                    best_eff_x2 = -1
                    for k in range(treesize + 1):
                        # if left side is too small or is invalid, skip
                        eff_size = params[i, n, 0]
                        eff_x = params[i, n, 1]
                        eff_x2 = params[i, n, 2]
                        index = j - i + 1 - minbins
                        if k == 0:
                            left = 0.0
                        elif index < 0:
                            continue
                        else:
                            left = scores[i, index, k - 1]
                            eff_size -= sub_params[i, index, k - 1, 0]
                            eff_x -= sub_params[i, index, k - 1, 1]
                            eff_x2 -= sub_params[i, index, k - 1, 2]
                        # if right side is too small or is invalid, skip
                        index = i + TADsize - j - minbins
                        if k == treesize:
                            right = 0.0
                        elif index < 0:
                            continue
                        else:
                            right = scores[j, index, k - 1]
                            eff_size -= sub_params[j, index, k - 1, 0]
                            eff_x -= sub_params[j, index, k - 1, 1]
                            eff_x2 -= sub_params[j, index, k - 1, 2]
                        # if possible, find score, otherwise mark invalid
                        if left == -Inf or right == -Inf or eff_size < 3:
                            continue
                        temp = eff_x / eff_size
                        combo_score = (pow(BIs[i, 0] * BIs[i + TADsize - 1, 1], gamma) -
                            (eff_x2 / eff_size - temp * temp)) + left + right
                        # determine if left-right tree sizes is superior
                        if combo_score > best_combo_score:
                            best_combo_score = combo_score
                            best_combo = k
                            best_eff_size = params[i, n, 0] - eff_size
                            best_eff_x = params[i, n, 1] - eff_x
                            best_eff_x2 = params[i, n, 2] - eff_x2
                    # determine in subdivision point is superior
                    if best_combo_score < best_subdiv_score:
                        best_subdiv_score = best_combo_score
                        best_subcombo = best_combo
                        best_subdiv = j

                # if valid subtree found, record score and parameters
                if best_subdiv == -1:
                    scores[i, n, treesize - 1] = -Inf
                #else:
                #    sub_params[i, n, k - 1, 1]
        """

    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def build_TAD_trees(
        np.ndarray[DTYPE_t, ndim=3] data,
        np.ndarray[DTYPE_t, ndim=2] fits,
        np.ndarray[DTYPE_t, ndim=2] deltas,
        np.ndarray[DTYPE_t, ndim=2] betas,
        np.ndarray[DTYPE_t, ndim=2] errors,
        np.ndarray[DTYPE_t, ndim=3] scores,
        np.ndarray[DTYPE_t, ndim=2] local_scores,
        np.ndarray[DTYPE_int_t, ndim=2] local_k,
        np.ndarray[DTYPE_int_t, ndim=2] local_t,
        np.ndarray[DTYPE_t, ndim=2] options,
        np.ndarray[DTYPE_t, ndim=5] intervals,
        int maxbins):
    cdef long long int i, j, k, l, m, n, p, t, tt, local_maxtree, best_m, best_t, pos, newpos
    cdef double old_score, old_delta, old_beta, best
    cdef int num_bins = data.shape[0]
    cdef int height = min(num_bins, maxbins)
    cdef int maxtree = scores.shape[2]
    with nogil:
        for i in range(2, height):
            for j in range(num_bins - i):
                k = i + j
                if fits[j, k] > 0:
                    continue
                if deltas[j, k] < 0:
                    for l in range(maxtree):
                        scores[j, k, l] = Inf
                else:
                    scores[j, k, 0] = fits[j, k]
                    local_maxtree = min(maxtree, i - 1)
                    for l in range(2, i + 1):
                        for t in range(1, min(maxtree, l - 1)):
                            options[0, 0] = local_scores[l - 1, t]
                            for m in range(2, l + 1):
                                for tt in range(min(m - 2, t)):
                                    old_delta = deltas[j, k]
                                    old_beta = betas[j, k]
                                    if old_delta > deltas[j + l - m, j + l]:
                                        options[m, tt] = Inf
                                    else:
                                        old_score = 0.0
                                        for n in range(m - 1):
                                            x = j + l - m + n
                                            for p in range(n + 1, m):
                                                y = j + l - m + p
                                                old_score += pow((old_delta * (p - n) + old_beta) *
                                                              data[x, y, 1] - data[x, y, 0], 2)
                                        options[m, tt] = (local_scores[l - m, t - tt - 1] +
                                                          scores[j + l - m, j + l, tt] - (old_score -
                                                          errors[j + l - m, j + l]))
                            best = Inf
                            for m in range(l + 1):
                                for tt in range(t):
                                    if options[m, tt] < best:
                                        best = options[m, tt]
                                        best_m = m
                                        best_t = tt
                            if best_m == 0 and best_t == 0:
                                local_k[l, t] = -1
                                local_t[l, t] = t
                            else:
                                local_k[l, t] = best_m
                                local_t[l, t] = best_t
                    for t in range(1, local_maxtree):
                        scores[j, k, t] = local_scores[i, t] + fits[j, k]
                        pos = i
                        tt = t
                        n = 0
                        while True:
                            if pos <= 2 or local_k[pos, tt] == 0:
                                break
                            if local_k[pos, tt] == -1:
                                pos -= 1
                            else:
                                newpos = pos - local_k[pos, tt]
                                intervals[j, k, t, n, 0] = j + newpos
                                intervals[j, k, t, n, 1] = j + pos
                                intervals[j, k, t, n, 2] = local_t[pos, tt]
                                tt -= local_t[pos, tt] + 1
                                pos = newpos
                                n += 1
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



















