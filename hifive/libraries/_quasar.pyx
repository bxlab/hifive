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
def filter_bins(
        np.ndarray[DTYPE_int64_t, ndim=2] raw,
        np.ndarray[DTYPE_int_t, ndim=1] sums,
        np.ndarray[DTYPE_int_t, ndim=1] valid,
        long long int width):
    cdef long long int i, bin1, bin2, span
    cdef long long int ndata = raw.shape[0]
    with nogil:
        for i in range(ndata):
            bin1 = raw[i, 0]
            bin2 = raw[i, 1]
            span = bin2 - bin1
            if valid[bin1] == 0 or valid[bin2] == 0 or span > width or span == 0:
                continue
            sums[bin1] += 1
            sums[bin2] += 1
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_bg_dist_norm(
        np.ndarray[DTYPE_int64_t, ndim=2] raw,
        np.ndarray[DTYPE_int_t, ndim=1] valid,
        np.ndarray[DTYPE_int64_t, ndim=2] temp,
        np.ndarray[DTYPE_64_t, ndim=1] bg,
        np.ndarray[DTYPE_64_t, ndim=2] dist,
        np.ndarray[DTYPE_64_t, ndim=2] corrs):
    cdef long long int i, j, bin1, bin2, span, count
    cdef long long int ndata = raw.shape[0]
    cdef long long int N = valid.shape[0]
    cdef long long int width = bg.shape[0]
    with nogil:
        for i in range(N):
            if valid[i] == 0:
                continue
            for j in range(i + 1, min(N, i + width + 1)):
                if valid[j] == 0:
                    continue
                span = j - i - 1
                temp[span, 1] += 1
                dist[i, span] = 1
        for i in range(ndata):
            bin1 = raw[i, 0]
            bin2 = raw[i, 1]
            span = bin2 - bin1 - 1
            if valid[bin1] == 0 or valid[bin2] == 0 or span >= width or span < 0:
                continue
            count = raw[i, 2]
            temp[span, 0] += count
            dist[bin1, span] = pow(1 + count, 0.5)
            corrs[bin1, width + span] = count
            corrs[bin2, width - span - 1] = count
        for i in range(width):
            if temp[i, 1] > 0 and temp[i, 0] > 0:
                bg[i] = temp[i, 0] / <double>temp[i, 1]
            else:
                bg[i] = 1.0
    return None


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_correlations(
        np.ndarray[DTYPE_64_t, ndim=2] norms,
        np.ndarray[DTYPE_int_t, ndim=1] vrows,
        np.ndarray[DTYPE_64_t, ndim=2] corrs,
        long long int start):
    cdef long long int i, j, k, pos, n, X, Y
    cdef double X1, X2, Y1, Y2, XY, temp, temp1
    cdef long long int N = vrows.shape[0]
    cdef long long int P = corrs.shape[0]
    cdef long long int M = corrs.shape[1]
    with nogil:
        for i in range(P):
            X = start + i
            if vrows[X] == 0:
                continue
            for Y in range(X + 1, min(N, X + M + 1)):
                if vrows[Y] == 0:
                    continue
                j = Y - start
                X1 = 0.0
                Y1 = 0.0
                X2 = 0.0
                Y2 = 0.0
                XY = 0.0
                n = 0
                for k in range(max(0, Y - M), X):
                    if vrows[k] == 0:
                        continue
                    pos = M - (X - k)
                    temp = norms[i, pos]
                    X1 += temp
                    X2 += temp * temp
                    pos = M - (Y - k)
                    temp1 = norms[j, pos]
                    Y1 += temp1
                    Y2 += temp1 * temp1
                    XY += temp * temp1
                    n += 1
                for k in range(X + 1, Y):
                    if vrows[k] == 0:
                        continue
                    pos = M + k - X - 1
                    temp = norms[i, pos]
                    X1 += temp
                    X2 += temp * temp
                    pos = M - (Y - k)
                    temp1 = norms[j, pos]
                    Y1 += temp1
                    Y2 += temp1 * temp1
                    XY += temp * temp1
                    n += 1
                for k in range(Y + 1, min(N, X + M + 1)):
                    if vrows[k] == 0:
                        continue
                    pos = M + k - X - 1
                    temp = norms[i, pos]
                    X1 += temp
                    X2 += temp * temp
                    pos = M + k - Y - 1
                    temp1 = norms[j, pos]
                    Y1 += temp1
                    Y2 += temp1 * temp1
                    XY += temp * temp1
                    n += 1
                if n < 3:
                    continue
                X1 /= n
                X2 = X2 / n - X1  * X1
                Y1 /= n
                Y2 = Y2 / n - Y1 * Y1
                if X2 <= 0.0 or Y2 <= 0.0:
                    continue
                XY /= n
                corrs[i, Y - X - 1] = (XY - X1 * Y1) / pow(X2 * Y2, 0.5)
    return None



