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
def find_correlations(
        np.ndarray[DTYPE_64_t, ndim=2] norms,
        np.ndarray[DTYPE_int_t, ndim=1] vrows,
        np.ndarray[DTYPE_64_t, ndim=2] corrs,
        long long int start,
        long long int stop):
    cdef long long int i, n, X, Y
    cdef double X1, X2, Y1, Y2, XY, temp, temp1
    cdef long long int N = norms.shape[0]
    cdef long long int M = corrs.shape[1]
    with nogil:
        for X in range(start, stop):
            if vrows[X] == 0:
                continue
            for Y in range(X + 1, min(N, X + M + 1)):
                if vrows[Y] == 0:
                    continue
                X1 = 0.0
                Y1 = 0.0
                X2 = 0.0
                Y2 = 0.0
                XY = 0.0
                n = 0
                for i in range(max(0, Y - M), X):
                    if vrows[i] == 0:
                        continue
                    temp = norms[i, X]
                    X1 += temp
                    X2 += temp * temp
                    temp1 = norms[i, Y]
                    Y1 += temp1
                    Y2 += temp1 * temp1
                    XY += temp * temp1
                    n += 1
                for i in range(X + 1, Y):
                    if vrows[i] == 0:
                        continue
                    temp = norms[X, i]
                    X1 += temp
                    X2 += temp * temp
                    temp1 = norms[i, Y]
                    Y1 += temp1
                    Y2 += temp1 * temp1
                    XY += temp * temp1
                    n += 1
                for i in range(Y + 1, min(N, X + M + 1)):
                    if vrows[i] == 0:
                        continue
                    temp = norms[X, i]
                    X1 += temp
                    X2 += temp * temp
                    temp1 = norms[Y, i]
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
                corrs[X - start, Y - X - 1] = (XY - X1 * Y1) / pow(X2 * Y2, 0.5)
    return None




