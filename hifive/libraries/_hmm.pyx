# distutils: language = c++

"""These functions provide increased speed in HMM processes.
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
cdef double PI = numpy.pi

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
def find_probabilities(
        np.ndarray[DTYPE_64_t, ndim=1] observations,
        np.ndarray[DTYPE_64_t, ndim=3] distributions,
        np.ndarray[DTYPE_64_t, ndim=3] probabilities):
    cdef int i, j, k
    cdef double temp, prob, sigma2, scalar
    cdef int num_states = probabilities.shape[0]
    cdef int num_dists = probabilities.shape[1] - 1
    cdef int num_obs = probabilities.shape[2]
    with nogil:
        for i in range(num_states):
            for j in range(num_dists):
                scalar = distributions[i, j, 0] * pow(2.0 * PI * distributions[i, j, 2], -0.5)
                sigma2 = 1.0 / (2.0 * distributions[i, j, 2])
                for k in range(num_obs):
                    temp = observations[k] - distributions[i, j, 1]
                    prob = scalar * exp(-temp * temp * sigma2)
                    probabilities[i, j, k] = prob
                    probabilities[i, num_dists, k] += prob
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_alphas(
        np.ndarray[DTYPE_64_t, ndim=3] probabilities,
        np.ndarray[DTYPE_64_t, ndim=1] pi,
        np.ndarray[DTYPE_64_t, ndim=2] transitions,
        np.ndarray[DTYPE_64_t, ndim=2] alphas,
        np.ndarray[DTYPE_64_t, ndim=1] scalars):
    cdef int i, j, k
    cdef int num_states = probabilities.shape[0]
    cdef int num_dists = probabilities.shape[1] - 1
    cdef int num_obs = probabilities.shape[2]
    with nogil:
        for i in range(num_states):
            alphas[i, 0] = probabilities[i, num_dists, 0] * pi[i]
            scalars[0] += alphas[i, 0]
        for i in range(num_states):
            alphas[i, 0] /= scalars[0]
        for k in range(1, num_obs):
            for i in range(num_states):
                for j in range(num_states):
                    alphas[i, k] += alphas[j, k - 1] * transitions[j, i]
                alphas[i, k] *= probabilities[i, num_dists, k]
                scalars[k] += alphas[i, k]
            for i in range(num_states):
                alphas[i, k] /= scalars[k]
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_betas(
        np.ndarray[DTYPE_64_t, ndim=3] probabilities,
        np.ndarray[DTYPE_64_t, ndim=2] transitions,
        np.ndarray[DTYPE_64_t, ndim=1] scalars,
        np.ndarray[DTYPE_64_t, ndim=2] betas):
    cdef int i, j, k, l
    cdef double total, temp
    cdef int num_states = probabilities.shape[0]
    cdef int num_dists = probabilities.shape[1] - 1
    cdef int num_obs = probabilities.shape[2]
    with nogil:
        total = 0.0
        for i in range(num_states):
            betas[i, num_obs - 1] = 1.0
        for l in range(1, num_obs):
            k = num_obs - l - 1
            for i in range(num_states):
                for j in range(num_states):
                    betas[i, k] += betas[j, k + 1] * probabilities[j, num_dists, k + 1] * transitions[i, j]
                betas[i, k] /= scalars[k + 1]
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_etas(
        np.ndarray[DTYPE_64_t, ndim=3] probabilities,
        np.ndarray[DTYPE_64_t, ndim=2] transitions,
        np.ndarray[DTYPE_64_t, ndim=2] alphas,
        np.ndarray[DTYPE_64_t, ndim=2] betas,
        np.ndarray[DTYPE_64_t, ndim=3] etas):
    cdef int i, j, k, l
    cdef double total, temp
    cdef int num_states = probabilities.shape[0]
    cdef int num_dists = probabilities.shape[1] - 1
    cdef int num_obs = probabilities.shape[2]
    with nogil:
        for k in range(num_obs - 1):
            total = 0.0
            for i in range(num_states):
                for j in range(num_states):
                    etas[i, j, k] = (alphas[i, k] * transitions[i, j] * probabilities[j, num_dists, k + 1] *
                                     betas[j, k + 1])
                    total += etas[i, j, k]
            for i in range(num_states):
                for j in range(num_states):
                    etas[i, j, k] /= total
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_gammas(
        np.ndarray[DTYPE_64_t, ndim=3] probabilities,
        np.ndarray[DTYPE_64_t, ndim=2] alphas,
        np.ndarray[DTYPE_64_t, ndim=2] betas,
        np.ndarray[DTYPE_64_t, ndim=3] gammas):
    cdef int i, j, k, l
    cdef double total, temp
    cdef int num_states = probabilities.shape[0]
    cdef int num_dists = probabilities.shape[1] - 1
    cdef int num_obs = probabilities.shape[2]
    with nogil:
        for k in range(num_obs):
            for i in range(num_states):
                temp = alphas[i, k] * betas[i, k] / probabilities[i, num_dists, k]
                for j in range(num_dists):
                    gammas[i, j, k] = temp * probabilities[i, j, k]
    return None

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def update_parameter_estimates(
        np.ndarray[DTYPE_64_t, ndim=1] observations,
        np.ndarray[DTYPE_64_t, ndim=1] pi,
        np.ndarray[DTYPE_64_t, ndim=2] transitions,
        np.ndarray[DTYPE_64_t, ndim=3] distributions,
        np.ndarray[DTYPE_64_t, ndim=1] new_pi,
        np.ndarray[DTYPE_64_t, ndim=2] new_transitions,
        np.ndarray[DTYPE_64_t, ndim=3] new_distributions,
        np.ndarray[DTYPE_64_t, ndim=3] etas,
        np.ndarray[DTYPE_64_t, ndim=3] gammas,
        int finalize):
    cdef int i, j, k
    cdef double total, temp, change
    cdef int num_states = gammas.shape[0]
    cdef int num_dists = gammas.shape[1]
    cdef int num_obs = gammas.shape[2]
    with nogil:
        change = 0.0
        for i in range(num_states):
            for j in range(num_states):
                new_pi[i] += etas[i, j, 0]
        for k in range(num_obs - 1):
            for i in range(num_states):
                for j in range(num_states):
                    new_transitions[i, j] += etas[i, j, k]
        for k in range(num_obs):
            for i in range(num_states):
                for j in range(num_dists):
                    new_distributions[i, j, 0] += gammas[i, j, k] 
                    new_distributions[i, j, 1] += gammas[i, j, k] * observations[k]
                    temp = observations[k] - distributions[i, j, 1]
                    new_distributions[i, j, 2] += gammas[i, j, k] * temp * temp
        if finalize == 1:
            total = 0.0
            for i in range(num_states):
                total += new_pi[i]
            for i in range(num_states):
                new_pi[i] /= total
                temp = new_pi[i] - pi[i]
                if temp > 0.0:
                    change = max(change, temp)
                else:
                    change = max(change, -temp)
                pi[i] = new_pi[i]
            for i in range(num_states):
                total = 0.0
                for j in range(num_states):
                    total += new_transitions[i, j]
                for j in range(num_states):
                    new_transitions[i, j] /= total
                    temp = new_transitions[i, j] - transitions[i, j]
                    if temp > 0.0:
                        change = max(change, temp)
                    else:
                        change = max(change, -temp)
                    transitions[i, j] = new_transitions[i, j]
            for i in range(num_states):
                total = 0.0
                for j in range(num_dists):
                    total += new_distributions[i, j, 0]
                    if new_distributions[i, j, 0] > 0:
                        new_distributions[i, j, 1] /= new_distributions[i, j, 0]
                        new_distributions[i, j, 2] /= new_distributions[i, j, 0]
                for j in range(num_dists):
                    if new_distributions[i, j, 0] > 0:
                        new_distributions[i, j, 0] /= total
                        temp = new_distributions[i, j, 0] - distributions[i, j, 0]
                        if temp > 0.0:
                            change = max(change, temp)
                        else:
                            change = max(change, -temp)
                        distributions[i, j, 0] = new_distributions[i, j, 0]
                        temp = new_distributions[i, j, 1] - distributions[i, j, 1]
                        if temp > 0.0:
                            change = max(change, temp)
                        else:
                            change = max(change, -temp)
                        distributions[i, j, 1] = new_distributions[i, j, 1]
                        temp = new_distributions[i, j, 2] - distributions[i, j, 2]
                        if temp > 0.0:
                            change = max(change, temp)
                        else:
                            change = max(change, -temp)
                        distributions[i, j, 2] = new_distributions[i, j, 2]
    return change

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def find_path(
        np.ndarray[DTYPE_64_t, ndim=1] observations,
        np.ndarray[DTYPE_64_t, ndim=1] pi,
        np.ndarray[DTYPE_64_t, ndim=2] transitions,
        np.ndarray[DTYPE_64_t, ndim=2] probs,
        np.ndarray[DTYPE_64_t, ndim=2] scores,
        np.ndarray[DTYPE_int_t, ndim=2] paths,
        np.ndarray[DTYPE_int_t, ndim=1] states):
    cdef int i, j, k
    cdef double prob
    cdef int num_states = probs.shape[0]
    cdef int num_obs = probs.shape[1]
    with nogil:
        for i in range(num_states):
            scores[i, 0] = probs[i, 0] + pi[i]
        for k in range(1, num_obs):
            for i in range(num_states):
                scores[i, k] = -Inf
                for j in range(num_states):
                    prob = scores[j, k - 1] + transitions[j, i]
                    if prob > scores[i, k]:
                        scores[i, k] = prob
                        paths[i, k - 1] = j
                scores[i, k] += probs[i, k]
        prob = scores[0, num_obs - 1]
        states[num_obs - 1] = 0
        for i in range(1, num_states):
            if scores[i, num_obs - 1] > prob:
                prob = scores[i, num_obs - 1]
                states[num_obs - 1] = i
        for i in range(2, num_obs + 1):
            j = num_obs - i
            states[j] = paths[states[j + 1], j]
    return prob

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def generate_sequence(
        np.ndarray[DTYPE_64_t, ndim=2] rand_nums,
        np.ndarray[DTYPE_64_t, ndim=1] pi,
        np.ndarray[DTYPE_64_t, ndim=2] transitions,
        np.ndarray[DTYPE_64_t, ndim=3] distributions,
        np.ndarray[DTYPE_64_t, ndim=1] pi_sums,
        np.ndarray[DTYPE_64_t, ndim=2] transition_sums,
        np.ndarray[DTYPE_64_t, ndim=2] weight_sums,
        np.ndarray[DTYPE_int_t, ndim=1] states,
        np.ndarray[DTYPE_int_t, ndim=1] distribution):
    cdef int i, j, state, dist
    cdef int num_obs = states.shape[0]
    cdef int num_states = pi.shape[0]
    cdef int num_dists = distributions.shape[1]
    with nogil:
        pi_sums[0] = pi[0]
        for i in range(1, num_states):
            pi_sums[i] = pi_sums[i - 1] + pi[i]
        for i in range(num_states):
            transition_sums[i, 0] = transitions[i, 0]
            for j in range(1, num_states):
                transition_sums[i, j] = transition_sums[i, j - 1] + transitions[i, j]
            weight_sums[i, 0] = distributions[i, 0, 0]
            for j in range(1, num_dists):
                weight_sums[i, j] = weight_sums[i, j - 1] + distributions[i, j, 0]
        state = 0
        while rand_nums[0, 0] > pi_sums[state]:
            state += 1
        states[0] = state
        dist = 0
        while rand_nums[0, 1] > weight_sums[state, dist]:
            dist += 1
        distribution[0] = dist
        for i in range(1, num_obs):
            state = 0
            while rand_nums[i, 0] > transition_sums[states[i - 1], state]:
                state += 1
            states[i] = state
            dist = 0
            while rand_nums[i, 1] > weight_sums[state, dist]:
                dist += 1
            distribution[i] = dist
    return None




