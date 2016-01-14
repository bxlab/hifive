#!/usr/bin/env python

import sys

import numpy
import _hmm

class HMM:
    """A simple HMM implementation with dependency only on numpy able to handle discrete and gaussian-mixture HMMs."""

    def __init__(self, num_states=2, num_distributions=1, transitions=None, pi=None, distributions=None, seed=None):
        """Initialize the model with values some or all of the arguments and randomly initializing the remainder."""
        self.num_states = num_states
        self.num_dists = num_distributions
        # There is no error checking to ensure that the number of states/distributions match the shape of passed values
        if seed is None:
            self.rng = numpy.random.RandomState()
        else:
            self.rng = numpy.random.RandomState(seed=int(seed))
        if transitions is not None:
            self.transitions = numpy.array(transitions, dtype=numpy.float64)
        else:
            self.transitions = self.rng.rand(num_states, num_states).astype(numpy.float64)
        self.transitions /= numpy.sum(self.transitions, axis=1).reshape(-1, 1)
        if pi is not None:
            self.pi = numpy.array(pi, dtype=numpy.float64)
        else:
            self.pi = self.rng.rand(num_states).astype(numpy.float64)
        self.pi /= numpy.sum(self.pi, axis=0)
        # distributions should be passed as a nested list or numpy array with the following characteristics:
        # shape - [num_states, num_mixtures, 3]
        # The last dimension contains the guassian components weight, mean, and variance
        if distributions is None:
            self.distributions = self.rng.rand(self.num_states, self.num_dists, 3).astype(numpy.float64)
        else:
            self.distributions = numpy.array(distributions).astype(numpy.float64)
        self.distributions[:, :, 0] /= numpy.sum(self.distributions[:, :, 0], axis=1).reshape(-1, 1)

    def find_probabilities(self, observations):
        """Calculate the probability of observations for each state and, if applicable, for each mixture component."""
        self.probs = numpy.zeros((self.num_states, self.num_dists + 1, observations.shape[0]), dtype=numpy.float64)
        _hmm.find_probabilities(observations, self.distributions, self.probs)

    def find_alphas(self):
        """Calculate the probability of being in a given state at a given all previous observations."""
        # self.alphas [num_states, num_observations] is the probability of being in that state after observations up
        # to that point. For example self.alphas[1, 10] is the probability of being in state1 after the first 10
        # observations
        # self.scalars [num_observations] rescales each component so there is no variable overflow.
        self.alphas = numpy.zeros((self.num_states, self.probs.shape[2]), dtype=numpy.float64)
        self.scalars = numpy.zeros(self.probs.shape[2], dtype=numpy.float64)
        _hmm.find_alphas(self.probs, self.pi, self.transitions, self.alphas, self.scalars)

    def find_betas(self):
        """Calculate the probability of reaching a state given all of the future observations."""
        # self.betas [num_states, num_observations] is the probability of being in that state after observations past
        # that point. For example self.betas[1, 10] is the probability of being in state1 followed by observations 11
        # and higher.
        self.betas = numpy.zeros((self.num_states, self.probs.shape[2]), dtype=numpy.float64)
        _hmm.find_betas(self.probs, self.transitions, self.scalars, self.betas)

    def find_etas(self):
        """Calculate the probability of being in state i and time t-1 and transitioning to state j at time t."""
        # self.etas [num_states, num_states, num_observations - 1] is the probability of observing the transition at
        # time t-1 to time t from state i to state j.
        self.etas = numpy.zeros((self.num_states, self.num_states, self.probs.shape[2] - 1), dtype=numpy.float64)
        _hmm.find_etas(self.probs, self.transitions, self.alphas, self.betas, self.etas)

    def find_gammas(self):
        """Calculate the probability of on observation coming from state i at time t and mixture component m (if continuous)."""
        # self.gamma [num_states, num_mixtures, num_observations] is the probability of an observation arising
        # from a particular component in a particular state.
        self.gammas = numpy.zeros((self.num_states, self.num_dists, self.probs.shape[2]), dtype=numpy.float64)
        _hmm.find_gammas(self.probs, self.alphas, self.betas, self.gammas)

    def train(self, observations, threshold=1e-6, max_iterations=1000):
        """Optimize the model parameters based on a set of training sequences."""
        # Training is stopped when either the maximum number of iterations is reached or the largest change in a
        # parameter is smaller than the threshold.
        change = numpy.inf
        iteration = 0
        new_pi = numpy.copy(self.pi)
        new_transitions = numpy.copy(self.transitions)
        new_distributions = numpy.copy(self.distributions)
        while change > threshold and iteration < max_iterations:
            new_pi.fill(0.0)
            new_transitions.fill(0.0)
            new_distributions.fill(0.0)
            for h in range(len(observations)):
                self.find_probabilities(observations[h])
                self.find_alphas()
                self.find_betas()
                self.find_etas()
                self.find_gammas()
                if h == len(observations) - 1:
                    finalize = 1
                else:
                    finalize = 0
                change = _hmm.update_parameter_estimates(observations[h], self.pi, self.transitions,
                                                         self.distributions, new_pi, new_transitions,
                                                         new_distributions, self.etas, self.gammas, finalize)
            print >> sys.stderr, ("\rIteration: %03i   Change: %.5f") % (iteration, change),
            iteration += 1
        print >> sys.stderr, ("\n"),
        return

    def find_path(self, observations):
        self.find_probabilities(observations)
        scores = numpy.zeros((self.num_states, observations.shape[0]), dtype=numpy.float64)
        paths = numpy.zeros((self.num_states, observations.shape[0]), dtype=numpy.int32)
        states = numpy.zeros(observations.shape[0], dtype=numpy.int32)
        log_pi = numpy.zeros(self.pi.shape, dtype=numpy.float64) - numpy.inf
        where = numpy.where(self.pi > 0.0)
        log_pi[where] = numpy.log(self.pi[where])
        log_transitions = numpy.zeros(self.transitions.shape, dtype=numpy.float64) - numpy.inf
        where = numpy.where(self.transitions > 0.0)
        log_transitions[where] = numpy.log(self.transitions[where])
        log_probs = numpy.zeros((self.probs.shape[0], self.probs.shape[2]), dtype=numpy.float64) - numpy.inf
        where = numpy.where(self.probs[:, -1, :] > 0.0)
        log_probs[where] = numpy.log(self.probs[where[0], -1, where[1]])
        ll = _hmm.find_path(observations, log_pi, log_transitions, log_probs, scores, paths, states)
        return states, ll

    def generate_sequence(self, seq_len=100):
        """Generate a set of hidden states and emmissions from the model."""
        states = numpy.zeros(seq_len, dtype=numpy.int32)
        distributions = numpy.zeros(seq_len, dtype=numpy.int32)
        observations = numpy.zeros(seq_len, dtype=numpy.float64)
        pi_sums = numpy.zeros(self.pi.shape, dtype=numpy.float64)
        transition_sums = numpy.zeros(self.transitions.shape, dtype=numpy.float64)
        weight_sums = numpy.zeros((self.num_states, self.num_dists), dtype=numpy.float64)
        rand_nums = self.rng.rand(seq_len, 2)
        _hmm.generate_sequence(rand_nums, self.pi, self.transitions, self.distributions, pi_sums,
                               transition_sums, weight_sums, states, distributions)
        observations[:] = self.rng.normal(self.distributions[states, distributions, 1],
                                          self.distributions[states, distributions, 2] ** 0.5)
        return states, observations

def test():
    hmm1 = HMM(
        seed=2001,
        num_distributions=2,
        num_states=2,
        pi=[0.5, 0.5],
        transitions=[[0.95, 0.05], [0.1, 0.9]],
        distributions=[[[0.5, 1.5, 0.25], [0.5, -1.0, 0.25]], [[0.25, 2.5, 0.25], [0.75, 4.5, 0.5]]],
        )
    observations = []
    for i in range(100):
        temp = hmm1.generate_sequence(100)
        observations.append(temp[1])
    hmm2 = HMM(
        seed=2002,
        num_distributions=2,
        num_states=2,
        )
    hmm2.train(observations, threshold=1e-6, max_iterations=400)
    print 'True Pi: %s    Learned Pi: %s' % (str(list(hmm1.pi)), str(list(hmm2.pi)))
    for i in range(hmm1.num_states):
        print 'State %i -' % i
        print 'True transitions: %s    Learned transitions: %s' % (str(list(hmm1.transitions[i, :])),
                                                                   str(list(hmm2.transitions[i, :])))
        print 'True distribution weights: %s    Learned distribution weights %s' % (
            str(list(hmm1.distributions[i, :, 0])),
            str(list(hmm2.distributions[i, :, 0])))
        print 'True distribution means: %s    Learned distribution means %s' % (
            str(list(hmm1.distributions[i, :, 1])),
            str(list(hmm2.distributions[i, :, 1])))
        print 'True distribution variances: %s    Learned distribution variances %s' % (
            str(list(hmm1.distributions[i, :, 2])),
            str(list(hmm2.distributions[i, :, 2])))
