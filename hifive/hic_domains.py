#!/usr/bin/env python

"""
This is a module contains scripts for finding domains in HiC interaction data.

Concepts
--------

These functions rely on the :class:`HiC` class in conjunction with the :class:`Fend` and :class:`HiCData` classes.


API Documentation
-----------------
"""

import sys

import numpy
import scipy
import h5py
try:
    from mpi4py import MPI
except:
    pass
try:
    from pyx import *
except:
    pass

import libraries._hic_domains as _hic_domains
from libraries.hmm import HMM
import hic_binning
import plotting


class TAD( object ):
    """
    """

    def __init__(self, hic, silent=False):
        self.hic = hic
        self.silent = silent

    def __getitem__(self, key):
        """Dictionary-like lookup."""
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        """Dictionary-like value setting."""
        self.__dict__[key] = value
        return None

    """
    def find_BI_TADs(self, binsize, width, minbins, maxbins, chroms=[]):
        if isinstance(chroms, str):
            chroms = [chroms]
        if len(chroms) == 0:
            chroms = list(self.hic.fends['chromosomes'][...])
        self.binsize = int(binsize)
        self.minsize = self.binsize * max(1, int(minbins))
        self.maxsize = self.binsize * (int(maxbins) + 1)
        self.width = int(width) * self.binsize
        maxsize = self.binsize * (int(maxbins) + int(width) + 1)
        self.TADs = {}
        for chrom in chroms:
            self.TADs[chrom] = []
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding heatmap for chromosome %s...") % (' ' * 80, chrom),
            temp = self.hic.cis_heatmap(chrom, binsize=binsize, datatype='enrichment', arraytype='compact',
                                        maxdistance=maxsize, returnmapping=True)
            if temp is None:
                continue
            data, mapping = temp
            print >> sys.stderr, ("\rFinding BI scores for chromosome %s...") % (chrom),
            BIs = numpy.zeros((data.shape[0], maxbins + 1, 2), dtype=numpy.float32)
            BIs.fill(-numpy.inf)
            _hic_domains.find_BIs(data, BIs, minbins, maxbins, width)
            path = numpy.zeros(data.shape[0], dtype=numpy.int32)
            path_scores = numpy.zeros(data.shape[0], dtype=numpy.float64)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding optimal domains for chromosome %s...") % (' ' * 80, chrom),
            _hic_domains.find_BI_path(BIs, path, path_scores, minbins, maxbins)
            i = path.shape[0] - 1
            domains = []
            while i > 0:
                if path[i] != 1:
                    domains.append([i - path[i] + 1, i])
                    self.TADs[chrom].append([mapping[domains[-1][0], 0], mapping[domains[-1][1], 1]])
                i -= path[i]
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinished finding TADs\n") % (' ' * 80),

    def find_arrowhead_TADs(self, binsize, minbins, maxbins, chroms=[]):
        if isinstance(chroms, str):
            chroms = [chroms]
        if len(chroms) == 0:
            chroms = list(self.hic.fends['chromosomes'][...])
        self.binsize = int(binsize)
        self.minsize = self.binsize * max(1, int(minbins))
        self.maxsize = self.binsize * (int(maxbins) + 1)
        self.TADs = {}
        for chrom in chroms:
            self.TADs[chrom] = []
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding heatmap for chromosome %s...") % (' ' * 80, chrom),
            temp = self.hic.cis_heatmap(chrom, binsize=binsize * 16, datatype='fend', arraytype='full', returnmapping=True)
            if temp is None:
                continue
            temp_data, mapping = temp
            heatmap = numpy.zeros((temp_data.shape[0] * 16, temp_data.shape[0] * 16, 2), dtype=numpy.float32)
            temp = numpy.zeros(heatmap.shape, dtype=numpy.float32)
            for i in range(16):
                for j in range(16):
                    heatmap[i::16, j::16, :] += temp_data
            temp_data = self.hic.cis_heatmap(chrom, binsize=binsize * 4, datatype='fend', arraytype='full', start=mapping[0, 0], stop=mapping[-1, 1])
            for i in range(4):
                for j in range(4):
                    temp[i::4, j::4, :] += temp_data
            where = numpy.where(temp[:, :, 0] > 0)
            heatmap[where[0], where[1], :] = temp[where[0], where[1], :]
            temp_data, mapping = self.hic.cis_heatmap(chrom, binsize=binsize, datatype='fend', arraytype='full', start=mapping[0, 0], stop=mapping[-1, 1], returnmapping=True)
            where = numpy.where(temp_data[:, :, 0] > 0)
            heatmap[where[0], where[1], :] = temp_data[where[0], where[1], :]
            data = numpy.zeros((heatmap.shape[0], maxbins - 1, 2), dtype=numpy.float32)
            for i in range(heatmap.shape[0] - 1):
                data[i, :min(data.shape[1], data.shape[0] - i - 1), :] = heatmap[i, (i + 1):min(data.shape[1] + i + 1, data.shape[0]), :]
            where = numpy.where(data[:, :, 1] > 0)
            data[where[0], where[1], 0] /= data[where[0], where[1], 1]
            scores = numpy.zeros((data.shape[0], data.shape[1]), dtype=numpy.float32)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding arrowhead transformation for chromosome %s...") % (' ' * 80, chrom),
            _hic_domains.find_arrowhead_transformation(data, scores, maxbins)
            sums = numpy.zeros(data.shape, dtype=numpy.float32)
            signs = numpy.zeros(data.shape, dtype=numpy.float32)
            variances = numpy.zeros((data.shape[0], data.shape[1], 2, 2), dtype=numpy.float32)
            domain_scores = numpy.zeros(scores.shape, dtype=numpy.float32)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding arrowhead scoring for chromosome %s...") % (' ' * 80, chrom),
            _hic_domains.find_arrowhead_scores(scores, sums, signs, variances, domain_scores, minbins)
            path = numpy.zeros(data.shape[0], dtype=numpy.int32)
            path_scores = numpy.zeros(data.shape[0], dtype=numpy.float64)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding optimal domains for chromosome %s...") % (' ' * 80, chrom),
            _hic_domains.find_arrowhead_path(domain_scores, path, path_scores, minbins, maxbins)
            i = path.shape[0] - 1
            domains = []
            while i > 0:
                if path[i] != 1:
                    domains.append([i - path[i], i])
                    self.TADs[chrom].append([mapping[domains[-1][0], 0], mapping[domains[-1][1], 1]])
                i -= path[i]
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinished finding TADs\n") % (' ' * 80),
    """

    def find_DI_TADs(self, binsize=20000, step=2500, window=500000, minsize=25000, maxsize=1500000, smoothing=6,
                     joindomains=True, chroms=[]):
        self.binsize = int(binsize)
        self.step = int(step)
        self.window = int(window)
        self.smoothing = int(smoothing)
        steps = binsize / step
        self.DIs = {}
        self.TADs = {}
        if isinstance(chroms, str):
            chroms = [chroms]
        if len(chroms) == 0:
            chroms = list(self.hic.fends['chromosomes'][...])
        pos_sum = 0.0
        pos_2sum = 0.0
        pos_count = 0
        neg_sum = 0.0
        neg_2sum = 0.0
        neg_count = 0
        training_seqs = []
        for chrom in chroms:
            temp = self.hic.cis_heatmap(chrom, binsize=self.step, maxdistance=self.window, datatype='fend',
                                        arraytype='compact', returnmapping=True)
            if temp is None:
                continue
            num_bins = temp[0].shape[1]
            scores = []
            positions = []
            for i in range(num_bins, temp[0].shape[0] - num_bins - steps + 1):
                observed0 = 0
                expected0 = 0.0
                observed1 = 0
                expected1 = 0.0
                for j in range(steps):
                    temp1 = numpy.arange(steps - j - 1, num_bins - j)
                    observed0 += numpy.sum(temp[0][i + j, temp1, 0])
                    expected0 += numpy.sum(temp[0][i + j, temp1, 1])
                    temp1 = numpy.arange(i + steps - num_bins - 1, i)
                    temp2 = i + j - temp1 - 1
                    observed1 += numpy.sum(temp[0][temp1, temp2, 0])
                    expected1 += numpy.sum(temp[0][temp1, temp2, 1])
                if observed0 > 0 and observed1 > 0:
                    scores.append(numpy.log(observed0 * expected1 / (observed1 * expected0)))
                    positions.append((temp[1][i, 0] + temp[1][i + steps - 1, 1]) / 2)
            self.DIs[chrom] = numpy.zeros(len(scores), dtype=numpy.dtype([('position', numpy.int32),
                                                                          ('score', numpy.float64)]))
            scores = numpy.array(scores)
            if self.smoothing > 1:
                scores *= self.smoothing
                for j in range(1, self.smoothing):
                    scores[j:] += scores[:-j] * (self.smoothing - j)
                    scores[:-j] += scores[j:] * (self.smoothing - j)
                scores[self.smoothing:-self.smoothing] /= self.smoothing ** 2
                for j in range(1, self.smoothing + 1):
                    scores[self.smoothing - j] /= self.smoothing ** 2 - (j * (j + 1)) / 2
                    scores[-self.smoothing + j - 1] /= self.smoothing ** 2 - (j * (j + 1)) / 2
            scores /= numpy.std(scores)
            self.DIs[chrom]['score'] = scores
            self.DIs[chrom]['position'] = positions
            training_seqs.append(numpy.array(scores, dtype=numpy.float64))
            where = numpy.where(self.DIs[chrom]['score'] > 0.0)[0]
            pos_sum += numpy.sum(self.DIs[chrom]['score'][where])
            pos_2sum += numpy.sum(self.DIs[chrom]['score'][where] ** 2)
            pos_count += where.shape[0]
            where = numpy.where(self.DIs[chrom]['score'] < 0.0)[0]
            neg_sum += numpy.sum(self.DIs[chrom]['score'][where])
            neg_2sum += numpy.sum(self.DIs[chrom]['score'][where] ** 2)
            neg_count += where.shape[0]
        pos_mean = pos_sum / pos_count
        pos_var = (pos_2sum / pos_count - pos_mean ** 2.0)
        neg_mean = neg_sum / neg_count
        neg_var = (neg_2sum / neg_count - neg_mean ** 2.0)
        seed = 2001
        pi = [0.5, 0.5]
        transitions = [[0.98, 0.02],
                       [0.02, 0.98]]
        distributions = [[[0.33, pos_mean * 1.25, pos_var / 3.0],
                          [0.33, pos_mean * 0.75, pos_var / 3.0],
                          [0.33, pos_mean * 0.25, pos_var / 3.0]],
                          [[0.33, neg_mean * 1.25, neg_var / 3.0],
                          [0.33, neg_mean * 0.75, neg_var / 3.0],
                          [0.33, neg_mean * 0.25, neg_var /3.0]]]
        hmm = HMM(
            seed=seed,
            num_states=2,
            num_distributions=2,
            pi=pi,
            transitions=transitions,
            distributions=distributions)
        hmm.train(training_seqs)
        for chrom in chroms:
            if chrom not in self.DIs:
                continue
            tads = []
            states = hmm.find_path(self.DIs[chrom]['score'])[0]
            i = 1
            while i < states.shape[0]:
                while i < states.shape[0] and not (states[i - 1] != 0 and states[i] == 0):
                    i += 1
                if i < states.shape[0]:
                    start = self.DIs[chrom]['position'][i] - self.step / 2
                    i += 1
                while i < states.shape[0] and not (states[i] != 1 and states[i - 1] == 1):
                    i += 1
                if i < states.shape[0]:
                    stop = self.DIs[chrom]['position'][i - 1] + self.step / 2
                    if stop - start >= minsize:
                        tads.append([start, stop])
                    start = None
            self.TADs[chrom] = numpy.array(tads, dtype=numpy.int32)
            where = numpy.where(self.TADs[chrom][1:, 0] < self.TADs[chrom][:-1, 1])[0]
            mids = (self.TADs[chrom][where, 1] + self.TADs[chrom][where + 1, 0]) / 2
            self.TADs[chrom][where, 1] = mids
            self.TADs[chrom][where + 1, 0] = mids
            if joindomains:
                temp = self.hic.cis_heatmap(chrom, binsize=self.step, maxdistance=maxsize, datatype='enrichment',
                                            arraytype='compact', returnmapping=True)
                num_tads = self.TADs[chrom].shape[0]
                new_num_tads = 0
                current_tads = numpy.zeros((num_tads, 3), dtype=numpy.float32)
                starts = numpy.searchsorted(temp[1][:, 0], self.TADs[chrom][:, 0])
                stops = numpy.searchsorted(temp[1][:, 0], self.TADs[chrom][:, 1], side='right')
                for i in range(num_tads):
                    for j in range(starts[i], stops[i] - 1):
                        current_tads[i, :2] += numpy.sum(temp[0][j, :(stops[i] - j - 1), :], axis=0)
                current_tads[:, 2] = numpy.log(current_tads[:, 0] / current_tads[:, 1])
                while new_num_tads != num_tads:
                    print num_tads
                    num_tads = self.TADs[chrom].shape[0]
                    join_scores = numpy.zeros((num_tads - 1, 3), dtype=numpy.float32)
                    for i in range(num_tads - 1):
                        if (stops[i] - starts[i]) * self.step > maxsize:
                            continue
                        for j in range(starts[i], stops[i]):
                            join_scores[i, :2] += numpy.sum(temp[0][j, (stops[i] - j - 1):(stops[i + 1] - j - 1), :], axis=0)
                    where = numpy.where(join_scores[:, 0] > 0)[0]
                    join_scores[where, 2] = numpy.log(join_scores[where, 0] / join_scores[where, 1])
                    new_tads = []
                    new_starts = []
                    new_stops = []
                    new_scores = []
                    i = 0
                    while i < num_tads:
                        if ((i < num_tads - 1 and (stops[i] - starts[i]) * self.step <= maxsize) and
                            ((stops[i] - i == 0 and join_scores[i, 2] >= max((current_tads[i, 2] +
                                current_tads[i + 1, 2]) / 2., join_scores[i + 1, 2])) or
                            (i > 0 and i < num_tads - 2 and join_scores[i, 2] >=
                                max((current_tads[i, 2] + current_tads[i + 1, 2]) / 2., join_scores[i + 1, 2],
                                join_scores[i - 1, 2])) or
                            (i == num_tads - 2 and join_scores[i, 2] >= max((current_tads[i, 2] + current_tads[i + 1, 2]) / 2.,
                                join_scores[i - 1, 2])))):
                            new_tads.append([self.TADs[chrom][i, 0], self.TADs[chrom][i + 1, 1]])
                            new_starts.append(starts[i])
                            new_stops.append(stops[i + 1])
                            new_scores.append([current_tads[i, 0] + current_tads[i + 1, 0] + join_scores[i, 0],
                                                current_tads[i, 1] + current_tads[i + 1, 1] + join_scores[i, 1], 0])
                            new_scores[-1][-1] = numpy.log(new_scores[-1][0] / new_scores[-1][1])
                            i += 2
                        else :
                            new_tads.append([self.TADs[chrom][i, 0], self.TADs[chrom][i, 1]])
                            new_starts.append(starts[i])
                            new_stops.append(stops[i])
                            new_scores.append(current_tads[i, :])
                            i += 1
                        #for i in range(num_tads - 1):
                        #    print current_tads[i, 2], join_scores[i,2]
                        #print current_tads[-1,2]
                    new_num_tads = len(new_tads)
                    self.TADs[chrom] = numpy.array(new_tads, dtype=numpy.int32)
                    current_tads = numpy.array(new_scores, dtype=numpy.float32)
                    starts = numpy.array(new_starts, dtype=numpy.int32)
                    stops = numpy.array(new_stops, dtype=numpy.int32)


    def plot_DI_tads(self, out_fname):
        if 'pyx' not in sys.modules:
            return None
        unit.set(defaultunit="cm")
        text.set(mode="latex")
        text.preamble(r"\usepackage{times}")
        text.preamble(r"\usepackage{sansmath}")
        text.preamble(r"\sansmath")
        text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
        painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-3)] )
        pages = []
        for chrom in self.DIs:
            start = self.DIs[chrom]['position'][0]
            stop = self.DIs[chrom]['position'][-1]
            max_val = numpy.amax(numpy.abs(self.DIs[chrom]['score'])) * 1.05
            g = graph.graphxy( width=40, height=5,
                x=graph.axis.lin(min=start, max=stop, title='', painter=painter), 
                y=graph.axis.lin(min=-max_val, max=max_val, title='', painter=painter), 
                x2=graph.axis.lin(min=0, max=1, parter=None),
                y2=graph.axis.lin(min=0, max=1, parter=None))
            g.text(40.1, 2.5, chrom, [text.halign.center, text.valign.bottom, trafo.rotate(-90), text.size(-3)])
            g.plot(graph.data.points(zip(self.DIs[chrom]['position'], self.DIs[chrom]['score']), x=1, y=2),
                [graph.style.line([style.linewidth.THIN])])
            g.stroke(path.line(0, 2.5, 40, 2.5), [style.linestyle.dotted, style.linewidth.THIN])
            for i in range(self.TADs[chrom].shape[0]):
                X0 = (self.TADs[chrom][i, 0] - start) / float(stop - start) * 40.0
                X1 = (self.TADs[chrom][i, 1] - start) / float(stop - start) * 40.0
                if i % 2 == 0:
                    Y = 1.25
                else:
                    Y = 3.75
                g.stroke(path.line(X0, Y, X1, Y), [style.linewidth.THICK])
                g.stroke(path.line(X0, 1.25, X0, 3.75), [style.linewidth.THIN, style.linestyle.dotted])
                if i == self.TADs[chrom].shape[0] - 1 or self.TADs[chrom][i, 1] != self.TADs[chrom][i + 1, 0]:
                    g.stroke(path.line(X1, 1.25, X1, 3.75), [style.linewidth.THIN, style.linestyle.dotted])
            pages.append(document.page(g))
        doc = document.document(pages)
        doc.writePDFfile(out_fname)    

    def write_TADs(self, fname):
        output = open(fname, 'w')
        chroms = self.TADs.keys()
        chroms.sort()
        for chrom in chroms:
            self.TADs[chrom].sort()
            for domain in self.TADs[chrom]:
                print >> output, "%s\t%i\t%i" % (chrom, domain[0], domain[1])
        output.close()


class Compartment( object ):
    """
    """
    def __init__(self, hic, binsize, chroms=[], minreads=3, out_fname=None, silent=False):
        self.hic = hic
        self.binsize = binsize
        self.minreads = minreads
        if 'mpi4py' in sys.modules.keys():
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.num_procs = self.comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.num_procs = 1
        if self.rank == 0:
            self.silent = silent
        else:
            self.silent = True
        if self.rank == 0:
            if out_fname is not None:
                storage = h5py.File(out_fname, 'a')
            else:
                storage = None
            if chroms == "" or (isinstance(chroms, list) and len(chroms) == 0):
                chroms = hic.fends['chromosomes'][...]
            elif isinstance(chroms, str):
                chroms = [chroms]
            self.chroms = chroms
            for i in range(1, self.num_procs):
                self.comm.send(chroms, dest=i, tag=11)
        else:
            self.chroms = self.comm.recv(source=0, tag=11)
        if hic.binned is None:
            chr_indices = hic.fends['chr_indices'][...]
            fends = hic.fends['fends'][...]
        else:
            chr_indices = hic.fends['bin_indices'][...]
            fends = hic.fends['bins'][...]
        if 'eigenv' not in self.__dict__.keys():
            self.eigenv = {}
            self.positions = {}
        for chrom in self.chroms:
            corrs = None
            chrint = hic.chr2int[chrom]
            startfend = chr_indices[chrint]
            while hic.filter[startfend] == 0:
                startfend += 1
            start = (fends['mid'][startfend] / binsize) * binsize
            stopfend = chr_indices[chrint + 1]
            while hic.filter[stopfend - 1] == 0:
                stopfend -= 1
            stop = ((fends['mid'][stopfend - 1] - 1) / binsize + 1) * binsize
            N = (stop - start) / binsize
            if self.rank == 0:
                if (storage is None or ('%s.eigenv' % chrom not in storage and
                    '%s.correlations' % chrom not in storage)):
                    indices = list(numpy.triu_indices(N, 1))
                    if storage is None or '%s.dbinned_expected' % chrom not in storage:
                        if storage is None or '%s.expected' not in storage:
                            expected = self.hic.cis_heatmap(chrom, binsize=binsize, start=start, stop=stop,
                                datatype='expected', arraytype='full', returnmapping=False, silent=True)[:, :, 1]
                            valid = (numpy.sum(expected > 0, axis=1) > 0).astype(numpy.int32)
                            positions = numpy.zeros((N, 2), dtype=numpy.int32)
                            positions[:, 0] = start + numpy.arange(N) * binsize
                            positions[:, 1] = positions[:, 0] + binsize
                            if storage is not None:
                                storage.create_dataset(name='%s.expected' % chrom, data=expected[indices])
                                storage.create_dataset(name='%s.positions' % chrom, data=positions)
                                storage.create_dataset(name='%s.valid' % chrom, data=valid)
                        else:
                            expected = numpy.zeros((N, N), dtype=numpy.float32)
                            expected[indices] = storage['%s.expected' % chrom][...]
                            expected[indices[1], indices[0]] = expected[indices]
                            valid = storage['%s.valid' % chrom][...]
                            positions = storage['%s.positions' % chrom][...]
                        if storage is None or '%s.counts' not in storage:
                            start_index = hic.data['cis_indices'][startfend]
                            stop_index = hic.data['cis_indices'][stopfend]
                            data = hic.data['cis_data'][start_index:stop_index, :].astype(numpy.int64)
                            data = data[numpy.where(hic.filter[data[:, 0]] * hic.filter[data[:, 1]])[0], :]
                            data[:, :2] = fends['mid'][data[:, :2]]
                            data[:, :2] -= start
                            data[:, :2] /= binsize
                            counts = numpy.bincount(data[:, 0] * N + data[:, 1], weights=data[:, 2],
                                                    minlength=(N * N)).reshape(N, N).astype(numpy.int32)
                            counts[indices[1], indices[0]] = counts[indices]
                            if storage is not None:
                                storage.create_dataset(name='%s.counts' % chrom, data=counts[indices])
                        else:
                            counts = numpy.zeros((N, N), dtype=numpy.int32)
                            counts[indices] = storage['%s.counts' % chrom][...]
                            counts[indices[1], indices[0]] = counts[indices]
                        for i in range(1, self.num_procs):
                            self.comm.send(0, dest=i, tag=11)
                        binned_c, binned_e = self._dynamically_bin(counts, expected, valid)
                        if storage is not None:
                            storage.create_dataset(name='%s.dbinned_expected' % chrom, data=binned_e[indices])
                            storage.create_dataset(name='%s.dbinned_counts' % chrom, data=binned_c[indices])
                    else:
                        binned_c = numpy.zeros((N, N), dtype=numpy.int32)
                        binned_c[indices] = storage['%s.dbinned_counts' % chrom][...]
                        binned_c[indices[1], indices[0]] = binned_c[indices]
                        binned_e = numpy.zeros((N, N), dtype=numpy.float32)
                        binned_e[indices] = storage['%s.dbinned_expected' % chrom][...]
                        binned_e[indices[1], indices[0]] = binned_e[indices]
                        valid = storage['%s.valid' % chrom][...]
                        positions = storage['%s.positions' % chrom][...]
                        for i in range(1, self.num_procs):
                            self.comm.send(1, dest=i, tag=11)
                    corrs = self._find_correlations(binned_c, binned_e, valid)
                    if storage is not None:
                        storage.create_dataset(name='%s.correlations' % chrom,
                                               data=corrs[numpy.triu_indices(numpy.sum(valid), 1)])
                else:
                    for i in range(1, self.num_procs):
                        self.comm.send(-1, dest=i, tag=11)
                if (corrs is None and storage is not None and
                    '%s.eigenv' % chrom not in storage and
                    '%s.correlations' % chrom in storage):
                    valid = storage['%s.valid' % chrom][...]
                    N = numpy.sum(valid)
                    indices = numpy.triu_indices(N, 1)
                    corrs = numpy.ones((N, N), dtype=numpy.float32)
                    corrs[indices] = storage['%s.correlations' % chrom][...]
                    corrs[indices[1], indices[0]] = corrs[indices]
                    positions = storage['%s.positions' % chrom][...]
                if storage is None or '%s.eigenv' % chrom not in storage:
                    self.eigenv[chrom] = scipy.sparse.linalg.eigs(corrs, k=1)[1][:, 0]
                    if storage is not None:
                        storage.create_dataset(name="%s.eigenv" % chrom, data=self.eigenv[chrom])
                    self.positions[chrom] = positions[numpy.where(valid)[0], :]
                else:
                    self.eigenv[chrom] = storage['%s.eigenv' % chrom][...]
                    valid = storage['%s.valid' % chrom][...]
                    self.positions[chrom] = storage['%s.positions' % chrom][...][numpy.where(valid)[0], :]
                for i in range(1, self.num_procs):
                    self.comm.send(self.eigenv[chrom], dest=i, tag=11)
                    self.comm.send(self.positions[chrom], dest=i, tag=11)
            else:
                job = self.comm.recv(source=0, tag=11)
                if job == 0:
                    self._dynamically_bin()
                if job >= 0:
                    self._find_correlations()
                self.eigenv[chrom] = self.comm.recv(source=0, tag=11)
                self.positions[chrom] = self.comm.recv(source=0, tag=11)
        if self.rank == 0:
            storage.close()

    def _dynamically_bin(self, counts=None, expected=None, valid=None):
        if self.rank == 0:
            if not self.silent:
                print >> sys.stderr, ("\rDynamically binning"),
            N = counts.shape[0]
            node_ranges = numpy.round(numpy.linspace(0, N * (N - 1) / 2, self.num_procs + 1)).astype(numpy.int32)
            indices = numpy.triu_indices(N, 1)
            for i in range(1, self.num_procs):
                self.comm.send(N, dest=i, tag=11)
                self.comm.send(counts, dest=i, tag=13)
                self.comm.send(expected, dest=i, tag=13)
                self.comm.send(valid, dest=i, tag=13)
                self.comm.send(indices[0][node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                self.comm.send(indices[1][node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
            binned_c = numpy.zeros(counts.shape, dtype=numpy.int32)
            binned_e = numpy.zeros(counts.shape, dtype=numpy.float32)
            _hic_domains.dynamically_bin(counts, expected, binned_c, binned_e, valid,
                indices[0][:node_ranges[1]].astype(numpy.int32), indices[1][:node_ranges[1]].astype(numpy.int32),
                self.minreads)
            for i in range(1, self.num_procs):
                binned_c[indices[0][node_ranges[i]:node_ranges[i + 1]],
                         indices[1][node_ranges[i]:node_ranges[i + 1]]] = self.comm.recv(source=i, tag=11)
                binned_e[indices[0][node_ranges[i]:node_ranges[i + 1]],
                         indices[1][node_ranges[i]:node_ranges[i + 1]]] = self.comm.recv(source=i, tag=11)
            binned_c[indices[1], indices[0]] = binned_c[indices]
            binned_e[indices[1], indices[0]] = binned_e[indices]
            if not self.silent:
                print >> sys.stderr, ("\r%s\r") % (' ' * 80),
            return binned_c, binned_e
        else:
            N = self.comm.recv(source=0, tag=11)
            counts = self.comm.recv(source=0, tag=13)
            expected = self.comm.recv(source=0, tag=13)
            valid = self.comm.recv(source=0, tag=13)
            indices0 = self.comm.recv(source=0, tag=11)
            indices1 = self.comm.recv(source=0, tag=11)
            binned_c = numpy.zeros((N, N), dtype=numpy.int32)
            binned_e = numpy.zeros((N, N), dtype=numpy.float32)
            _hic_domains.dynamically_bin(counts, expected, binned_c, binned_e, valid, indices0.astype(numpy.int32),
                            indices1.astype(numpy.int32), self.minreads)
            self.comm.send(binned_c[indices0, indices1], dest=0, tag=11)
            self.comm.send(binned_e[indices0, indices1], dest=0, tag=11)
        return None

    def _find_correlations(self, counts=None, expected=None, valid=None):
        if self.rank == 0:
            if not self.silent:
                print >> sys.stderr, ("\rFinding correlations"),
            valid2 = numpy.where(valid)[0]
            M = valid2.shape[0]
            data = numpy.ones((M, M), dtype=numpy.float32)
            indices = numpy.triu_indices(M, 1)
            data[indices] = counts[valid2, :][:, valid2][indices]
            data[indices] /= expected[valid2, :][:, valid2][indices]
            data[indices[1], indices[0]] = data[indices]
            data = numpy.log(data)
            data -= numpy.mean(data, axis=1).reshape(-1, 1)
            data /= numpy.std(data, axis=1).reshape(-1, 1)
            node_ranges = numpy.round(numpy.linspace(0, indices[0].shape[0], self.num_procs + 1)).astype(numpy.int32)
            for i in range(1, self.num_procs):
                self.comm.send(M, dest=i, tag=11)
                self.comm.send(data, dest=i, tag=11)
                self.comm.send(indices[0][node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                self.comm.send(indices[1][node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
            indices0 = indices[0][:node_ranges[1]]
            indices1 = indices[1][:node_ranges[1]]
        else:
            M = self.comm.recv(source=0, tag=11)
            data = self.comm.recv(source=0, tag=11)
            indices0 = self.comm.recv(source=0, tag=11)
            indices1 = self.comm.recv(source=0, tag=11)
        correlations = numpy.ones((M, M), dtype=numpy.float32)
        for i in range(indices0.shape[0]):
            correlations[indices0[i], indices1[i]] = numpy.mean(data[indices0[i], :] * data[indices1[i], :])
        if self.rank == 0:
            for i in range(1, self.num_procs):
                correlations[indices[0][node_ranges[i]:node_ranges[i + 1]],
                      indices[1][node_ranges[i]:node_ranges[i + 1]]] = self.comm.recv(source=i, tag=11)
            correlations[indices[1], indices[0]] = correlations[indices]
            if not self.silent:
                print >> sys.stderr, ("\r%s\r") % (' ' * 80),
            return correlations
        else:
            self.comm.send(correlations[indices0, indices1], dest=0, tag=11)
        return None

    def orient_eigenvectors(self, fname, storage_fname=None, higher='A'):
        if self.rank != 0:
            self.eigenv = self.comm.recv(source=0, tag=11)
            self.positions = self.comm.recv(source=0, tag=11)
            return None
        if not self.silent:
            print >> sys.stderr, ("\r%s\rOrienting eigenvector scores") % (' '*80),
        chroms = self.eigenv.keys()
        scores = {}
        for chrom in chroms:
            scores[chrom] = []
        for line in open(fname):
            temp = line.rstrip('\n').split('\t')
            if line[0] == '#':
                continue
            temp[0] = temp[0].strip('chr')
            if temp[0] not in scores:
                continue
            scores[temp[0]].append((int(temp[1]), int(temp[2]), float(temp[3])))
        temp = numpy.zeros(0, dtype=numpy.float64)
        for chrom in chroms:
            scores[chrom].sort()
            scores[chrom] = numpy.array(scores[chrom], dtype=numpy.dtype([('start', numpy.int32),
                                        ('stop', numpy.int32), ('score', numpy.float64)]))
            temp = numpy.hstack((scores[chrom]['score'], temp))
        if storage_fname is not None:
            storage = h5py.File(storage_fname, 'a')
        else:
            storage = None
        temp.sort()
        minval = temp[int(0.01 * temp.shape[0])]
        maxval = temp[int(0.99 * temp.shape[0])]
        for chrom in chroms:
            scores[chrom]['score'] = numpy.maximum(minval, numpy.minimum(maxval, numpy.array(scores[chrom]['score'])))
            signs = numpy.sign(self.eigenv[chrom])
            changes = numpy.r_[0, numpy.where(signs[1:] != signs[:-1])[0] + 1, signs.shape[0]]
            bounds = numpy.zeros((changes.shape[0] - 1, 3), dtype=numpy.int32)
            bounds[:, 0] = self.positions[chrom][changes[:-1], 0]
            bounds[:, 1] = self.positions[chrom][changes[1:] - 1, 1]
            bounds[:, 2] = signs[changes[:-1]]
            A = numpy.zeros(2, dtype=numpy.float64)
            B = numpy.zeros(2, dtype=numpy.float64)
            mids = (scores[chrom]['start'] + scores[chrom]['stop']) / 2
            starts = numpy.searchsorted(mids, bounds[:, 0])
            stops = numpy.searchsorted(mids, bounds[:, 1])
            where = numpy.where(bounds[:, 2] == 1)[0]
            for i in where:
                B[0] += numpy.sum(scores[chrom]['score'][starts[i]:stops[i]])
                B[1] += numpy.sum(scores[chrom]['stop'][starts[i]:stops[i]] -
                                  scores[chrom]['start'][starts[i]:stops[i]])
            where = numpy.where(bounds[:, 2] == -1)[0]
            for i in where:
                A[0] += numpy.sum(scores[chrom]['score'][starts[i]:stops[i]])
                A[1] += numpy.sum(scores[chrom]['stop'][starts[i]:stops[i]] -
                                  scores[chrom]['start'][starts[i]:stops[i]])
            if A[0] / A[1] > B[0] / B[1] and higher != 'A':
                self.eigenv[chrom] = -self.eigenv[chrom]
                if storage is not None:
                    storage['%s.eigenv' % chrom][:] = self.eigenv[chrom]
        for i in range(1, self.num_procs):
            self.comm.send(self.eigenv, dest=i, tag=11)
            self.comm.send(self.positions, dest=i, tag=11)
        if storage is not None:
            storage.close()
        if not self.silent:
            print >> sys.stderr, ("\r%s\r") % (' '*80),
        return None

    def find_likelihood_scores(self, fname=None, max_iterations=200, burnin_iterations=20, min_reads=10000, min_distance=1000000, min_interactions=5, update_fraction=0.5, chroms=[]):
        if self.rank == 0:
            if fname is not None:
                storage = h5py.File(fname, 'a')
            else:
                storage = None
            if len(chroms) == 0:
                chroms = self.eigenv.keys()
            elif isinstance(chroms, str) and chroms in self.eigenv:
                chroms = [chroms]
            else:
                for i in range(len(chroms))[::-1]:
                    if chroms[i] not in self.eigenv:
                        del chroms[i]
            chroms.sort()
            for i in range(1, self.num_procs):
                self.comm.send(chroms, dest=i, tag=11)
        self.likelihood_scores = {}
        minbin = min_distance / self.binsize
        self.minreads = min_reads
        self.update_fraction = max(0.0, min(1.0, update_fraction))
        for chrom in chroms:
            if self.rank == 0:
                binsize = self.binsize
                if storage is not None and '%s.counts' % chrom in storage:
                    start = storage['%s.positions' % chrom][0, 0]
                    N = valid.shape[0]
                    counts = numpy.zeros((N, N), numpy.int32)
                    counts[numpy.triu_indices(N, 1)] = storage['%s.counts' % chrom][...]
                else:
                    if self.hic.binned is None:
                        chr_indices = self.hic.fends['chr_indices'][...]
                        fends = self.hic.fends['fends'][...]
                    else:
                        chr_indices = self.hic.fends['bin_indices'][...]
                        fends = self.hic.fends['bins'][...]
                    chrint = self.hic.chr2int[chrom]
                    startfend = chr_indices[chrint]
                    while self.hic.filter[startfend] == 0:
                        startfend += 1
                    start = (fends['mid'][startfend] / binsize) * binsize
                    stopfend = chr_indices[chrint + 1]
                    while self.hic.filter[stopfend - 1] == 0:
                        stopfend -= 1
                    stop = ((fends['mid'][stopfend - 1] - 1) / binsize + 1) * binsize
                    N = (stop - start) / binsize
                    start_index = self.hic.data['cis_indices'][startfend]
                    stop_index = self.hic.data['cis_indices'][stopfend]
                    data = self.hic.data['cis_data'][start_index:stop_index, :].astype(numpy.int64)
                    data = data[numpy.where(self.hic.filter[data[:, 0]] * self.hic.filter[data[:, 1]])[0], :]
                    data[:, :2] = fends['mid'][data[:, :2]]
                    data[:, :2] -= start
                    data[:, :2] /= binsize
                    counts = numpy.bincount(data[:, 0] * N + data[:, 1], weights=data[:, 2],
                                            minlength=(N * N)).reshape(N, N).astype(numpy.int32)
                counts += counts.T
                N = counts.shape[0]
                stop = start + N * binsize
                signs = numpy.sign(self.eigenv[chrom])
                changes = numpy.r_[0, numpy.where(signs[1:] != signs[:-1])[0] + 1, signs.shape[0]]
                bounds = numpy.zeros((changes.shape[0] - 1, 3), dtype=numpy.int32)
                bounds[:, 0] = self.positions[chrom][changes[:-1], 0]
                bounds[:, 1] = self.positions[chrom][changes[1:] - 1, 1]
                bounds[:, 2] = signs[changes[:-1]]
                indices = numpy.searchsorted(numpy.linspace(start, stop, N + 1), numpy.r_[bounds[:, 0], bounds[-1, 1]])
                states = numpy.zeros(N, dtype=numpy.int32)
                for i in range(bounds.shape[0]):
                    states[indices[i]:indices[i + 1]] = bounds[i, 2]
                expected = self.hic.cis_heatmap(chrom, binsize=binsize, start=start, stop=stop,
                                datatype='fend', arraytype='full', returnmapping=False, silent=True)[:, :, 1]
                indices = numpy.triu_indices(N, min_distance / binsize)
                valid = numpy.bincount(indices[0], weights=(counts[indices] > 0), minlength=N)
                valid += numpy.bincount(indices[1], weights=(counts[indices] > 0), minlength=N)
                valid = (valid >= min_interactions).astype(numpy.int32)
                expected = expected[indices]
                counts = counts[indices]
                indices = list(indices)
                indices[0] = indices[0].astype(numpy.int32)
                indices[1] = indices[1].astype(numpy.int32)
                node_ranges = numpy.round(numpy.linspace(0, indices[0].shape[0],
                                                         self.num_procs + 1)).astype(numpy.int32)
                for i in range(1, self.num_procs):
                    self.comm.send(counts[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    self.comm.send(expected[node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    self.comm.send(indices[0][node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    self.comm.send(indices[1][node_ranges[i]:node_ranges[i + 1]], dest=i, tag=11)
                    self.comm.send(valid, dest=i, tag=11)
                    self.comm.send(states, dest=i, tag=11)
                counts = counts[:node_ranges[1]]
                expected = expected[:node_ranges[1]]
                indices0 = indices[0][:node_ranges[1]]
                indices1 = indices[1][:node_ranges[1]]
                del indices
            else:
                counts = self.comm.recv(source=0, tag=11)
                expected = self.comm.recv(source=0, tag=11)
                indices0 = self.comm.recv(source=0, tag=11)
                indices1 = self.comm.recv(source=0, tag=11)
                valid = self.comm.recv(source=0, tag=11)
                states = self.comm.recv(source=0, tag=11)
            N = valid.shape[0]
            valid1 = numpy.where(valid > 0)[0]
            score_sums = numpy.zeros(N, dtype=numpy.float64)
            new_states = numpy.copy(states)
            new_scores = numpy.zeros(N, dtype=numpy.float64)
            score_changes = numpy.zeros(N, dtype=numpy.float64)
            new_score_changes = numpy.zeros(N, dtype=numpy.float64)
            where = numpy.zeros(2)
            iteration = 0
            while iteration < burnin_iterations or (iteration < max_iterations and not where.shape[0] == 0):
                prev_score_changes = score_changes
                states = new_states
                score_changes = new_score_changes
                scores = new_scores
                paramsA, paramsB, paramsAB = self._find_distance_relationship(counts, expected, states, valid,
                                                                              indices0, indices1, minbin)
                new_states, new_scores = self._optimize_bound_iteration(counts, expected, states, valid, indices0,
                                                                        indices1, paramsA, paramsB, paramsAB)
                if iteration >= burnin_iterations:
                    score_sums += new_scores
                where = valid1[numpy.where(states[valid1] != new_states[valid1])[0]]
                new_score_changes = new_scores[where]
                if where.shape[0] > 0 and self.rank == 0 and not self.silent:
                    print >> sys.stderr, ("\nIteration: %i\tDifferences: %i\tDelta: %0.6f     ") % (iteration, where.shape[0], numpy.mean(numpy.abs(new_scores[where] - scores[where]))),
                iteration += 1
                if where.shape[0] == 0:
                    score_sums = new_scores
                    iteration = burnin_iterations
                    break
                if numpy.array_equal(prev_score_changes, new_score_changes):
                    score_sums = scores + new_scores
                    iteration = burnin_iterations + 1
                    break
            score_sums /= iteration - burnin_iterations + 1
            self.likelihood_scores[chrom] = score_sums[valid1]
            if self.rank == 0 and storage is not None:
                if '%s.likelihood' % chrom not in storage:
                    storage.create_dataset(name='%s.likelihood' % chrom, data=self.likelihood_scores[chrom])
                else:
                    storage['%s.likelihood' % chrom][:] = self.likelihood_scores[chrom]
        if self.rank == 0 and storage is not None:
            storage.close()
        return None

    def write_likelihood_scores(self, fname):
        if self.rank > 0:
            return None
        outfile = open(fname, 'w')
        chroms = self.likelihood_scores.keys()
        chroms.sort()
        for chrom in chroms:
            for i in range(self.likelihood_scores[chrom].shape[0]):
                print >> outfile, "%s\t%i\t%i\t%f" % (chrom, self.positions[chrom][i, 0], self.positions[chrom][i, 1],
                                                      self.likelihood_scores[chrom][i])
        outfile.close()
        return None

    def write_eigen_scores(self, fname):
        if self.rank > 0:
            return None
        outfile = open(fname, 'w')
        chroms = self.eigenv.keys()
        chroms.sort()
        for chrom in chroms:
            for i in range(self.eigenv[chrom].shape[0]):
                print >> outfile, "%s\t%i\t%i\t%f" % (chrom, self.positions[chrom][i, 0], self.positions[chrom][i, 1],
                                                      self.eigenv[chrom][i])
        outfile.close()
        return None

    def _find_distance_relationship(self, counts, expected, states, valid, indices0, indices1, minbin):
        N = states.shape[0]
        dist_counts = numpy.zeros((N, 3), dtype=numpy.int32)
        dist_expected = numpy.zeros((N, 3), dtype=numpy.float64)
        binsizes = numpy.zeros((N, 3), dtype=numpy.int32)
        _hic_domains.find_distance_parameters(
            counts,
            expected,
            valid,
            states,
            indices0,
            indices1,
            dist_counts,
            dist_expected,
            binsizes)
        if self.rank == 0:
            for i in range(1, self.num_procs):
                dist_counts += self.comm.recv(source=i, tag=11)
                dist_expected += self.comm.recv(source=i, tag=11)
                binsizes += self.comm.recv(source=i, tag=11)
            distances = numpy.log(numpy.r_[0.5, numpy.arange(1, N)].astype(numpy.float64))
            params = []
            for i in range(3):
                data = []
                pos = minbin
                c_sum = 0
                e_sum = 0.0
                d_sum = 0.0
                n_sum = 0
                while pos < N:
                    c_sum += dist_counts[pos, i]
                    e_sum += dist_expected[pos, i]
                    n_sum += binsizes[pos, i]
                    d_sum += distances[pos] * binsizes[pos, i]
                    if c_sum >= self.minreads:
                        data.append([d_sum / n_sum, numpy.log(c_sum / e_sum)])
                        c_sum = 0
                        e_sum = 0.0
                        d_sum = 0.0
                        n_sum = 0
                    pos  += 1
                if n_sum > 0 and e_sum > 0 and c_sum > self.minreads:
                    data.append([d_sum / n_sum, numpy.log(c_sum / e_sum)])
                data = numpy.array(data, dtype=numpy.float64)
                slopes = (data[1:, 1] - data[:-1, 1]) / (data[1:, 0] - data[:-1, 0])
                intercepts = data[1:, 1] - data[1:, 0] * slopes
                indices = numpy.searchsorted(data[1:-1, 0], distances)
                distance_parameters = numpy.exp(distances * slopes[indices] + intercepts[indices]).astype(numpy.float64)
                params.append(distance_parameters)
            for i in range(1, self.num_procs):
                self.comm.send(params, dest=i, tag=11)
        else:
            self.comm.send(dist_counts, dest=0, tag=11)
            self.comm.send(dist_expected, dest=0, tag=11)
            self.comm.send(binsizes, dest=0, tag=11)
            params = self.comm.recv(source=0, tag=11)
        return params

    def _optimize_bound_iteration(self, counts, expected, states, valid, indices0, indices1, paramsA, paramsB, paramsAB):
        N = valid.shape[0]
        scores = numpy.zeros((N, 2), dtype=numpy.float64)
        _hic_domains.optimize_compartment_bounds(
            counts,
            expected,
            valid,
            states,
            indices0,
            indices1,
            paramsA,
            paramsB,
            paramsAB,
            scores)
        if self.rank == 0:
            for i in range(1, self.num_procs):
                scores += self.comm.recv(source=i, tag=11)
            scores = scores[:, 0] - scores[:, 1]
            temp_states = numpy.sign(scores)
            valid1 = numpy.where(valid > 0)[0]
            where2 = valid1[numpy.where(states[valid1] != temp_states[valid1])[0]]
            temp = numpy.abs(numpy.copy(scores[where2]))
            temp.sort()
            new_states = numpy.copy(states)
            if temp.shape[0] > 0:
                cutoff = temp[int(numpy.floor((1.0 - self.update_fraction) * temp.shape[0]))]
                where3 = where2[numpy.where(numpy.abs(scores[where2]) >= cutoff)[0]]
                new_states[where3] = numpy.copy(temp_states[where3])
            scores[numpy.where(valid == 0)] = numpy.nan
            for i in range(1, self.num_procs):
                self.comm.Send(scores, dest=i, tag=13)
                self.comm.Send(new_states, dest=i, tag=13)
        else:
            self.comm.send(scores, dest=0, tag=11)
            scores = numpy.zeros(N, dtype=numpy.float64)
            new_states = numpy.zeros(N, dtype=numpy.int32)
            self.comm.Recv(scores, source=0, tag=13)
            self.comm.Recv(new_states, source=0, tag=13)
        return new_states, scores

    def find_clusters(self):
        means = numpy.zeros(2, dtype=numpy.float64)
        counts = numpy.zeros(2, dtype=numpy.int32)
        observations = []
        for chrom in self.eigenv:
            pos = numpy.where(self.eigenv[chrom] >= 0)[0]
            means[0] += numpy.sum(self.eigenv[chrom][pos])
            counts[0] += pos.shape[0]
            neg = numpy.where(self.eigenv[chrom] < 0)[0]
            means[1] += numpy.sum(self.eigenv[chrom][neg])
            counts[1] += neg.shape[0]
            observations.append(numpy.copy(self.eigenv[chrom]).astype(numpy.float64))
        means /= counts
        distributions = [
            [[0.33, means[0] * 0.5, means[0] / 4.0],
             [0.34, means[0], means[0] / 4.0],
             [0.33, means[0] * 1.5, means[0] / 4.0]],
            [[0.33, means[1] * 0.5, -means[1] / 4.0],
             [0.34, means[1], -means[1] / 4.0],
             [0.33, means[1] * 1.5, -means[1] / 4.0]]]
        transitions = [[0.99, 0.01], [0.01, 0.99]]
        hmm = HMM(num_states=2, num_distributions=3, pi=[0.5, 0.5], distributions=distributions,
                  transitions=transitions, seed=2001)
        hmm.train(observations)
        self.clusters = {}
        for chrom in self.eigenv:
            self.clusters[chrom] = hmm.find_path(self.eigenv[chrom].astype(numpy.float64))[0]

    def write_compartments(self, out_fname):
        outfile = open(out_fname, 'w')
        chroms = self.clusters.keys()
        chroms.sort()
        for chrom in chroms:
            bounds = numpy.r_[0, numpy.where(self.clusters[chrom][1:] != self.clusters[chrom][:-1])[0] + 1,
                              self.clusters[chrom].shape[0]]
            for i in range(bounds.shape[0] - 1):
                j = bounds[i]
                k = bounds[i + 1] - 1
                print >> outfile, "%s\t%i\t%i\t%i" % (chrom, self.positions[chrom][j, 0], self.positions[chrom][k, 1],
                                                      self.clusters[chrom][j])
        outfile.close()

class Boundary( object ):
    """
    """

    def __init__(self, hic, silent=False):
        self.hic = hic
        self.silent = silent

    def find_band_boundaries(self, minband=10000, maxband=82000, bandstep=4000, step=1000, minwidth=6000,
                             maxwidth=40000, chroms=[]):
        minband = (minband / step) * step
        maxband = (maxband / step) * step
        n = (maxband - minband) / bandstep + 1
        self.bands = numpy.zeros((n, 2), dtype=numpy.int32)
        self.bands[:, 0] = numpy.arange(n) * bandstep + minband
        self.bands[:, 1] = self.bands[:, 0] + (numpy.round(10 ** (numpy.linspace(numpy.log10(minwidth),
                           numpy.log10(maxwidth), n))).astype(numpy.int32) / (step * 2)) * step * 2
        self.step = step
        self.band_scores = {}
        if isinstance(chroms, str):
            chroms = [chroms]
        if len(chroms) == 0:
            chroms = list(self.hic.fends['chromosomes'][...])
        for chrom in chroms:
            temp = self.find_band_scores(chrom)
            if temp is not None:
                self.band_scores[chrom] = temp
            temp = numpy.zeros((self.band_scores[chrom].shape[0], self.band_scores[chrom].shape[1], 2), dtype=numpy.float32)
            where = numpy.where(numpy.logical_not(numpy.isnan(self.band_scores[chrom])))
            temp[where[0], where[1], 0] = self.band_scores[chrom][where]
            temp[where[0], where[1], 1] = 1
            img = plotting.plot_full_array(temp[:, ::-1, :], symmetricscaling=False, logged=False)
            img.save('test.png')

    def __getitem__(self, key):
        """Dictionary-like lookup."""
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        """Dictionary-like value setting."""
        self.__dict__[key] = value
        return None

    def find_BI_scores(self, binsize=10000, width=50000, window=250000, chroms=[]):
        if isinstance(chroms, str):
            chroms = [chroms]
        if len(chroms) == 0:
            chroms = list(self.hic.fends['chromosomes'][...])    
        self.binsize = binsize
        self.width = width
        self.window = window
        widthB = width / binsize
        windowB = window / binsize
        if 'scores' not in self.__dict__.keys():
            self.scores = {}
        for chrom in chroms:
            chrint = self.hic.chr2int[chrom]
            if 'binned' in self.hic.__dict__.keys() and self.hic['binned'] is not None:
                fends = self.hic.fends['bins'][...]
                chr_indices = self.hic.fends['bin_indices'][...]
            else:
                fends = self.hic.fends['fends'][...]
                chr_indices = self.hic.fends['chr_indices'][...]
            start_fend = chr_indices[chrint]
            stop_fend = chr_indices[chrint + 1]
            while self.hic.filter[start_fend] == 0:
                start_fend += 1
            while self.hic.filter[stop_fend - 1] == 0:
                stop_fend -= 1
            start = (fends['mid'][start_fend] / binsize) * binsize
            stop = ((fends['mid'][stop_fend - 1] - 1) / binsize + 1) * binsize
            N = (stop - start) / binsize
            temp = self.hic.cis_heatmap(chrom, binsize=binsize, datatype='expected', arraytype='compact',
                                        maxdistance=window, returnmapping=False)
            if temp is None:
                continue
            temp = temp[:, :, 1]
            expected = numpy.zeros((N, N), dtype=numpy.float32)
            for i in range(N - 1):
                expected[i, (i + 1):min(N, i + temp.shape[1] + 1)] = temp[i, :min(N - i - 1, temp.shape[1])]
            expected += expected.T
            start_index = self.hic.data['cis_indices'][start_fend]
            stop_index = self.hic.data['cis_indices'][stop_fend]
            data = self.hic.data['cis_data'][start_index:stop_index, :].astype(numpy.int64)
            data = data[numpy.where(self.hic.filter[data[:, 0]] * self.hic.filter[data[:, 1]])[0], :]
            data[:, :2] = fends['mid'][data[:, :2]]
            data[:, :2] -= start
            data[:, :2] /= binsize
            counts = numpy.bincount(data[:, 0] * N + data[:, 1], weights=data[:, 2],
                                    minlength=(N * N)).reshape(N, N).astype(numpy.int32)
            counts += counts.T
            M = N - widthB + 1
            hm = numpy.zeros((M, M, 2), dtype=numpy.float32)
            for i in range(widthB):
                hm[:, :, 0] += counts[i:(M + i), :M]
                hm[:, :, 1] += expected[i:(M + i), :M]
            hm2 = numpy.copy(hm)
            for i in range(1, widthB):
                hm[:, :(M - i), :] += hm2[:, i:, :]
                for j in range(widthB):
                    hm[:, (M - i):, 0] += counts[j:(M + j), M:(M + i)]
                    hm[:, (M - i):, 1] += expected[j:(M + j), M:(M + i)]
            del counts
            del expected
            del hm2
            N = hm.shape[1]
            mapping = numpy.zeros((N, 2), dtype=numpy.int32)
            mapping[:, 0] = numpy.arange(N) * binsize + start
            mapping[:, 1] = numpy.arange(1, N + 1) * binsize + start
            scores = []
            for i in range(widthB, N):
                if i - windowB < 0:
                    start = i - (i / binsize) * binsize
                else:
                    start = i - windowB
                if i + windowB > N:
                    stop = i + ((N - i) / binsize) * binsize - widthB + 1
                else:
                    stop = i + windowB - widthB + 1
                valid = numpy.r_[
                            numpy.where((hm[i - widthB, start:(i - 2 * widthB + 1), 0] > 0) *
                                        (hm[i, start:(i - 2 * widthB + 1), 0] > 0))[0] + start,
                            numpy.where((hm[i - widthB, (i + widthB):stop, 0] > 0) *
                                        (hm[i, (i + widthB):stop, 0] > 0))[0] + i + widthB]
                if valid.shape[0] > 5:
                    set1 = numpy.log(hm[i - widthB, valid, 0] / hm[i - widthB, valid, 1])
                    set2 = numpy.log(hm[i, valid, 0] / hm[i, valid, 1])
                    scores.append((mapping[i, 0], numpy.corrcoef(set1, set2)[0, 1]))
            scores.sort()
            scores = numpy.array(scores, dtype=numpy.dtype([('position', numpy.int32), ('score', numpy.float32)]))
            self.scores[chrom] = scores
        return None

    def write_scores(self, fname):
        chroms = self.scores.keys()
        chroms.sort()
        outfile = open(fname, 'w')
        for chrom in chroms:
            for i in self.scores[chrom].shape[0]:
                if not numpy.isnan(self.scores[chrom]['score'][i]):
                    print >> outfile, "%s\t%i\t%f" % (chrom, self.scores[chrom]['position'][i],
                                                      self.scores[chrom]['score'][i])
        outfile.close()

    def find_band_scores(self, chrom):
        temp = self.hic.cis_heatmap(chrom, binsize=self.step, datatype='fend', arraytype='compact',
                                    maxdistance=self.bands[-1, 1], returnmapping=True)
        if temp is None:
            return None
        hm, mapping = temp
        scores = numpy.zeros((mapping.shape[0] + 1, self.bands.shape[0]), dtype=numpy.float32)
        scores.fill(numpy.nan)
        for i in range(scores.shape[1]):
            _hic_domains.find_band_score(hm, scores, self.bands[i, 0] / self.step, self.bands[i, 1] / self.step, i)
        return scores






