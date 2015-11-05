#!/usr/bin/env python

"""
This is a module contains scripts for finding domains in HiC interaction data.

Concepts
--------

These functions rely on the :class:`HiC` class in conjunction with the :class:`Fend` and :class:`HiCData` classes.


API Documentation
-----------------
"""

import os
import sys

import numpy
import h5py
try:
    from mpi4py import MPI
except:
    pass

import libraries._hic_tads as _hic_tads
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
    def learn_TADs(self, maxsize=2000000, maxtreesize=25, p=3, q=12, gamma=0.5, chroms=[], binsize=10000,
                   minsize=50000):
        if isinstance(chroms, list) and len(chroms) == 0:
            chroms = list(self.hic.fends['chromosomes'][...])
        elif isinstance(chroms, str):
            chroms = [chroms]
        if maxsize % binsize != 0:
            print >> sys.stderr, ("Maximum TAD size must be a multiple of bin size.\n"),
            return
        self.parameters = {
            'maxsize': int(maxsize),
            'minsize': int(minsize),
            'maxtreesize': int(maxtreesize),
            'p': int(p),
            'q': int(q),
            'gamma': float(gamma),
            'binsize': int(binsize),
        }
        self.chromosomes = numpy.array(chroms)
        for chrom in chroms:
            self[chrom] = self.find_TAD(chrom)

    def find_TAD(self, chrom):
        p, q, gamma = self.parameters['p'], self.parameters['q'], self.parameters['gamma']
        maxdist = max(max(p, q) * self.parameters['binsize'], self.parameters['maxsize'])
        data = self.hic.cis_heatmap(chrom, binsize=self.parameters['binsize'], arraytype='compact',
                                    datatype='fend', include_diagonal=False, maxdistance=maxdist)
        #where = numpy.where(data[:, :, 0] == 0)
        #data[where[0], where[1], 1] = 0
        #where = numpy.where(data[:, :, 0] > 0)
        #data[where[0], where[1], 0] = numpy.log(data[where[0], where[1], 0] / data[where[0], where[1], 1])
        #data[where[0], where[1], 1] = 1
        maxbins = self.parameters['maxsize'] / self.parameters['binsize']
        minbins = self.parameters['minsize'] / self.parameters['binsize']
        n = data.shape[0]
        m = data.shape[1]
        print >> sys.stderr, ("\rFinding BI scores for chromosome %s...") % (chrom),
        BIs = numpy.zeros((n, m, 2), dtype=numpy.float32)
        BIs.fill(-numpy.inf)
        _hic_tads.find_BIs(data, BIs, p, minbins)
        where = numpy.where(BIs[:, :, 1] == 0)
        print numpy.amin(BIs[where[0], where[1], 0]), numpy.mean(BIs[where[0], where[1], 0]), numpy.amax(BIs[where[0], where[1], 0])        
        where = numpy.where(BIs[:, :, 1] == 1)
        print numpy.amin(BIs[where[0], where[1], 0]), numpy.mean(BIs[where[0], where[1], 0]), numpy.amax(BIs[where[0], where[1], 0])
        BI_scores = numpy.zeros((n, maxbins - minbins + 1), dtype=numpy.float32)
        for i in range(minbins, maxbins + 1):
            BI_scores[:(n - i + 1), i - minbins] = BIs[:(n - i + 1), 0] * BIs[(i - 1):, 1]
        where = numpy.where(BI_scores <= 0.0)
        BI_scores[where] = -numpy.inf
        where = numpy.where(BI_scores > -numpy.inf)
        BI_scores[where] = BI_scores[where] ** gamma
        where = numpy.where(data[:, :, 0] == 0)
        data[where[0], where[1], 1] = 0
        where = numpy.where(data[:, :, 0] > 0)
        data[where[0], where[1], 0] = numpy.log(data[where[0], where[1], 0] / data[where[0], where[1], 1])
        data[where[0], where[1], 1] = 1
        print >> sys.stderr, ("\rFinding TAD parameters for chromosome %s...") % (chrom),
        scores = numpy.zeros((n, maxbins + 1), numpy.float32)
        scores.fill(numpy.inf)
        std_params = numpy.zeros((n, maxbins + 1, 3), dtype=numpy.float32)
        _hic_tads.find_initial_TAD_std_params(data, BIs, scores, std_params, maxbins, minbins, gamma)
        paths = numpy.zeros((n, maxbins - minbins + 2), dtype=numpy.int32)
        paths.fill(-1)
        path_scores = numpy.zeros((n + 1, maxbins - minbins + 2, 2), dtype=numpy.float32)
        #path_scores[:, :, 0].fill(numpy.inf)
        final_path = numpy.zeros(n, dtype=numpy.int32)
        _hic_tads.find_TAD_path(scores, paths, path_scores, final_path, minbins, maxbins)
        print list(final_path)
        for i in range(10):
            print list(path_scores[i, :, 0])
        where = numpy.where(numpy.abs(scores) < numpy.inf)
        print numpy.amin(scores[where]), numpy.mean(scores[where]), numpy.amax(scores[where])
        where = numpy.where(numpy.abs(BI_scores) < numpy.inf)
        print numpy.amin(BI_scores[where]), numpy.mean(BI_scores[where]), numpy.amax(BI_scores[where])
        #numpy.savetxt('temp.txt',paths, fmt="%i", delimiter='\t')
        #subTAD_scores = numpy.zeros((n, maxbins - minbins + 1, self.parameters['maxtreesize']), dtype=numpy.float32)
        #subTAD_params = numpy.zeros((n, maxbins - minbins + 1, self.parameters['maxtreesize'], 3), dtype=numpy.float32)
        #_hic_tads.find_TAD_subparts(subTAD_scores, subTAD_params, BIs, std_params, minbins, gamma)
        #where = numpy.where(std_params[:, :, 0] >= 3)
        #errors = numpy.zeros((std_params.shape[0], std_params.shape[1], 2), dtype=numpy.float32)
        ##errors[where[0], where[1], 0] = (std_params[where[0], where[1], 2] / std_params[where[0], where[1], 0]) - (std_params[where[0], where[1], 1] / std_params[where[0], where[1], 0]) ** 2 - numpy.maximum(0, BIs[where[0], 0] * BIs[where[0] + where[1] + minbins - 1, 1])**gamma
        #print numpy.amax(where[0] + where[1] + minbins), BIs.shape
        #where1 = numpy.where(errors[where[0], where[1], 0] < 0)
        ##errors[where[0][where1], where[1][where1], 1] = 1
        #print numpy.mean(errors), numpy.amax(errors)
        #print numpy.amax(BIs), numpy.mean(BIs)
        #_hic_tads.find_betadeltas(data, betas, deltas, fits, errors, maxbins)
        #where = numpy.where(data[:, :, 0] == 0)
        #data[where[0], where[1], 1] = 0
        #where = numpy.where(data[:, :, 0] > 0)
        #data[where[0], where[1], 0] = numpy.log(data[where[0], where[1], 0] / data[where[0], where[1], 1])
        #data[where[0], where[1], 1] = 1
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                if BIs[i] < 0 or BIs[j] < 0:
                    fits[i,j] = numpy.inf
                else:
                    fits[i,j] -= gamma * (BIs[i] + BIs[j])
        print >> sys.stderr, ("\rFinding TAD trees for chromosome %s...") % (chrom),
        scores = numpy.zeros((n, n, self.parameters['maxtreesize']), dtype=numpy.float32)
        _hic_tads.build_TAD_trees(data, fits, deltas, betas, errors, scores, maxbins)
        alldata = numpy.zeros((n, n, 2), dtype=numpy.float32)
        for i in range(data.shape[0] - 1):
            alldata[(i + 1):min(i + data.shape[1] + 1, n), i, :] = data[i, :min(data.shape[1], n - i - 1), :]
        for i in range(scores.shape[0] - 1 - minbins):
            where = numpy.where(numpy.abs(scores[i, :min(scores.shape[1], n - i - minbins - 1)]) < numpy.inf)[0]
            alldata[i, where + i + minbins, 0] = scores[i, where]
            alldata[i, where + i + minbins, 1] = 1
        for i in range(BIs.shape[0] - 1):
            alldata[i, (i + 1):min(i + data.shape[1] + 1, n), :] = BIs[i, :min(data.shape[1], n - i - 1), :]
            #where = numpy.where(numpy.abs(BIs[i, :min(BIs.shape[1], n - i - 1), 1]) < numpy.inf)[0]
            #alldata[i, where + i, 0] = BIs[i, where, 1]
            #alldata[i, where + i, 1] = 1
        indices = numpy.triu_indices(n, 1)
        where = numpy.where(alldata[indices[0], indices[1], 1] > 0)[0]
        alldata[indices[0][where], indices[1][where], 0] -= numpy.amin(alldata[indices[0][where], indices[1][where], 0])
        alldata[indices[0][where], indices[1][where], 0] /= numpy.amax(alldata[indices[0][where], indices[1][where], 0])
        where = numpy.where(alldata[indices[1], indices[0], 1] > 0)[0]
        alldata[indices[1][where], indices[0][where], 0] -= numpy.amin(alldata[indices[1][where], indices[0][where], 0])
        alldata[indices[1][where], indices[0][where], 0] /= numpy.amax(alldata[indices[1][where], indices[0][where], 0])
        for i in range(n):
            if final_path[i] != 0:
                indices = numpy.triu_indices(final_path[i], 1)
                temp = alldata[i + final_path[i] - 1, i, 0]
                alldata[indices[1] + i, indices[0] + i, :] = 1.0
                alldata[i + final_path[i] - 1, i, 0] = temp
        img = plotting.plot_full_array(alldata, symmetricscaling=False, logged=False)
        img.save('BIs_%s.png' % chrom)
    """

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
            #temp = self.hic.cis_heatmap(chrom, binsize=binsize, maxdistance=self.maxsize, datatype='fend',
            #                                     arraytype='compact', returnmapping=True, start=mapping[0, 0], stop=mapping[-1, 1])
            #if temp is None:
            #    continue
            #data, mapping = temp[:2]
            data = numpy.zeros((heatmap.shape[0], maxbins - 1, 2), dtype=numpy.float32)
            for i in range(heatmap.shape[0] - 1):
                data[i, :min(data.shape[1], data.shape[0] - i - 1), :] = heatmap[i, (i + 1):min(data.shape[1] + i + 1, data.shape[0]), :]
            #heatmap = plotting.plot_compact_array(data, logged=True, symmetricscaling=False, silent=True)
            where = numpy.where(data[:, :, 1] > 0)
            data[where[0], where[1], 0] /= data[where[0], where[1], 1]
            scores = numpy.zeros((data.shape[0], data.shape[1]), dtype=numpy.float32)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding arrowhead transformation for chromosome %s...") % (' ' * 80, chrom),
            _hic_tads.find_arrowhead_transformation(data, scores, maxbins)
            sums = numpy.zeros(data.shape, dtype=numpy.float32)
            signs = numpy.zeros(data.shape, dtype=numpy.float32)
            variances = numpy.zeros((data.shape[0], data.shape[1], 2, 2), dtype=numpy.float32)
            domain_scores = numpy.zeros(scores.shape, dtype=numpy.float32)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding arrowhead scoring for chromosome %s...") % (' ' * 80, chrom),
            _hic_tads.find_arrowhead_scores(scores, sums, signs, variances, domain_scores, minbins)
            """
            where = numpy.where(heatmap[:, :, 0] == 0)
            heatmap[where[0], where[1], 1] = 0
            where = numpy.where(heatmap[:, :, 0] > 0)
            heatmap[where[0], where[1], 0] = numpy.log(heatmap[where[0], where[1], 0] / heatmap[where[0], where[1], 1])
            heatmap[where[0], where[1], 0] -= numpy.amin(heatmap[where[0], where[1], 0])
            heatmap[where[0], where[1], 0] /= numpy.amax(heatmap[where[0], where[1], 0])
            heatmap[where[0], where[1], 1] = 1
            scores = numpy.zeros((domain_scores.shape[0], domain_scores.shape[1], 2), dtype=numpy.float32)
            scores[:, :, 0] = domain_scores
            scores[:, :, 0] -= numpy.amin(domain_scores)
            scores[:, :, 0] /= numpy.amax(scores[:, :, 0])
            scores[:, :, 1] = 1
            temp = hic_binning._compact_to_upper(scores)
            indices = numpy.triu_indices(heatmap.shape[0], 1)
            heatmap[indices[1], indices[0], :] = temp
            img = plotting.plot_full_array(heatmap, symmetricscaling = False, logged=False)
            img.save("arrowhead_scoring_%s.png" % chrom)
            """

            path = numpy.zeros(data.shape[0], dtype=numpy.int32)
            path_scores = numpy.zeros(data.shape[0], dtype=numpy.float64)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding optimal domains for chromosome %s...") % (' ' * 80, chrom),
            _hic_tads.find_arrowhead_path(domain_scores, path, path_scores, minbins, maxbins)
            i = path.shape[0] - 1
            domains = []
            while i > 0:
                if path[i] != 1:
                    domains.append([i - path[i], i])
                    self.TADs[chrom].append([mapping[domains[-1][0], 0], mapping[domains[-1][1], 1]])
                i -= path[i]

            """
            c = pyx.canvas.canvas()
            c.insert(pyx.bitmap.bitmap(0, 0, img, width=20))
            span = float(mapping[-1, 1] - mapping[0, 0])
            for i in range(len(self.TADs[chrom])):
                start = (self.TADs[chrom][i][0] - mapping[0, 0]) / span * 20.0
                stop = (self.TADs[chrom][i][1] - mapping[0, 0]) / span * 20.0
                c.stroke(pyx.path.rect(start, 20 - start, stop-start, start-stop), [pyx.style.linewidth.THIN])
            c.writePDFfile('arrowhead_scoring_%s.pdf' % chrom)

            c = pyx.canvas.canvas()
            c0 = pyx.canvas.canvas([pyx.canvas.clip(pyx.path.path(
                pyx.path.moveto(0, 20),
                pyx.path.lineto(20, 20),
                pyx.path.lineto(20, 0),
                pyx.path.closepath()))])
            c1 = pyx.canvas.canvas([pyx.canvas.clip(pyx.path.path(
                pyx.path.moveto(0, 20),
                pyx.path.lineto(0, 0),
                pyx.path.lineto(20, 0),
                pyx.path.closepath()))])
            width = 20.0 * (2.0 ** 0.5)
            offset = data.shape[1] / float(data.shape[0]) * width / 2.0
            c = pyx.canvas.canvas([pyx.canvas.clip(pyx.path.rect(0, -offset, width, 2.0 * offset))])
            c0.insert(pyx.bitmap.bitmap(0, 0, heatmap, width=20))
            c.insert(c0, [pyx.trafo.translate(0, -20), pyx.trafo.rotate(45)])
            c1.insert(pyx.bitmap.bitmap(0, 0, domain_img, width=20))
            c.insert(c1, [pyx.trafo.translate(0, -20), pyx.trafo.rotate(45)])
            for d in domains:
                x0 = d[0] / float(data.shape[0]) * width
                x1 = d[1] / float(data.shape[0]) * width
                y = (x1 - x0) / 2.0
                c.stroke(pyx.path.path(
                    pyx.path.moveto(x0, 0),
                    pyx.path.lineto(x0 + y, y),
                    pyx.path.lineto(x1, 0),
                    pyx.path.lineto(x0 + y, -y),
                    pyx.path.closepath()),
                    [pyx.style.linewidth.Thin])
            c.writePDFfile('arrowhead_heatmap_%s.pdf' % chrom)
            """
        if not self.silent:
            print >> sys.stderr, ("\r%s\rFinished finding TADs\n") % (' ' * 80),

    def write_TADs(self, fname):
        output = open(fname, 'w')
        chroms = self.TADs.keys()
        chroms.sort()
        for chrom in chroms:
            self.TADs[chrom].sort()
            for domain in self.TADs[chrom]:
                print >> output, "%s\t%i\t%i" % (chrom, domain[0], domain[1])
        output.close()







