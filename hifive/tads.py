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

import libraries._hic_tads as _hic_tads
import hic_binning

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
            #temp_data, mapping = self.hic.cis_heatmap(chrom, binsize=binsize * 16, datatype='fend', arraytype='full', returnmapping=True)
            temp = self.hic.cis_heatmap(chrom, binsize=binsize, datatype='enrichment', arraytype='compact',
                                        maxdistance=maxsize, returnmapping=True)#, start=mapping[0, 0], stop=mapping[-1, 1])
            if temp is None:
                continue
            data, mapping = temp
            print >> sys.stderr, ("\rFinding BI scores for chromosome %s...") % (chrom),
            BIs = numpy.zeros((data.shape[0], maxbins + 1, 2), dtype=numpy.float32)
            BIs.fill(-numpy.inf)
            _hic_tads.find_BIs(data, BIs, minbins, maxbins, width)
            path = numpy.zeros(data.shape[0], dtype=numpy.int32)
            path_scores = numpy.zeros(data.shape[0], dtype=numpy.float64)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding optimal domains for chromosome %s...") % (' ' * 80, chrom),
            _hic_tads.find_BI_path(BIs, path, path_scores, minbins, maxbins)
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
            _hic_tads.find_arrowhead_transformation(data, scores, maxbins)
            sums = numpy.zeros(data.shape, dtype=numpy.float32)
            signs = numpy.zeros(data.shape, dtype=numpy.float32)
            variances = numpy.zeros((data.shape[0], data.shape[1], 2, 2), dtype=numpy.float32)
            domain_scores = numpy.zeros(scores.shape, dtype=numpy.float32)
            if not self.silent:
                print >> sys.stderr, ("\r%s\rFinding arrowhead scoring for chromosome %s...") % (' ' * 80, chrom),
            _hic_tads.find_arrowhead_scores(scores, sums, signs, variances, domain_scores, minbins)
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







