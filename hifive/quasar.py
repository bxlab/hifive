#!/usr/bin/env python

"""A class for determining HiC data quality."""

import os
import sys
from math import ceil, floor

import numpy
import h5py
from scipy.optimize import curve_fit
try:
    from mpi4py import MPI
except:
    pass
try:
    from pyx import *
    unit.set(defaultunit="cm")
    text.set(mode="latex")
except:
    pass

import libraries._quasar as _quasar


class Quasar(object):

    """This class performs subsampling and QuASAR transformations for calculating HiC quality.

    .. note::
      This class is also available as hifive.Quasar

    When initialized, this class creates an h5dict in which to store all data associated with this object.
    
    :param filename: The file name of the h5dict to store the QuASAR-transformed data in.
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`Quasar` class object.
    """

    def __init__(self, filename, mode='a', silent=False):
        """Create a :class:`Quasar` object."""
        try:
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.num_procs = self.comm.Get_size()
        except:
            self.comm = None
            self.rank = 0
            self.num_procs = 1
        self.file = os.path.abspath(filename)
        self.hic_fname = None
        self.silent = silent
        self.filetype = 'quasar'
        self.strict_qcutoff = 0.05798399
        self.loose_qcutoff = 0.04345137
        self.strict_rcutoff = 0.85862067
        self.loose_rcutoff = 0.80026913
        if mode != "w":
            self.load(mode)
        else:
            if self.rank == 0:
                self.storage = h5py.File(self.file, 'w')
            else:
                self.storage = None
        return None

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

    def save(self):
        """
        Save analysis parameters to h5dict.

        :returns: None
        """
        if self.rank != 0:
            return None
        for key in self.__dict__.keys():
            if key in ['file', 'storage', 'silent', 'comm', 'rank', 'num_procs']:
                continue
            elif self[key] is None:
                continue
            elif isinstance(self[key], numpy.ndarray):
                if key in self.storage:
                    del self.storage[key]
                self.storage.create_dataset(key, data=self[key])
            elif isinstance(self[key], list):
                if isinstance(self[key][0], numpy.ndarray):
                    for i in range(len(self[key])):
                        self.storage.create_dataset("%s.%s" % (key, chroms[i]), data=self[key][i])
            elif not isinstance(self[key], dict):
                self.storage.attrs[key] = self[key]
        return None

    def load(self, mode='a'):
        """
        Load data from h5dict specified at object creation.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :param mode: The mode to open the h5dict with.
        :type mode: str.

        :returns: None
        """
        if self.rank != 0:
            self.storage = None
            return None
        self.storage = h5py.File(self.file, mode)
        for key in self.storage.keys():
            if key.split('.')[0] in ['valid', 'dist', 'corr']:
                continue
            self[key] = numpy.copy(self.storage[key])
        for key in self.storage['/'].attrs.keys():
            self[key] = self.storage['/'].attrs[key]
        return None

    def close(self):
        """
        Close h5dict file.

        :returns: None
        """
        if self.rank == 0:
            self.storage.close()
        return None

    def find_transformation(self, hic, chroms=[], resolutions=[1000000, 200000, 40000, 10000],
        coverages=[0, 40000000, 20000000, 10000000, 5000000, 2000000, 1000000], seed=None):
        """
        Find QuASAR transformation from the specified HiC project.

        :param hic: The HiC project from which to calculate QuASAR transformations from. If this function has been previously called with a different HiC project, the current transformed matrices will be deleted prior to calculating new matrices.
        :type hic: class:`HiC` class object.
        :param chroms: A list of chromosome names to calculate transformed matrices from. If this is an empty list, all chromosomes from the HiC object will be used.
        :type mode: list
        :param resolutions: A list of binning resolutions to find transformed matrices for.
        :type resolutions: list
        :param coverages: A list of cis read counts to downsample to prior to finding transformed matrices. A value of 0 indicates to use all reads. Coverages are calculated across only chromosomes specified in the 'chroms' argument.
        :type coverages: list
        :param seed: An integer to use as the initialization value for the random number generator.

        :returns: :class:`Quasar` class object.
        """
        if self.rank == 0:
            coverages = numpy.array(coverages, dtype=numpy.int64)
            resolutions = numpy.array(resolutions, dtype=numpy.int64)
            resolutions.sort()
            hic_fname = hic.file
            if self.hic_fname is not None and self.hic_fname != hic_fname:
                for key in self.storage.keys():
                    if key.split('.')[0] in ['valid', 'dist', 'corr']:
                        del self.storage[key]
                    if 'chromosomes' in self.storage:
                        del self.storage['chromosomes']
            self.hic_fname = hic_fname
            if seed is not None:
                RNG = numpy.random.RandomState(seed=seed)
            else:
                RNG = numpy.random.RandomState()

            # load partition information
            if 'binned' in hic.__dict__ and hic.binned is not None:
                temp_mids = hic.fends['bins']['mid'][...]
                chr_indices = hic.fends['bin_indices'][...]
            else:
                temp_mids = hic.fends['fends']['mid'][...]
                chr_indices = hic.fends['chr_indices'][...]

            # fill in chromosome list if empty. Otherwise check that all specified chromosomes exist.
            if not isinstance(chroms, list) or len(chroms) == 0:
                chroms = hic.fends['chromosomes'][...]
                valid = numpy.ones(chroms.shape[0], dtype=numpy.bool)
                for i in range(chroms.shape[0]):
                    if chr_indices[i + 1] - chr_indices[i] == 0:
                        valid[i] = False
                    elif hic.data['cis_indices'][chr_indices[i + 1]] - hic.data['cis_indices'][chr_indices[i]] == 0:
                        valid[i] = False
                chroms = chroms[valid]

            # Load raw counts
            bounds = numpy.zeros((len(chroms), 2), numpy.int64)
            for i, chrom in enumerate(chroms):
                chrint = hic.chr2int[chrom]
                bounds[i, 0] = hic.data['cis_indices'][chr_indices[chrint]]
                bounds[i, 1] = hic.data['cis_indices'][chr_indices[chrint + 1]]
            raw = numpy.zeros((numpy.sum(bounds[:, 1] - bounds[:, 0]), 3), dtype=numpy.int64)
            indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            mids = {}
            starts = numpy.zeros(len(chroms), dtype=numpy.int32)
            for i, chrom in enumerate(chroms):
                chrint = hic.chr2int[chrom]
                indices[i + 1] = indices[i] + bounds[i, 1] - bounds[i, 0]
                temp = hic.data['cis_data'][bounds[i, 0]:bounds[i, 1], :]
                temp[:, :2] -= chr_indices[chrint]
                raw[indices[i]:indices[i + 1], :] = temp
                mids[chrom] = temp_mids[chr_indices[chrint]:chr_indices[chrint + 1]]
                starts[i] = mids[chrom][0]

            # only consider coverage levels that are less than or equal to the number of cis reads
            coverages = coverages[numpy.where(numpy.sum(raw[:, 2]) >= coverages)]
            coverages.sort()
            if coverages[0] == 0:
                coverages[:-1] = coverages[1:]
                coverages[-1] = 0
            store_coverages = numpy.copy(coverages)
            total_reads = numpy.sum(raw[:, 2])
            if coverages.shape[0] > 0 and coverages[-1] == 0:
                coverages[-1] = total_reads
            coverages = coverages[::-1]
        else:
            coverages = None
            resolutions = None
        if self.comm is not None:
            coverages = self.comm.bcast(coverages, root=0)
            resolutions = self.comm.bcast(resolutions, root=0)
            chroms = self.comm.bcast(chroms, root=0)
        if coverages.shape[0] == 0:
            return None

        if self.rank == 0:
            # write arguements to h5dict
            if 'chromosomes' in self.storage:
                del self.storage['chromosomes']
            self.storage.create_dataset(name='chromosomes', data=numpy.array(chroms))
            if 'resolutions' in self.storage:
                del self.storage['resolutions']
            self.storage.create_dataset(name='resolutions', data=numpy.array(resolutions))
            if 'coverages' in self.storage:
                del self.storage['coverages']
            self.storage.create_dataset(name='coverages', data=numpy.array(store_coverages))
            if 'starts' in self.storage:
                del self.storage['starts']
            self.storage.create_dataset(name='starts', data=starts)
            self.storage.attrs['total_reads'] = total_reads

            # rebin data to highest resolution for faster processing
            remapped = {}
            new_mids = {}
            new_indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                start = (starts[i] / resolutions[0]) * resolutions[0]
                stop = ((mids[chrom][-1] - 1) / resolutions[0] + 1) * resolutions[0]
                N = (stop - start) / resolutions[0]
                mapping = (mids[chrom] - start) / resolutions[0]
                raw[indices[i]:indices[i + 1], 0] = mapping[raw[indices[i]:indices[i + 1], 0]]
                raw[indices[i]:indices[i + 1], 1] = mapping[raw[indices[i]:indices[i + 1], 1]]
                new_index = numpy.unique(raw[indices[i]:indices[i + 1], 0] * N + raw[indices[i]:indices[i + 1], 1])
                index = numpy.searchsorted(new_index, raw[indices[i]:indices[i + 1], 0] * N +
                                                      raw[indices[i]:indices[i + 1], 1])
                remapped[chrom] = numpy.zeros((new_index.shape[0], 3), dtype=numpy.int64)
                remapped[chrom][:, 0] = new_index / N
                remapped[chrom][:, 1] = new_index % N
                remapped[chrom][:, 2] = numpy.bincount(index, weights=raw[indices[i]:indices[i + 1], 2])
                new_indices[i + 1] = new_index.shape[0] + new_indices[i]
                new_mids[chrom] = (start + resolutions[0] / 2 + numpy.arange(N) *
                                   resolutions[0]).astype(numpy.int32)
            indices = new_indices.astype(numpy.int64)
            mids = new_mids
            raw = numpy.zeros((indices[-1], 3), dtype=numpy.int64)
            for i, chrom in enumerate(chroms):
                raw[indices[i]:indices[i + 1], :] = remapped[chrom]
            del remapped

        # Transfer mids data
        chrom_ranges = numpy.round(numpy.linspace(0, len(chroms), self.num_procs + 1)).astype(numpy.int32)
        if self.rank == 0:
            for i in range(1, self.num_procs):
                for j in range(chrom_ranges[i], chrom_ranges[i + 1]):
                    chrom = chroms[j]
                    self.comm.send(mids[chrom].shape[0], dest=i, tag=3)
                    self.comm.Send(mids[chrom], dest=i, tag=4)
        else:
            mids = {}
            for i in range(chrom_ranges[self.rank], chrom_ranges[self.rank + 1]):
                chrom = chroms[i]
                N = self.comm.recv(source=0, tag=3)
                mids[chrom] = numpy.zeros(N, dtype=numpy.int32)
                self.comm.Recv(mids[chrom], source=0, tag=4)

        # cycle through coverages
        for c, cov in enumerate(coverages):
            raw_counts = {}
            if self.rank == 0:
                if not self.silent:
                    print >> sys.stderr, ("\r%s\rDownsampling to %i coverage") % (' ' * 120, cov),
                raw, indices = self._downsample(raw, indices, cov, RNG)
                for i in range(chrom_ranges[1]):
                    chrom = chroms[i]
                    raw_counts[chrom] = raw[indices[i]:indices[i + 1], :]
                for i in range(1, self.num_procs):
                    for j in range(chrom_ranges[i], chrom_ranges[i + 1]):
                        chrom = chroms[j]
                        self.comm.send(indices[j + 1] - indices[j], dest=i, tag=1)
                        self.comm.Send(raw[indices[j]:indices[j + 1], :], dest=i, tag=2)
            else:
                for i in range(chrom_ranges[self.rank], chrom_ranges[self.rank + 1]):
                    chrom = chroms[i]
                    N = self.comm.recv(source=0, tag=1)
                    raw_counts[chrom] = numpy.zeros((N, 3), dtype=numpy.int64)
                    self.comm.Recv(raw_counts[chrom], source=0, tag=2)

            # cycle through resolutions
            for r, res in enumerate(resolutions):
                dists = {}
                norms = {}
                valids = {}
                corrs = {}

                # For each chromosome, normalize and find distance-corrected matrix
                if self.rank == 0:
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rCoverage %i Resolution %i - Normalizing counts") % (
                            ' ' * 120, cov, res),
                for chrom in raw_counts:
                    norm, dist, valid_rows = self._normalize(chrom, raw_counts[chrom], mids[chrom], res)
                    dists[chrom] = dist
                    norms[chrom] = norm
                    valids[chrom] = valid_rows
                if self.rank == 0:
                    for i in range(1, self.num_procs):
                        self._transfer_dict(dists, 0, i)
                        self._transfer_dict(norms, 0, i)
                        self._transfer_dict(valids, 0, i)
                else:
                    self._transfer_dict(dists, 0, self.rank)
                    self._transfer_dict(norms, 0, self.rank)
                    self._transfer_dict(valids, 0, self.rank)
                    dists = None

                # cycle through chromosomes finding correlation matrices
                for i, chrom in enumerate(chroms):
                    task = True
                    if self.rank == 0:
                        if cov == total_reads:
                            key = '%s.0C.%iR' % (chrom, res)
                        else:
                            key = '%s.%iC.%iR' % (chrom, cov, res)
                        if 'corr.%s' % key in self.storage:
                            task = False
                        if self.comm is not None:
                            task = self.comm.bcast(task, root=0)
                        if task:
                            if not self.silent:
                                print >> sys.stderr, ("\r%s\rCoverage %i Resolution %i - Correlating chrom %s") % (
                                    ' ' * 120, cov, res, chrom),
                            corrs[chrom] = self._find_correlations(norms[chrom], valids[chrom])
                    else:
                        task = self.comm.bcast(task, root=0)
                        if task:
                            self._find_correlations()

                # write resulting matrices to hdf5 file
                if self.rank == 0:
                    if not self.silent:
                        print >> sys.stderr, ("\r%s\rCoverage %i Resolution %i - Writing results") % (
                            ' ' * 120, cov, res),
                    for chrom in chroms:
                        if cov == total_reads:
                            key = '%s.0C.%iR' % (chrom, res)
                        else:
                            key = '%s.%iC.%iR' % (chrom, cov, res)
                        if valids[chrom] is None:
                            self.storage.attrs['%s.invalid' % (key)] = True
                        elif chrom in corrs:
                            self.storage.create_dataset(name="valid.%s" % (key),
                                data=valids[chrom])
                            self.storage.create_dataset(name="dist.%s" % (key), data=dists[chrom])
                            self.storage.create_dataset(name="corr.%s" % (key), data=corrs[chrom])
        if self.rank == 0 and not self.silent:
            print >> sys.stderr, ("\r%s\r") % (' ' * 120),
        return None

    def find_quality_scores(self, force=False):
        """
        Find QuASAR quality scores across whole dataset.

        :param force: Recalculate scores even if they have already been stored in the QuASAR object.
        :type force: bool.

        :returns: A structured numpy array with the fields resolution', 'coverage', and 'score'.
        """
        if self.rank > 0:
            return None
        if 'chromosomes' not in self.storage:
            return None
        chroms = self.storage['chromosomes'][...]
        coverages = self.storage['coverages'][...]
        tcoverages = numpy.copy(coverages)
        if tcoverages[-1] == 0:
            tcoverages[-1] = self.storage.attrs['total_reads']
        resolutions = self.storage['resolutions'][...]
        if ('quality_results' in self.storage and not force and
            self.storage['quality_results'].shape == (chroms.shape[0] + 1, resolutions.shape[0], coverages.shape[0])):
            scores = self.storage['quality_results'][...]
            results = numpy.zeros(scores.shape[1] * scores.shape[2], dtype=numpy.dtype([
                ('resolution', numpy.int64), ('coverage', numpy.int64), ('score', numpy.float64)]))
            results['score'] = scores[-1, :, :].ravel(order='C')
            results['resolution'] = numpy.repeat(resolutions, coverages.shape[0])
            results['coverage'] = numpy.tile(tcoverages, resolutions.shape[0])
            return results
        elif 'quality_results' in self.storage:
            del self.storage['quality_results']
        scores = numpy.zeros((chroms.shape[0] + 1, resolutions.shape[0], coverages.shape[0]), dtype=numpy.float64)
        scores.fill(numpy.nan)
        for j, res in enumerate(resolutions):
            for k, cov in enumerate(coverages):
                temp = numpy.zeros(4, dtype=numpy.float64)
                for i, chrom in enumerate(chroms):
                    key = '%s.%iC.%iR' % (chrom, cov, res)
                    if 'valid.%s' % key in self.storage:
                        valid_rows = self.storage['valid.%s' % key][...]
                        dists = self.storage['dist.%s' % key][...]
                        corrs = self.storage['corr.%s' % key][...]
                        valid = numpy.zeros(corrs.shape, dtype=numpy.bool)
                        N, M = corrs.shape
                        valid = numpy.zeros((N, M), dtype=numpy.int32)
                        for l in range(min(N - 1, M)):
                            P = N - l - 1
                            valid[:P, l] = valid_rows[(l + 1):] * valid_rows[:P]
                        valid[numpy.where((numpy.abs(dists) == numpy.inf) | (numpy.abs(corrs) == numpy.inf))] = 0
                        where = numpy.where(valid)
                        N = where[0].shape[0]
                        if N > 0:
                            trans = numpy.sum(corrs[where] * dists[where])
                            corrs = numpy.sum(corrs[where])
                            dists = numpy.sum(dists[where])
                            scores[i, j, k] = trans / dists - corrs / N
                            temp += [trans, dists, corrs, N]
                scores[-1, j, k] = temp[0] / temp[1] - temp[2] / temp[3]
        self.storage.create_dataset(name='quality_results', data=scores)
        results = numpy.zeros(scores.shape[1] * scores.shape[2], dtype=numpy.dtype([
            ('resolution', numpy.int64), ('coverage', numpy.int64), ('score', numpy.float64)]))
        results['score'] = scores[-1, :, :].ravel(order='C')
        results['resolution'] = numpy.repeat(resolutions, coverages.shape[0])
        results['coverage'] = numpy.tile(coverages, resolutions.shape[0])
        return results

    def find_replicate_scores(self, replicate, force=False):
        """
        Find QuASAR replicate scores across whole dataset.

        :param replicate: A class:`Quasar` object to calculate replicate scores with. If this function has been previously called with a different sample, the current transformed matrices will be deleted prior to calculating new matrices.
        :type replicate: class:`Quasar` class object.
        :param force: Recalculate scores even if they have already been stored in the QuASAR object.
        :type force: bool.

        :returns: A structured numpy array with the fields 'resolution', 'coverage', and 'score'.
        """
        if self.rank > 0:
            return None
        if 'chromosomes' not in self.storage or 'chromosomes' not in replicate.storage:
            return None
        chroms = numpy.intersect1d(self.storage['chromosomes'][...], replicate.storage['chromosomes'][...])
        resolutions = numpy.intersect1d(self.storage['resolutions'][...], replicate.storage['resolutions'][...])
        coverages = numpy.intersect1d(self.storage['coverages'][...], replicate.storage['coverages'][...])
        if coverages[0] == 0:
            coverages[:-1] = coverages[1:]
            coverages[-1] = 0
        tcoverages = numpy.copy(coverages)
        if tcoverages[-1] == 0:
            tcoverages[-1] = (self.storage.attrs['total_reads'] + replicate.storage.attrs['total_reads']) / 2
        if 'replicate_results' in self.storage and 'replicate_results' in replicate.storage and not force:
            if (numpy.array_equal(self.storage['replicate_results'][...],
                replicate.storage['replicate_results'][...]) and
                self.storage['replicate_results'].shape ==
                (chroms.shape[0] + 1, resolutions.shape[0], coverages.shape[0])):
                scores = self.storage['replicate_results'][...]
                results = numpy.zeros(scores.shape[1] * scores.shape[2], dtype=numpy.dtype([
                    ('resolution', numpy.int64), ('coverage', numpy.int64), ('score', numpy.float64)]))
                results['score'] = scores[-1, :, :].ravel(order='C')
                results['resolution'] = numpy.repeat(resolutions, coverages.shape[0])
                results['coverage'] = numpy.tile(tcoverages, resolutions.shape[0])
                return results
            else:
                if 'replicate_results' in self.storage:
                    del self.storage['replicate_results']
                if 'replicate_results' in replicate.storage:
                    del replicate.storage['replicate_results']
        else:
            if 'replicate_results' in self.storage:
                del self.storage['replicate_results']
            if 'replicate_chromosomes' in self.storage:
                del self.storage['replicate_chromosomes']
            if 'replicate_resolutions' in self.storage:
                del self.storage['replicate_resolutions']
            if 'replicate_coverages' in self.storage:
                del self.storage['replicate_coverages']
            if 'replicate_results' in replicate.storage:
                del replicate.storage['replicate_results']
            if 'replicate_chromosomes' in replicate.storage:
                del replicate.storage['replicate_chromosomes']
            if 'replicate_resolutions' in replicate.storage:
                del replicate.storage['replicate_resolutions']
            if 'replicate_coverages' in replicate.storage:
                del replicate.storage['replicate_coverages']
        starts1 = numpy.zeros(chroms.shape[0], dtype=numpy.int64)
        starts2 = numpy.zeros(chroms.shape[0], dtype=numpy.int64)
        for i, chrom in enumerate(chroms):
            starts1[i] = numpy.where(self.storage['chromosomes'][...] == chrom)[0][0]
            starts2[i] = numpy.where(replicate.storage['chromosomes'][...] == chrom)[0][0]
        scores = numpy.zeros((chroms.shape[0] + 1, resolutions.shape[0], coverages.shape[0]), dtype=numpy.float64)
        scores.fill(numpy.nan)
        for j, res in enumerate(resolutions):
            for k, cov in enumerate(coverages):
                temp = numpy.zeros(6, dtype=numpy.float64)
                for i, chrom in enumerate(chroms):
                    key = '%s.%iC.%iR' % (chrom, cov, res)
                    if 'valid.%s' % key in self.storage and 'valid.%s' % key in replicate.storage:
                        valid_rows1 = self.storage['valid.%s' % key][...]
                        dists1 = self.storage['dist.%s' % key][...]
                        corrs1 = self.storage['corr.%s' % key][...]
                        valid1 = numpy.zeros(corrs1.shape, dtype=numpy.bool)
                        N, M = corrs1.shape
                        for l in range(min(N - 1, M)):
                            P = N - l - 1
                            valid1[:P, l] = valid_rows1[(l + 1):] * valid_rows1[:P]
                        valid1[numpy.where((numpy.abs(dists1) == numpy.inf) | (numpy.abs(corrs1) == numpy.inf))] = 0
                        valid_rows2 = replicate.storage['valid.%s' % key][...]
                        dists2 = replicate.storage['dist.%s' % key][...]
                        corrs2 = replicate.storage['corr.%s' % key][...]
                        valid2 = numpy.zeros(corrs2.shape, dtype=numpy.bool)
                        N, M = corrs2.shape
                        for l in range(min(N - 1, M)):
                            P = N - l - 1
                            valid2[:P, l] = valid_rows2[(l + 1):] * valid_rows2[:P]
                        valid2[numpy.where((numpy.abs(dists2) == numpy.inf) | (numpy.abs(corrs2) == numpy.inf))] = 0
                        start1 = starts1[i]
                        start2 = starts2[i]
                        if start2 > start1:
                            start1 = (start2 - start1) / res
                            start2 = 0
                        elif start2 > start1:
                            start2 = (start1 - start2) / res
                            start1 = 0
                        else:
                            start1 = 0
                            start2 = 0
                        stop1 = corrs1.shape[0] - start1
                        stop2 = corrs2.shape[0] - start2
                        if stop2 > stop1:
                            stop2 = stop1 + start2
                            stop1 = start1 + stop1
                        else:
                            stop1 = stop2 + start1
                            stop2 = start2 + stop2
                        stop3 = min(corrs1.shape[1], corrs2.shape[1])
                        valid1 = valid1[start1:stop1, :stop3]
                        valid2 = valid2[start2:stop2, :stop3]
                        valid = numpy.where(valid1 & valid2)
                        if valid[0].shape[0] == 0:
                            continue
                        trans1 = corrs1[start1:stop1, :stop3][valid] * dists1[start1:stop1, :stop3][valid]
                        trans2 = corrs2[start2:stop2, :stop3][valid] * dists2[start2:stop2, :stop3][valid]
                        X = numpy.sum(trans1)
                        Y = numpy.sum(trans2)
                        X2 = numpy.sum(trans1 ** 2.0)
                        Y2 = numpy.sum(trans2 ** 2.0)
                        XY = numpy.sum(trans1 * trans2)
                        N = valid[0].shape[0]
                        if N == 0:
                            continue
                        temp += [X, Y, X2, Y2, XY, N]
                        Xmu = X / N
                        Ymu = Y / N
                        X2mu = X2 / N
                        Y2mu = Y2 / N
                        XYmu = XY / N
                        if Xmu ** 2.0 > X2mu or Ymu ** 2.0 > Y2mu:
                            continue
                        Xstd = (X2mu - Xmu ** 2.0) ** 0.5
                        Ystd = (Y2mu - Ymu ** 2.0) ** 0.5
                        if Xstd == 0 or Ystd == 0:
                            continue
                        scores[i, j, k] = (XYmu - Xmu * Ymu) / (Xstd * Ystd)
                Xmu = temp[0] / temp[5]
                Ymu = temp[1] / temp[5]
                X2mu = temp[2] / temp[5]
                Y2mu = temp[3] / temp[5]
                XYmu = temp[4] / temp[5]
                Xstd = (X2mu - Xmu ** 2.0) ** 0.5
                Ystd = (Y2mu - Ymu ** 2.0) ** 0.5
                scores[-1, j, k] = (XYmu - Xmu * Ymu) / (Xstd * Ystd)
        self.storage.create_dataset(name='replicate_results', data=scores)
        self.storage.create_dataset(name='replicate_chromosomes', data=chroms)
        self.storage.create_dataset(name='replicate_resolutions', data=resolutions)
        self.storage.create_dataset(name='replicate_coverages', data=tcoverages)
        if 'replicate_results' not in replicate.storage:
            replicate.storage.create_dataset(name='replicate_results', data=scores)
            replicate.storage.create_dataset(name='replicate_chromosomes', data=chroms)
            replicate.storage.create_dataset(name='replicate_resolutions', data=resolutions)
            replicate.storage.create_dataset(name='replicate_coverages', data=tcoverages)
        results = numpy.zeros(scores.shape[1] * scores.shape[2], dtype=numpy.dtype([
            ('resolution', numpy.int64), ('coverage', numpy.int64), ('score', numpy.float64)]))
        results['score'] = scores[-1, :, :].ravel(order='C')
        results['resolution'] = numpy.repeat(resolutions, coverages.shape[0])
        results['coverage'] = numpy.tile(tcoverages, resolutions.shape[0])
        return results

    def print_report(self, filename):
        """
        Write QuASAR scores to output file.

        :param filename: The location to write the report to. The suffix will be used to determine the output format.
        :type filename: str.

        :returns: None
        """
        if self.rank > 0:
            return None
        qscores = None
        rscores = None
        if 'quality_results' in self.storage:
            qscores = self.storage['quality_results'][...]
        if 'replicate_results' in self.storage:
            rscores = self.storage['replicate_results'][...]
        if qscores is None and rscores is None:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rNo scores found. Calculate quality and/or replicate scores first.\n") % (
                    ' ' * 80),
            return None
        format = filename.split('.')[-1]
        if format not in ['pdf', 'txt']:
            if not self.silent:
                print >> sys.stderr, ("\r%s\rThe output format was not recognized.\n") % (
                    ' ' * 80),
            return None
        if format == 'txt':
            self._print_txt_report(filename, qscores, rscores)
        elif format == 'pdf':
            if 'pyx' in sys.modules.keys():
                self._print_pdf_report(filename, qscores, rscores)
            elif not self.silent:
                print >> sys.stderr, ("\r%s\rThe pyx package is needed for writing PDFs.\n") % (' ' * 80),
        else:
            self._print_html_report(filename, qscores, rscores)
        return None

    def _transfer_dict(self, data, dest, source):
        if self.rank == dest:
            key = self.comm.recv(source=source, tag=5)
            while key:
                shape, dtype = self.comm.recv(source=source, tag=6)
                if shape is None:
                    data[key] = None
                else:
                    data[key] = numpy.zeros(shape, dtype=dtype)
                    self.comm.Recv(data[key], source=source, tag=7)
                key = self.comm.recv(source=source, tag=5)
        elif self.rank == source:
            for key, value in data.iteritems():
                self.comm.send(key, dest=dest, tag=5)
                if value is None:
                    self.comm.send([None, None], dest=dest, tag=6)
                else:
                    self.comm.send([value.shape, value.dtype], dest=dest, tag=6)
                    self.comm.Send(value, dest=dest, tag=7)
            self.comm.send(False, dest=dest, tag=5)
        return None

    def _downsample(self, data, indices, target_count, rng=None):
        if target_count == 0:
            return numpy.copy(data), numpy.copy(indices)
        elif numpy.sum(data[:, 2]) == target_count:
            return numpy.copy(data), numpy.copy(indices)
        if rng is None:
            rng = numpy.random.RandomState()
        initial_count = numpy.sum(data[:, 2])
        percent = target_count / float(initial_count)

        # select which reads to keep, based on percent of reads to keep
        keep = rng.rand(initial_count) < percent

        # adjust mismatch between selected read count and target read count by selecting reads to add/remove
        kept = numpy.sum(keep)
        if kept > target_count:
            pool_size = kept
            adjust_size = kept - target_count
            remove = True
        elif kept < target_count:
            pool_size = initial_count - kept
            adjust_size = target_count - kept
            remove = False
        else:
            adjust_size = 0
        reads = {}
        while len(reads) < adjust_size:
            temp_rand = rng.randint(0, pool_size, adjust_size * 2)
            for i in range(temp_rand.shape[0]):
                reads[temp_rand[i]] = None
                if len(reads) == adjust_size:
                    break
        if adjust_size > 0:
            if remove:
                where = numpy.where(keep)[0]
                for i in where[reads.keys()]:
                    keep[i] = False
            else:
                where = numpy.where(numpy.logical_not(keep))[0]
                for i in where[reads.keys()]:
                    keep[i] = True

        # adjust read counts in data
        counts = numpy.repeat(numpy.arange(data.shape[0]), data[:, 2])
        new_data = numpy.copy(data)
        new_data[:, 2] = numpy.bincount(counts, weights=keep, minlength=data.shape[0])
        new_indices = numpy.zeros(indices.shape[0], dtype=numpy.int64)
        for i in range(1, new_indices.shape[0]):
            new_indices[i] = numpy.sum(new_data[indices[i - 1]:indices[i], 2] > 0) + new_indices[i - 1]
        new_data = new_data[numpy.where(new_data[:, 2] > 0)[0], :]
        return new_data, new_indices

    def _normalize(self, chrom, raw, mids, binsize):
        width = 100
        # convert into square matrix
        start = (mids[0] / binsize) * binsize
        stop = ((mids[-1] - 1) / binsize + 1) * binsize
        N = ((stop - start) / binsize).astype(numpy.int64)
        mapping = (mids - start) / binsize
        data = numpy.bincount(mapping[raw[:, 0]] * N + mapping[raw[:, 1]], minlength=(N * N),
                              weights=raw[:, 2]).reshape(N, N)
        indices = numpy.triu_indices(N, 0)
        data[indices[1], indices[0]] = data[indices]
        data[numpy.arange(N), numpy.arange(N)] *= 2
        data = data.astype(numpy.float64)

        # filter rows with too few observations
        row_counts = numpy.bincount(mapping, minlength=N)
        prev_valid_rows = N + 1
        # we need several non-zero observations for finding correlations
        sums = numpy.zeros(N, dtype=numpy.int32)
        for i in range(min(width, N - 1)):
            temp1 = numpy.arange(N - i - 1)
            temp2 = numpy.arange(i + 1, N)
            sums[:(N - i - 1)] += data[temp1, temp2] > 0
            sums[(i + 1):N] += data[temp1, temp2] > 0
        valid_rows = sums >= 1
        while prev_valid_rows > numpy.sum(valid_rows):
            data[numpy.logical_not(valid_rows), :] = 0
            data[:, numpy.logical_not(valid_rows)] = 0
            prev_valid_rows = numpy.sum(valid_rows)
            sums = numpy.zeros(N, dtype=numpy.int32)
            for i in range(min(width, N - 1)):
                temp1 = numpy.arange(N - i - 1)
                temp2 = numpy.arange(i + 1, N)
                sums[:(N - i - 1)] += data[temp1, temp2] > 0
                sums[(i + 1):N] += data[temp1, temp2] > 0
            valid_rows = sums >= 1
        valid = numpy.where(valid_rows)[0]
        if valid.shape[0] == 0:
            return None, None, None

        # find distance dependent signal
        bg = numpy.ones(N, dtype=numpy.float64)
        for i in range(N):
            temp1 = numpy.arange(N - i)
            temp2 = numpy.arange(i, N)
            obs = numpy.sum(data[temp1, temp2])
            pos = float(numpy.sum(valid_rows[temp1] & valid_rows[temp2]))
            if obs > 0:
                bg[i] = obs / pos

        # find compact data array for weighting correlation matrix
        dist = numpy.zeros((N, min(width, N - 1)), dtype=numpy.float64)
        for i in range(dist.shape[1]):
            M = N - i - 1
            temp1 = numpy.arange(M)
            temp2 = numpy.arange(i + 1, N)
            dist[:M, i] = (data[temp1, temp2] + 1) ** 0.5

        # find fill data array for finding correlations
        norm = numpy.zeros((N, N), dtype=numpy.float64)
        for i in valid:
            norm[i, valid] = data[i, valid]
        for i in range(N):
            temp1 = numpy.arange(N - i)
            temp2 = numpy.arange(i, N)
            norm[temp1, temp2] /= bg[i]
            norm[temp2, temp1] = norm[temp1, temp2]
        return norm, dist, valid_rows

    def _find_correlations(self, norm=None, vrows=None):
        width = 100
        if self.rank == 0:
            if norm is not None:
                N = norm.shape[0]
            else:
                N = None
        else:
            N = None
        if self.comm is not None:
            N = self.comm.bcast(N, root=0)
        if N is None:
            return None
        if self.rank != 0:
            norm = numpy.zeros((N, N), dtype=numpy.float64)
            vrows = numpy.zeros(N, dtype=numpy.bool)
        if self.comm is not None:
            self.comm.Bcast(norm, root=0)
            self.comm.Bcast(vrows, root=0)
        valid = numpy.where(vrows)[0]
        node_ranges = numpy.round(numpy.linspace(0, N, self.num_procs + 1)).astype(numpy.int32)
        M = node_ranges[self.rank + 1] - node_ranges[self.rank]
        corr = numpy.zeros((M, min(N - 1, width)), dtype=numpy.float64)
        corr.fill(numpy.inf)
        _quasar.find_correlations(
            norm,
            vrows.astype(numpy.int32),
            corr,
            node_ranges[self.rank],
            node_ranges[self.rank + 1])
        if self.rank == 0:
            corrs = numpy.zeros((N, min(N - 1, width)), dtype=numpy.float64)
            corrs[:node_ranges[1], :] = corr
            for j in range(1, self.num_procs):
                self.comm.Recv(corrs[node_ranges[j]:node_ranges[j + 1], :], source=j, tag=9)
            corrs[numpy.where(numpy.isnan(corrs))] = numpy.inf
            return corrs
        else:
            self.comm.Send(corr, dest=0, tag=9)
        return None

    def _print_txt_report(self, filename, qscores, rscores):
        output = open(filename, 'w')
        if qscores is not None:
            chromosomes = self.storage['chromosomes'][...]
            resolutions = self.storage['resolutions'][...]
            coverages = self.storage['coverages'][...]
            if numpy.where(coverages == 0)[0].shape[0] > 0:
                coverages[numpy.where(coverages == 0)[0][0]] = self.storage.attrs['total_reads']
            print >> output, 'Quality Score Results\n'
            temp = ['Resolution', 'Coverage', 'All'] + list(chromosomes)
            print >> output, '\t'.join(temp)
            for i, res in enumerate(resolutions):
                if res < 1000:
                    label = "%i bp" % res
                elif res < 1000000:
                    label = "%i Kb" % (res / 1000)
                else:
                    label = "%i Mb" % (res / 1000000)
                for j, cov in enumerate(coverages):
                    temp = [label, self._num2str(cov), "%0.6f" % qscores[-1, i, j]]
                    for k in range(chromosomes.shape[0]):
                        temp.append("%0.6f" % qscores[k, i, j])
                    print >> output, '\t'.join(temp)
            print >> output, '\n'

            # if there are sufficient data points, calculate estimated maximum
            if coverages.shape[0] > 3:

                def f(x, x0, k, L):
                    return L / (1 + numpy.exp(-k * (x - x0)))

                Xs = numpy.log10(coverages)
                print >> output, "Estimated quality-coverage curves (Maximum Quality / (1 + e^(-Scale * (x - Inflection)))"
                for i, res in enumerate(resolutions):
                    Ys = qscores[-1, i, :]
                    if res < 1000:
                        label = "%i bp" % res
                    elif res < 1000000:
                        label = "%i Kb" % (res / 1000)
                    else:
                        label = "%i Mb" % (res / 1000000)
                    try:
                        params = curve_fit(f, Xs, Ys, p0=(Xs[-1], 0.5, Ys[-1] * 2), maxfev=(5000*Xs.shape[0]),
                             bounds=((-numpy.inf, -numpy.inf, 0), (numpy.inf, numpy.inf, 2)))[0]
                        Y1s = f(Xs, *params)
                        print >> output, "Resolution: %s, Maximum Quality: %f, Scale: %f, Inflection: %f, Mean error: %e" % (label, params[2], params[1], params[0], numpy.mean((Ys - Y1s) ** 2.0))
                    except:
                        print >> output, "Resolution: %s, curve could not be estimated" % label
                print >> output, '\n'

            # if there are multiple coverages, estimate maximum usable resolution
            if resolutions.shape[0] > 1:
                Xs = numpy.log10(resolutions)
                Ys = qscores[-1, :, -1]
                pos = 0
                while pos < Xs.shape[0] - 1 and self.strict_qcutoff > Ys[pos + 1]:
                    pos += 1
                if self.strict_qcutoff < Ys[0]:
                    print >> output, "Strict quality maximum resolution is above the resolutions tested."
                elif self.strict_qcutoff > numpy.amax(Ys):
                    print >> output, "Strict quality maximum resolution is below the resolutions tested."
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (Xs[pos + 1] - Xs[pos])
                    B = Ys[pos] - Xs[pos] * A
                    print >> output, "Strict quality maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.strict_qcutoff - B) / A))))
                pos = 0
                while pos < Xs.shape[0] - 1 and self.loose_qcutoff > Ys[pos + 1]:
                    pos += 1
                if self.loose_qcutoff < Ys[0]:
                    print >> output, "Loose quality maximum resolution is above the resolutions tested."
                elif self.loose_qcutoff > numpy.amax(Ys):
                    print >> output, "Loose quality maximum resolution is below the resolutions tested."
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (Xs[pos + 1] - Xs[pos])
                    B = Ys[pos] - Xs[pos] * A
                    print >> output, "Loose quality maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.loose_qcutoff - B) / A))))
                print >> output, '\n'

        if rscores is not None:
            chromosomes = self.storage['replicate_chromosomes'][...]
            resolutions = self.storage['replicate_resolutions'][...]
            coverages = self.storage['replicate_coverages'][...]
            print >> output, 'Replicate Score Results\n'
            temp = ['Resolution', 'Coverage', 'All'] + list(chromosomes)
            print >> output, '\t'.join(temp)
            for i, res in enumerate(resolutions):
                if res < 1000:
                    label = "%i bp" % res
                elif res < 1000000:
                    label = "%i Kb" % (res / 1000)
                else:
                    label = "%i Mb" % (res / 1000000)
                for j, cov in enumerate(coverages):
                    temp = [label, self._num2str(cov), "%0.6f" % rscores[-1, i, j]]
                    for k in range(chromosomes.shape[0]):
                        temp.append("%0.6f" % rscores[k, i, j])
                    print >> output, '\t'.join(temp)
            print >> output, '\n'

            # if there are multiple coverages, estimate maximum usable resolution
            if resolutions.shape[0] > 1:
                Xs = numpy.log10(resolutions)
                Ys = rscores[-1, :, -1]
                pos = 0
                while pos < Xs.shape[0] - 1 and self.strict_rcutoff > Ys[pos + 1]:
                    pos += 1
                if self.strict_rcutoff < Ys[0]:
                    print >> output, "Strict replicate maximum resolution is above the resolutions tested."
                elif self.strict_rcutoff > numpy.amax(Ys):
                    print >> output, "Strict replicate maximum resolution is below the resolutions tested."
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (Xs[pos + 1] - Xs[pos])
                    B = Ys[pos] - Xs[pos] * A
                    print >> output, "Strict replicate maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.strict_rcutoff - B) / A))))
                pos = 0
                while pos < Xs.shape[0] - 1 and self.loose_rcutoff > Ys[pos + 1]:
                    pos += 1
                if self.loose_rcutoff < Ys[0]:
                    print >> output, "Loose replicate maximum resolution is above the resolutions tested."
                elif self.loose_rcutoff > numpy.amax(Ys):
                    print >> output, "Loose replicate maximum resolution is below the resolutions tested."
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (Xs[pos + 1] - Xs[pos])
                    B = Ys[pos] - Xs[pos] * A
                    print >> output, "Loose replicate maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.loose_rcutoff - B) / A))))
                print >> output, '\n'
        output.close()
        return None

    def _print_pdf_report(self, filename, qscores, rscores):
        c = canvas.canvas()
        H = 0
        if qscores is not None:
            resolutions = self.storage['resolutions'][...]
            coverages = self.storage['coverages'][...]
            chromosomes = self.storage['chromosomes'][...]
            if numpy.where(coverages == 0)[0].shape[0] > 0:
                coverages[numpy.where(coverages == 0)[0][0]] = self.storage.attrs['total_reads']
            c.text(H, 0, "Quality Score Results", [text.halign.left, text.valign.bottom, text.size(0)])
            H -= 0.6
            hoffset = 1.7
            for i in range(coverages.shape[0] + 1):
                c.stroke(path.line(hoffset + i * 2.0, H + 0.3, hoffset + i * 2.0, H - 0.3 * resolutions.shape[0]))
                if i < coverages.shape[0]:
                    c.text(hoffset + (i + 0.5) * 2.0, H + 0.05, self._num2str(coverages[i]),
                        [text.halign.center, text.valign.bottom, text.size(-2)])
            for i in range(resolutions.shape[0] + 1):
                c.stroke(path.line(0, H - 0.3 * i, hoffset + coverages.shape[0] * 2.0, H - 0.3 * i))
            c.text(hoffset + 1.0 * coverages.shape[0], H + 0.35, "Coverage",
                [text.halign.center, text.valign.bottom, text.size(-2)])
            c.text(hoffset - 0.1, H + 0.05, "Resolution",
                [text.halign.right, text.valign.bottom, text.size(-2)])
            H -= 0.3
            for i, res in enumerate(resolutions):
                if res < 1000:
                    label = "%i bp" % res
                elif res < 1000000:
                    label = "%i Kb" % (res / 1000)
                else:
                    label = "%i Mb" % (res / 1000000)
                c.text(hoffset - 0.1, H + 0.05, label, [text.halign.right, text.valign.bottom, text.size(-2)])
                for j, cov in enumerate(coverages):
                    c.text(hoffset + (j + 0.5) * 2.0, H + 0.05, '%0.6f' % qscores[-1, i, j],
                        [text.halign.center, text.valign.bottom, text.size(-2)])
                H -= 0.3
            H -= 0.6
            hoffset = 0.9

            # if there are sufficient data points, calculate estimated maximum quality
            if coverages.shape[0] > 3:
                width = 17.0 / min(3, resolutions.shape[0]) - 0.3

                def f(x, x0, k, L):
                    return L / (1 + numpy.exp(-k * (x - x0)))

                Xs = coverages
                lXs = numpy.log10(coverages)
                c.text(0, H, r"Estimated quality-coverage curves $\frac{Maximum Quality}{1 + e^{-Scale * (x - Inflection)}}$", [text.halign.left, text.valign.bottom, text.size(-2)])
                H -= 0.3
                passed = False
                for i, res in enumerate(resolutions):
                    Ys = qscores[-1, i, :]
                    c1 = self._plot_graph(Xs, Ys, width, 'q')
                    if res < 1000:
                        label = "%i bp" % res
                    elif res < 1000000:
                        label = "%i Kb" % (res / 1000)
                    else:
                        label = "%i Mb" % (res / 1000000)
                    try:
                        params = curve_fit(f, lXs, Ys, p0=(lXs[-1], 0.5, Ys[-1] * 2), maxfev=(5000*Xs.shape[0]),
                             bounds=((-numpy.inf, -numpy.inf, 0), (numpy.inf, numpy.inf, 2)))[0]
                        Y1s = f(lXs, *params)
                        passed = True
                        minX = numpy.log2(numpy.amin(Xs))
                        maxX = numpy.log2(numpy.amax(Xs))
                        minY = numpy.amin(Ys)
                        maxY = numpy.amax(Ys)
                        spanY = maxY - minY
                        minY -= 0.05 * spanY
                        maxY += 0.05 * spanY
                        c3 = canvas.canvas()
                        c3.insert(self._plot_line(Xs, Y1s, width, minX, maxX, minY, maxY, color.gray(0.5)))
                        c3.insert(c1)
                        c1 = c3
                        c1.text(hoffset, -0.15, "Resolution: %s" % label,
                            [text.halign.left, text.valign.top, text.size(-2)])
                        c1.text(hoffset, -0.45, "Maximum Quality: %f" % params[2],
                            [text.halign.left, text.valign.top, text.size(-2)])
                        c1.text(hoffset, -0.75, "Scale: %f," % params[1],
                            [text.halign.left, text.valign.top, text.size(-2)])
                        c1.text(hoffset, -1.05, "Inflection: %f" % params[0],
                            [text.halign.left, text.valign.top, text.size(-2)])
                        c1.text(hoffset, -1.35, "Mean Error: %e" % numpy.mean((Ys - Y1s) ** 2.0),
                            [text.halign.left, text.valign.top, text.size(-2)])
                    except:
                        c1.text(hoffset, -0.15, "Resolution: %s" % label,
                            [text.halign.left, text.valign.top, text.size(-2)])
                        c1.text(hoffset, -0.45, "curve could not be estimated",
                            [text.halign.left, text.valign.top, text.size(-2)])
                    c.insert(c1, [trafo.translate((width + 0.3) * i, H - width)])
                if passed:
                    H -= width + 2.0
                else:
                    H -= width + 1.1

            # if there are multiple coverages, estimate maximum usable resolution
            if resolutions.shape[0] > 1:
                Xs = resolutions
                Ys = qscores[-1, :, -1]
                width = 17.0 / 3.0 - 0.3
                c1 = self._plot_graph(Xs, Ys, width, 'q', 'res')
                pos = 0
                while pos < Xs.shape[0] - 1 and self.strict_qcutoff > Ys[pos + 1]:
                    pos += 1
                if self.strict_qcutoff < Ys[0]:
                    c1.text(0, -0.15, "Strict quality maximum resolution is above the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                elif self.strict_qcutoff > numpy.amax(Ys):
                    c1.text(0, -0.15, "Strict quality maximum resolution is below the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (numpy.log10(Xs[pos + 1]) - numpy.log10(Xs[pos]))
                    B = Ys[pos] - numpy.log10(Xs[pos]) * A
                    c1.text(0, -0.15, "Strict quality maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.strict_qcutoff - B) / A)))),
                        [text.halign.left, text.valign.top, text.size(-2)])
                pos = 0
                while pos < Xs.shape[0] - 1 and self.loose_qcutoff > Ys[pos + 1]:
                    pos += 1
                if self.loose_qcutoff < Ys[0]:
                    c1.text(0, -0.45, "Loose quality maximum resolution is above the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                elif self.loose_qcutoff > numpy.amax(Ys):
                    c1.text(0, -0.45, "Loose quality maximum resolution is below the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (numpy.log10(Xs[pos + 1]) - numpy.log10(Xs[pos]))
                    B = Ys[pos] - numpy.log10(Xs[pos]) * A
                    c1.text(0, -0.45, "Loose quality maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.loose_qcutoff - B) / A)))),
                        [text.halign.left, text.valign.top, text.size(-2)])
                c.insert(c1, [trafo.translate(0, H - width)])
                H -= width + 1.1

        if rscores is not None:
            H -= 0.4
            chromosomes = self.storage['replicate_chromosomes'][...]
            resolutions = self.storage['replicate_resolutions'][...]
            coverages = self.storage['replicate_coverages'][...]
            c.text(0, H, "Replicate Score Results", [text.halign.left, text.valign.bottom, text.size(0)])
            H -= 0.6
            hoffset = 1.7
            for i in range(coverages.shape[0] + 1):
                c.stroke(path.line(hoffset + i * 2.0, H + 0.3, hoffset + i * 2.0, H - 0.3 * resolutions.shape[0]))
                if i < coverages.shape[0]:
                    c.text(hoffset + (i + 0.5) * 2.0, H + 0.05, self._num2str(coverages[i]),
                        [text.halign.center, text.valign.bottom, text.size(-2)])
            for i in range(resolutions.shape[0] + 1):
                c.stroke(path.line(0, H - 0.3 * i, hoffset + coverages.shape[0] * 2.0, H - 0.3 * i))
            c.text(hoffset + 1.0 * coverages.shape[0], H + 0.35, "Coverage",
                [text.halign.center, text.valign.bottom, text.size(-2)])
            c.text(hoffset - 0.1, H + 0.05, "Resolution",
                [text.halign.right, text.valign.bottom, text.size(-2)])
            H -= 0.3
            for i, res in enumerate(resolutions):
                if res < 1000:
                    label = "%i bp" % res
                elif res < 1000000:
                    label = "%i Kb" % (res / 1000)
                else:
                    label = "%i Mb" % (res / 1000000)
                c.text(hoffset - 0.1, H + 0.05, label, [text.halign.right, text.valign.bottom, text.size(-2)])
                for j, cov in enumerate(coverages):
                    c.text(hoffset + (j + 0.5) * 2.0, H + 0.05, '%0.6f' % rscores[-1, i, j],
                        [text.halign.center, text.valign.bottom, text.size(-2)])
                H -= 0.3
            H -= 0.4
            hoffset = 0.9

            # if there are sufficient data points, plot coverage vs score
            if coverages.shape[0] >= 2:
                width = 17.0 / min(3, resolutions.shape[0]) - 0.3
                Xs = coverages
                for i, res in enumerate(resolutions):
                    if res < 1000:
                        label = "%i bp" % res
                    elif res < 1000000:
                        label = "%i Kb" % (res / 1000)
                    else:
                        label = "%i Mb" % (res / 1000000)
                    Ys = rscores[-1, i, :]
                    c1 = self._plot_graph(Xs, Ys, width, 'r')
                    c1.text(hoffset, -0.15, "Resolution: %s" % label,
                        [text.halign.left, text.valign.top, text.size(-2)])
                    c.insert(c1, [trafo.translate((width + 0.3) * i, H - width)])
                H -= width + 0.8

            # if there are multiple coverages, estimate maximum usable resolution
            if resolutions.shape[0] > 1:
                Xs = resolutions
                lXs = numpy.log10(resolutions)
                Ys = rscores[-1, :, -1]
                width = 17.0 / 3.0 - 0.3
                c1 = self._plot_graph(Xs, Ys, width, 'r', 'res')
                pos = 0
                while pos < Xs.shape[0] - 1 and self.strict_rcutoff > Ys[pos + 1]:
                    pos += 1
                if self.strict_rcutoff < Ys[0]:
                    c1.text(0, -0.15, "Strict replicate maximum resolution is above the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                elif self.strict_rcutoff > numpy.amax(Ys):
                    c1.text(0, -0.15, "Strict replicate maximum resolution is below the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (lXs[pos + 1] - lXs[pos])
                    B = Ys[pos] - lXs[pos] * A
                    c1.text(0, -0.15, "Strict replicate maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.strict_rcutoff - B) / A)))),
                        [text.halign.left, text.valign.top, text.size(-2)])
                pos = 0
                while pos < Xs.shape[0] - 1 and self.loose_rcutoff > Ys[pos + 1]:
                    pos += 1
                if self.loose_rcutoff < Ys[0]:
                    c1.text(0, -0.45, "Loose replicate maximum resolution is above the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                elif self.loose_rcutoff > numpy.amax(Ys):
                    c1.text(0, -0.45, "Loose replicate maximum resolution is below the resolutions tested.",
                        [text.halign.left, text.valign.top, text.size(-2)])
                else:
                    A = (Ys[pos + 1] - Ys[pos]) / (lXs[pos + 1] - lXs[pos])
                    B = Ys[pos] - lXs[pos] * A
                    c1.text(0, -0.45, "Loose replicate maximum resolution: %s bp" % self._num2str(
                        int(round(10.0 ** ((self.loose_rcutoff - B) / A)))),
                        [text.halign.left, text.valign.top, text.size(-2)])
                c.insert(c1, [trafo.translate(0, H - width)])
        c.writePDFfile(filename)
        return None

    def _plot_graph(self, Xs, Ys, width, dtype, xaxis='cov'):
        hoffset = 0.9
        voffset = 0.6
        pwidth = width - hoffset
        pheight = width - voffset
        c = canvas.canvas()
        minX = numpy.log2(numpy.amin(Xs))
        maxX = numpy.log2(numpy.amax(Xs))
        if xaxis == 'cov':
            minY = numpy.amin(Ys)
            maxY = numpy.amax(Ys)
        else:
            if dtype == 'q':
                minY = min(numpy.amin(Ys), self.loose_qcutoff)
                maxY = max(numpy.amax(Ys), self.strict_qcutoff)
            else:
                minY = min(numpy.amin(Ys), self.loose_rcutoff)
                maxY = max(numpy.amax(Ys), self.strict_rcutoff)
        spanX = maxX - minX
        spanY = maxY - minY
        minY -= 0.05 * spanY
        maxY += 0.05 * spanY
        spanY = maxY - minY
        c1 = self._plot_line(Xs, Ys, width, minX, maxX, minY, maxY)
        if dtype == 'q':
            ylab = 'Quality Score'
        else:
            ylab = 'Replicate Score'
        c.text(0, voffset + pheight * 0.5, ylab,
            [text.halign.center, text.valign.top, text.size(-2), trafo.rotate(90)])
        if xaxis == 'cov':
            if numpy.amin(Xs) < 1000000:
                start = -1
            else:
                start = 0
            stop = int(floor(maxX - numpy.log2(1000000)))
            for i in range(start, stop + 1, 2):
                val = 2 ** i
                X = (i + numpy.log2(1000000) - minX) / spanX * pwidth + hoffset
                c.stroke(path.line(X, voffset, X, voffset - 0.08))
                c.stroke(path.line(X, voffset, X, width), [color.gray(0.9)])
                if val < 0.1:
                    label = '%0.2f' % val
                elif val < 1.0:
                    label = '%0.1f' % val
                else:
                    label = '%i' % val
                c.text(X, voffset - 0.1, label, [text.halign.center, text.valign.top, text.size(-2)])
            for i in range(start - 2, stop + 3, 2):
                xs = numpy.log2(numpy.linspace(2 ** i * 1000000, 2 ** (i + 2) * 1000000, 9)[1:-1])
                xs = (xs - minX) / spanX * pwidth + hoffset
                for x in xs:
                    if x > hoffset and x < width:
                        c.stroke(path.line(x, voffset, x, voffset - 0.05))
            label = "Millions of Reads"
        else:
            for x in Xs:
                X = (numpy.log2(x) - minX) / spanX * pwidth + hoffset
                if x < 1000000:
                    label = '%iK' % (x / 1000)
                else:
                    label = '%iM' % (x / 1000000)
                c.stroke(path.line(X, voffset, X, voffset - 0.05))
                c.text(X, voffset - 0.1, label, [text.halign.center, text.valign.top, text.size(-2)])
                if X > hoffset and X < width:
                    c.stroke(path.line(X, voffset, X, width), [color.gray(0.9)])
            label = "Resolution"
        c.text(hoffset + pwidth * 0.5, 0, label, [text.halign.center, text.valign.bottom, text.size(-2)])
        scale = 1.0
        while maxY * 10 ** scale < 1.0:
            scale += 1
        step = 1.0
        while (floor(maxY * 10 ** scale) - ceil(minY * 10 ** scale)) / step > 5:
            step += 1.0
        N = (floor(maxY * 10 ** scale) - ceil(minY * 10 ** scale)) / step
        for i in numpy.linspace(floor(maxY * 10 ** scale) / 10 ** scale, ceil(minY * 10 ** scale) / 10 ** scale,
            int(floor((floor(maxY * 10 ** scale) - ceil(minY * 10 ** scale)) / step)) + 1):
            Y = (i - minY) / spanY * pheight + voffset
            c.stroke(path.line(hoffset, Y, hoffset - 0.05, Y))
            if Y > voffset and Y < width:
                c.stroke(path.line(hoffset, Y, width, Y), [color.gray(0.9)])
            if i < 0.01:
                label = '%0.3f' % i
            if i < 0.1:
                label = '%0.2f' % i
            else:
                label = '%0.1f' % i
            c.text(hoffset - 0.1, Y, label, [text.halign.right, text.valign.middle, text.size(-2)])
        if xaxis == 'res':
            if dtype == 'q':
                Y1 = (self.strict_qcutoff - minY) / spanY * pheight + voffset
                Y2 = (self.loose_qcutoff - minY) / spanY * pheight + voffset
            else:
                Y1 = (self.strict_rcutoff - minY) / spanY * pheight + voffset
                Y2 = (self.loose_rcutoff - minY) / spanY * pheight + voffset
            c.stroke(path.line(hoffset, Y1, width, Y1), [style.linestyle.dashed])
            c.stroke(path.line(hoffset, Y2, width, Y2), [style.linestyle.dotted])
        c.insert(c1)
        c.stroke(path.rect(hoffset, voffset, pwidth, pheight))
        return c

    def _plot_line(self, Xs, Ys, width, minX, maxX, minY, maxY, pcolor=None):
        if pcolor is None:
            pcolor = color.rgb.black
        spanX = maxX - minX
        spanY = maxY - minY
        hoffset = 0.9
        voffset = 0.6
        pwidth = width - hoffset
        pheight = width - voffset
        c = canvas.canvas([canvas.clip(path.rect(hoffset, voffset, pwidth, pheight))])
        order = numpy.argsort(Xs)
        xs = (numpy.log2(Xs[order]) - minX) / spanX * pwidth + hoffset
        ys = (Ys[order] - minY) / spanY * pheight + voffset
        lpath = path.path(path.moveto(xs[0], ys[0]))
        for i in range(1, xs.shape[0]):
            lpath.append(path.lineto(xs[i], ys[i]))
        c.stroke(lpath, [pcolor])
        return c

    def _num2str(self, n):
        s = []
        n1 = str(n)
        while len(n1) > 3:
            s = [n1[-3:]] + s
            n1 = n1[:-3]
        if len(n1) > 0:
            s = [n1] + s
        return ','.join(s)
