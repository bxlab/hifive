#!/usr/bin/env python

import sys
import os

import numpy
from math import ceil, log, exp
import h5py
try:
    from mpi4py import MPI
except:
    pass

import hifive


def main():
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    if len(sys.argv) < 6 and rank == 0:
        print "Usage: python combine_BIs.py BI_FILE_1 BI_FILE_2 OUT_FILE SMOOTHING CHROM_1[,CHROM_2...,CHROM_N]"
        print "BI_FILE_1           h5dicts created by the hifive.BI class"
        print "BI_FILE_2           h5dict created by the hifive.BI class"
        print "OUT_FILE            file name for the new h5dict created by this script"
        print "SMOOTHING           integer specifying the width of smoothing weights"
        print "CHROM1[,CHROM2...]  a comma-separated list of chromosome names to include in processing"
        print "This function is MPI compatible."
        return None
    elif len(sys.argv < 6):
        return None
    bi_fname1, bi_fname2, out_fname, smoothing, chroms = sys.argv[1:6]
    chroms = chroms.split(',')
    smoothing = int(smoothing)
    if rank == 0:
        bi1 = hifive.BI()
        bi1.load(bi_fname1)
        bi2 = hifive.BI()
        bi2.load(bi_fname2)
        BI = hifive.BI(bi1.width, bi1.window, bi1.height, bi1.mincount)
        BI.chromosomes = numpy.array(chroms)
        BI.chr2int = {}
        for i, chrom in enumerate(BI.chromosomes):
            BI.chr2int[chrom] = i
        all_BI = numpy.zeros(0, dtype=numpy.dtype([('chromosome', numpy.int32),
                             ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32),
                             ('score', numpy.float32)]))
        for chrom in chroms:
            if chrom not in bi1.chr2int:
                continue
            chrint = bi1.chr2int[chrom]
            start = bi1.chr_indices[chrint]
            stop = bi1.chr_indices[chrint + 1]
            temp_BI = numpy.zeros(stop - start, dtype=numpy.dtype([('chromosome', numpy.int32),
                             ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32),
                             ('score', numpy.float32)]))
            temp_BI['chromosome'][:] = BI.chr2int[chrom]
            temp_BI['start'][:] = bi1.BI['start'][start:stop]
            temp_BI['stop'][:] = bi1.BI['stop'][start:stop]
            temp_BI['mid'][:] = bi1.BI['mid'][start:stop]
            if 'original' in bi1.BI.dtype.names:
                temp_score = bi1.BI['original'][start:stop]
            else:
                temp_score = bi1.BI['score'][start:stop]
            temp_score -= numpy.mean(temp_score)
            temp_score /= numpy.std(temp_score)
            temp_BI['score'][:] = temp_score
            all_BI = numpy.hstack((all_BI, temp_BI))
        for chrom in chroms:
            if chrom not in bi2.chr2int:
                continue
            chrint = bi2.chr2int[chrom]
            start = bi2.chr_indices[chrint]
            stop = bi2.chr_indices[chrint + 1]
            temp_BI = numpy.zeros(stop - start, dtype=numpy.dtype([('chromosome', numpy.int32),
                             ('start', numpy.int32), ('stop', numpy.int32), ('mid', numpy.int32),
                             ('score', numpy.float32)]))
            temp_BI['chromosome'][:] = BI.chr2int[chrom]
            temp_BI['start'][:] = bi2.BI['start'][start:stop]
            temp_BI['stop'][:] = bi2.BI['stop'][start:stop]
            temp_BI['mid'][:] = bi2.BI['mid'][start:stop]
            if 'original' in bi2.BI.dtype.names:
                temp_score = bi2.BI['original'][start:stop]
            else:
                temp_score = bi2.BI['score'][start:stop]
            temp_score -= numpy.mean(temp_score)
            temp_score /= numpy.std(temp_score)
            temp_BI['score'][:] = temp_score
            all_BI = numpy.hstack((all_BI, temp_BI))
        order = numpy.lexsort((all_BI['mid'], all_BI['chromosome']))
        all_BI = all_BI[order]
        where = numpy.where((all_BI['chromosome'][:-1] == all_BI['chromosome'][1:]) *
                            (all_BI['mid'][:-1] == all_BI['mid'][1:]))[0]
        all_BI['start'][where + 1] = (all_BI['start'][where] + all_BI['start'][where + 1]) / 2
        all_BI['stop'][where + 1] = (all_BI['stop'][where] + all_BI['stop'][where + 1]) / 2
        all_BI['score'][where + 1] = (all_BI['score'][where] + all_BI['score'][where + 1]) / 2
        valid = [0] + list(numpy.where(all_BI['mid'][:-1] != all_BI['mid'][1:])[0] + 1)
        all_BI = all_BI[valid]
        BI.BI = all_BI
        BI.chr_indices = numpy.zeros(BI.chromosomes.shape[0] + 1, dtype=numpy.int32)
        BI.chr_indices[1:] = numpy.bincount(BI.BI['chromosome'][:], minlength=BI.chromosomes.shape[0])
        for i in range(1, BI.chr_indices.shape[0]):
            BI.chr_indices[i] += BI.chr_indices[i - 1]
        BI.save(out_fname)
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
    else:
        comm.recv(source=0, tag=11)
        BI = hifive.bi.BI()
        BI.load(out_fname)
    if smoothing > 0:
        BI.smooth_bi(smoothing)
    if rank == 0:
        BI.save(out_fname)
    return None


if __name__ == "__main__":
    main()
