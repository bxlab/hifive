#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys

import numpy
try:
    from mpi4py import MPI
except:
    pass

import hifive


def find_BI(hic_fname, BI_fname, width=10000, window=2500000, height=0, mincount=10, smoothing=10000, chroms=[]):
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    hic = hifive.analysis.HiC(hic_fname, 'r')
    BI = hifive.bi.BI(width=width, window=window, height=height, mincount=10)
    BI.find_bi_from_hic(hic,  datatype='enrichment', chroms=chroms)
    if smoothing > 0:
        BI.smooth_bi(smoothing)
    if rank == 0:
        BI.save(BI_fname)
    return None


if __name__ == '__main__':
    hic_fname, BI_fname, width, window, height, mincount, smoothing = sys.argv[1:8]
    width, window, height, mincount, smoothing = int(width), int(window), int(height), int(mincount), int(smoothing)
    if len(sys.argv) > 8:
        chroms = sys.argv[8].split(',')
    else:
        chroms = []
    find_BI(hic_fname, BI_fname, width, window, height, mincount, smoothing, chroms)
