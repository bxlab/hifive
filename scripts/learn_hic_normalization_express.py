#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

import sys

try:
    from mpi4py import MPI
except:
    pass

import hifive


if 'mpi4py' in sys.modules.keys():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()
else:
    comm = None
    rank = 0
    num_procs = 1
hic_fname, iterations, mininteractions, mindistance, usereads, removedistance, recalc = sys.argv[1:8]
if removedistance in ['1', 'true', 'True', 'TRUE']:
    removedistance = True
else:
    removedistance = False
hic = hifive.analysis.HiC(hic_fname, 'r')
hic.find_express_fend_corrections(iterations=int(iterations), mindistance=int(mindistance),
                                  mininteractions=int(mininteractions), usereads=usereads,
                                  removedistance=removedistance, recalculate_distance=int(recalc))
if rank == 0:
    hic.save()
