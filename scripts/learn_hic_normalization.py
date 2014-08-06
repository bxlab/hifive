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
hic_fname, burnin_iterations, annealing_iterations, correctionmaxdistance, recalc, learning_rate, display = sys.argv[1:8]
hic = hifive.analysis.HiC(hic_fname, 'r')
hic.find_fend_corrections(display=int(display), maxdistance=int(correctionmaxdistance),
                          learningrate=float(learning_rate), burnin_iterations=int(burnin_iterations),
                          annealing_iterations=int(annealing_iterations), recalculate_distance=int(recalc))
if rank == 0:
    hic.save()
