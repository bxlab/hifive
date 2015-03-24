#!/usr/bin/env python

import sys

try:
    from mpi4py import MPI
except:
    pass

from ..hic import HiC


def run(args):
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    if rank == 0:
        hic = HiC(args.output, 'w', silent=args.silent)
        hic.load_data(args.data)
        hic.filter_fends(mininteractions=args.minint, mindistance=args.mindist, maxdistance=args.maxdist)
        hic.save()
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
    else:
        comm.recv(source=0, tag=11)
        hic = HiC(args.output, 'r', silent=True)
    hic.find_distance_parameters(minsize=args.minbin, numbins=args.numbins)
    if rank == 0:
        hic.save()
    return None
