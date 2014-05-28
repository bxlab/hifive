#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

import sys
import os

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
    data_fname, hic_fname = sys.argv[1:3]
    if rank == 0:
        hic = hifive.analysis.HiC(hic_fname, 'w')
        hic.load_data(data_fname)
        if len(sys.argv) >= 4:
            mininteractions = int(sys.argv[3])
            maxdistance = int(sys.argv[4])
            hic.filter_fends(mininteractions=mininteractions, maxdistance=maxdistance)
        else:
            hic.filter_fends()
        hic.save()
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
    else:
        comm.recv(source=0, tag=11)
        hic = hifive.analysis.HiC(hic_fname, 'r')
    if len(sys.argv) >= 8:
        minsize = int(sys.argv[5])
        numbins = int(sys.argv[6])
        smoothed = int(sys.argv[7])
        hic.find_distance_means(minsize=minsize, smoothed=smoothed, numbins=numbins)
    else:
        hic.find_distance_means()
    if rank == 0:
        hic.save()
    return None

    
if __name__ == "__main__":
    main()
