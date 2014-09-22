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
    if len(sys.argv) < 8:
        print "Usage python create_hic_set.py DATA_FILE HIC_FILE MIN_INTERACTIONS MAX_DIST MIN_SIZE NUM_BINS SMOOTHED"
        print "DATA_FILE         File name of HiCData h5dict to link with analysis."
        print "OUT_FILE          File name to write HiC h5dict to."
        print "MIN_INTERACTIONS  Minimum number of interactions needed for valid fends."
        print "MAX_DIST          The largest interaction distance to be included for filtering fends."
        print "MIN_SIZE          The smallest interaction distance bin size for distance function."
        print "NUM_BINS          The number of bins to partion interaction distance range into for distance function."
        print "SMOOTHED          Number of adjacent bins to include for smoothing of distance function line."
        print "This function is MPI compatible."
        return None
    elif len(sys.argv) < 8:
        return None
    data_fname, hic_fname, mininteractions, maxdistance, minsize, numbins, smoothed = sys.argv[1:8]
    mininteractions, maxdistance, minsize, numbins, smoothed = (
        int(mininteractions), int(maxdistance), int(minsize), int(numbins), int(smoothed))
    if rank == 0:
        hic = hifive.HiC(hic_fname, 'w')
        hic.load_data(data_fname)
        hic.filter_fends(mininteractions=mininteractions, maxdistance=maxdistance)
        hic.save()
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
    else:
        comm.recv(source=0, tag=11)
        hic = hifive.HiC(hic_fname, 'r')
    hic.find_distance_means(minsize=minsize, smoothed=smoothed, numbins=numbins)
    if rank == 0:
        hic.save()
    return None

    
if __name__ == "__main__":
    main()
