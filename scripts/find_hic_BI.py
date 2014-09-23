#!/usr/bin/env python

import sys

import numpy
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
    if len(sys.argv) < 8 and rank == 0:
        print "Usage: python find_hic_BI.py HIC_FILE OUT_FILE WIDTH WINDOW HEIGHT MINCOUNT SMOOTHING [CHROM_1,CHROM_2...,CHROM_N]"
        print "HIC_FILE            h5dict created by the hifive.HiC class"
        print "OUT_FILE            file name for the new h5dict created by this script"
        print "WIDTH               integer specifying the width about each boundary point"
        print "HEIGHT              integer specifying the height of bins extending across each window"
        print "WINDOW              integer specifying the window around each boundary point"
        print "MINCOUNT            minimum number of valid bin pairs needed to find BI value"
        print "SMOOTHING           integer specifying the width of smoothing weights"
        print "[CHROM1,CHROM2...]  a comma-separated list of chromosome names to include in processing"
        print "This function is MPI compatible."
        return None
    elif len(sys.argv) < 8:
        return None
    hic_fname, BI_fname, width, height, window, mincount, smoothing = sys.argv[1:8]
    width, window, height, mincount, smoothing = int(width), int(window), int(height), int(mincount), int(smoothing)
    if len(sys.argv) > 8:
        chroms = sys.argv[8].split(',')
    else:
        chroms = []
    hic = hifive.HiC(hic_fname, 'r')
    BI = hifive.BI(width=width, window=window, height=height, mincount=10)
    BI.find_bi_from_hic(hic,  datatype='enrichment', chroms=chroms)
    if smoothing > 0:
        BI.smooth_bi(smoothing)
    if rank == 0:
        BI.save(BI_fname)
    return None


if __name__ == '__main__':
    main()
