#!/usr/bin/env python

import sys

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
        print "Usage python learn_hic_normalization_express.py HIC_FILE ITERATIONS MIN_INT MIN_DIST USE_READS REMOVE_DISTANCE RECALC"
        print "HIC_FILE         File name of HiC h5dict to analyze."
        print "ITERATIONS       Number of iterations to run learning for."
        print "MIN_INT          Minimum number of interactions for fend filtering, if refiltering is required."
        print "MIN_DIST         Minimum interaction distance to include for learning."
        print "USE_READS        Which set of reads, 'cis', 'trans', or 'both', to use for learning."
        print "REMOVE_DISTANCE  Specifies whether to remove distance-dependent portion of the signal prior to learning."
        print "RECALC           Number of iterations to wait between recalculating distance function parameters."
        print "This function is MPI compatible."
        return None
    elif len(sys.argv) < 8:
        return None
    hic_fname, iterations, mininteractions, mindistance, usereads, remove_distance, recalc = sys.argv[1:8]
    if remove_distance in ['1', 'true', 'True', 'TRUE']:
        remove_distance = True
    else:
        remove_distance = False
    hic = hifive.HiC(hic_fname, 'r')
    hic.find_express_fend_corrections(iterations=int(iterations), mindistance=int(mindistance),
                                      mininteractions=int(mininteractions), usereads=usereads,
                                      remove_distance=remove_distance, recalculate_distance=int(recalc))
    if rank == 0:
        hic.save()


if __name__ == "__main__":
    main()
