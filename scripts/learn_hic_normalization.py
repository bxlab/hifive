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
        print "Usage python learn_hic_normalization.py HIC_FILE BURNIN ANNEALING MAX_DIST RECALC RATE DISPLAY"
        print "HIC_FILE    File name of HiC h5dict to analyze."
        print "BURNIN      Number of iterations to run burn-in phase for."
        print "ANNEALING   Number of iterations to run annealing phase for."
        print "MAX_DIST    Maximum interaction distance to include in learning."
        print "RECALC      Number of iterations to wait between recalculating distance function parameters."
        print "RATE        Percent of gradient to use for updating parameter values."
        print "DISPLAY     Number of iterations to wait before explicitly calculating cost and updating display."
        print "This function is MPI compatible."
        return None
    elif len(sys.argv) < 8:
        return None
    hic_fname, burnin_iterations, annealing_iterations, correctionmaxdistance, recalc, learning_rate, display = (
        sys.argv[1:8])
    hic = hifive.HiC(hic_fname, 'r')
    hic.find_fend_corrections(display=int(display), maxdistance=int(correctionmaxdistance),
                              learningrate=float(learning_rate), burnin_iterations=int(burnin_iterations),
                              annealing_iterations=int(annealing_iterations), recalculate_distance=int(recalc))
    if rank == 0:
        hic.save()


__name__ = "__main__":
    main()
