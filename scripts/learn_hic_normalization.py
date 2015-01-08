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
        print "Usage python learn_hic_normalization.py HIC_FILE BURNIN ANNEALING MIN_CHANGE MIN_DIST MAX_DIST RATE DISPLAY CHROMS"
        print "HIC_FILE    File name of HiC h5dict to analyze."
        print "BURNIN      Number of iterations to run burn-in phase for."
        print "ANNEALING   Number of iterations to run annealing phase for."
        print "MIN_CHANGE  The minimum mean change in fend correction parameter values needed to keep running past 'burnin_iterations' number of iterations during burn-in phase."
        print "MIN_DIST    Minimum interaction distance to include in learning."
        print "MAX_DIST    Maximum interaction distance to include in learning."
        print "RATE        Percent of gradient to use for updating parameter values."
        print "DISPLAY     Number of iterations to wait before explicitly calculating cost and updating display."
        print "CHROMS      A comma-separated list of chromosome names to learn corrections for"
        print "This function is MPI compatible."
        return None
    elif len(sys.argv) < 9:
        return None
    hic_fname, burnin_iterations, annealing_iterations, minchange, mindist, maxdist, learning_rate, display = (
        sys.argv[1:9])
    if len(sys.argv) > 9:
        chroms = sys.argv[9].split(',')
    else:
        chroms = []
    hic = hifive.HiC(hic_fname, 'r')
    hic.find_fend_corrections(display=int(display), mindistance=int(mindist), maxdistance=int(maxdist), chroms=chroms,
                              learningrate=float(learning_rate), burnin_iterations=int(burnin_iterations),
                              annealing_iterations=int(annealing_iterations), minchange=float(minchange))
    if rank == 0:
        hic.save()


if __name__ == "__main__":
    main()
