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
    chroms = args.chroms.split(',')
    if len(chroms) == 1 and chroms[0] == '':
        chroms = []
    if not args.model is None:
        model = args.model.split(',')
        for par in model:
            if par not in ['gc', 'len', 'distance']:
                if rank == 0:
                    print sys.stderr, ("Not all arguments in -v/--model are valid.")
                return 1
        modelbins = args.modelbins.split(',')
        for i in range(len(modelbins)):
            try:
                modelbins[i] = int(modelbins[i])
            except:
                if rank == 0:
                    print sys.stderr, ("Not all arguments in -n/--modelbins could be converted to integers.")
                return 1
        if len(model) != len(modelbins):
            if rank == 0:
                print sys.stderr, ("-v/--model and -n/--modelbins not equal lengths.")
            return 1
    hic = HiC(args.project, 'r', silent=args.silent)
    precorrect = False
    if args.algorithm in ['regression', 'regression-express', 'regression-probability']:
        hic.find_regression_fend_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                             chroms=chroms, num_bins=modelbins, model=model, usereads=args.regreads,
                                             learning_threshold=args.threshold, max_iterations=args.regiter)
        precorrect = True
    if args.algorithm in ['probability', 'regression-probability']:
        hic.find_probability_fend_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                              minchange=args.change, burnin_iterations=args.burnin,
                                              annealing_iterations=args.anneal, learningrate=args.rate,
                                              display=args.display, chroms=chroms,
                                              precalculate=args.precalc, precorrect=precorrect)
    elif args.algorithm in ['express', 'regression-express']:
        hic.find_express_fend_corrections(iterations=args.expiter, mindistance=args.mindist,
                                          maxdistance=args.maxdist, remove_distance=args.nodist,
                                          usereads=args.expreads, mininteractions=args.minint,
                                          chroms=chroms, precorrect=precorrect)
    if rank == 0:
        hic.save(args.outfname)
