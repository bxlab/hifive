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
    if args.chroms is None:
        chroms = []
    else:
        chroms = args.chroms.split(',')
        if len(chroms) == 1 and chroms[0] == '':
            chroms = []
    if args.algorithm.count('binning') > 0:
        model = args.model.split(',')
        modelbins = args.modelbins.split(',')
        parameters = args.parameters.split(',')
        for i in range(len(modelbins)):
            try:
                modelbins[i] = int(modelbins[i])
            except:
                if rank == 0:
                    print sys.stderr, ("Not all arguments in -n/--modelbins could be converted to integers.")
                return 1
        if len(model) != len(modelbins):
            if rank == 0:
                print sys.stderr, ("-v/--model, -n/--modelbins, and -u/--parameter-types must be equal lengths.")
            return 1
    hic = HiC(args.project, 'r', silent=args.silent)
    precorrect = False
    if args.algorithm in ['binning', 'binning-express', 'binning-probability']:
        hic.find_binning_fend_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                          chroms=chroms, num_bins=modelbins, model=model, parameters=parameters,
                                          usereads=args.binreads, learning_threshold=args.threshold,
                                          max_iterations=args.biniter, pseudocounts=args.pseudo)
        precorrect = True
    if args.algorithm in ['probability', 'binning-probability']:
        hic.find_probability_fend_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                              minchange=args.change, max_iterations=args.probiter,
                                              learningstep=args.step, chroms=chroms,
                                              precalculate=args.precalc, precorrect=precorrect,
                                              model=args.probmodel)
    elif args.algorithm in ['express', 'binning-express']:
        hic.find_express_fend_corrections(iterations=args.expiter, mindistance=args.mindist,
                                          maxdistance=args.maxdist, remove_distance=args.nodist,
                                          usereads=args.expreads, mininteractions=args.minint,
                                          chroms=chroms, minchange=args.change, precorrect=precorrect,
                                          binary=args.binary, kr=args.kr)
    if rank == 0:
        hic.save(args.output)
