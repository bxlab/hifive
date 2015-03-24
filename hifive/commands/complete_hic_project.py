#!/usr/bin/env python

import sys

try:
    from mpi4py import MPI
except:
    pass

from ..fend import Fend
from ..hic_data import HiCData
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
        modelbins = args.modelbins.split(',')
        parameters.split(',')
        for i in range(len(modelbins)):
            try:
                modelbins[i] = int(modelbins[i])
            except:
                if rank == 0:
                    print sys.stderr, ("Not all arguments in -n/--modelbins could be converted to integers.")
                return 1
        if len(model) != len(modelbins) or len(model) != len(parameters):
            if rank == 0:
                print sys.stderr, ("-v/--model, -n/--modelbins, and -u/--parameter-types must be equal lengths.")
            return 1
    if args.prefix is None:
        fend_fname, data_fname, project_fname = args.output
    else:
        fend_fname = "%s.fends" % args.prefix
        data_fname = "%s.hcd" % args.prefix
        project_fname = "%s.hcp" % args.prefix
    if rank == 0:
        fends = Fend(fend_fname, 'w', silent=args.silent)
        if args.bed is None:
            fends.load_fends(args.fend, genome_name=args.genome, re_name=args.re, format="fend")
        else:
            fends.load_fends(args.bed, genome_name=args.genome, re_name=args.re, format="bed")
        fends.save()
        del fends
        data = HiCData(data_fname, 'w', silent=args.silent)
        if not args.bam is None: 
            data.load_data_from_bam(args.fend, args.bam, args.insert)
        elif not args.raw is None: 
            data.load_data_from_raw(args.fend, args.raw, args.insert)
        elif not args.mat is None: 
            data.load_data_from_mat(args.fend, args.mat, args.insert)
        data.save()
        del data
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
    else:
        comm.recv(source=0, tag=11)
    hic = HiC(project_fname, 'w', silent=args.silent)
    hic.load_data(data_fname)
    hic.filter_fends(mininteractions=args.minint, mindistance=args.mindist, maxdistance=args.maxdist)
    hic.find_distance_parameters(minsize=args.minbin, numbins=args.numbins)
    precorrect = False
    if args.algorithm in ['binning', 'binning-express', 'binning-probability']:
        hic.find_binning_fend_corrections(mindistance=args.mindist, maxdistance=args.maxdist, parameters=parameters,
                                             chroms=chroms, num_bins=modelbins, model=model, usereads=args.binreads,
                                             learning_threshold=args.threshold, max_iterations=args.biniter)
        precorrect = True
    if args.algorithm in ['probability', 'binning-probability']:
        hic.find_probability_fend_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                              minchange=args.change, burnin_iterations=args.burnin,
                                              annealing_iterations=args.anneal, learningrate=args.rate,
                                              display=args.display, chroms=chroms,
                                              precalculate=args.precalc, precorrect=precorrect)
    elif args.algorithm in ['express', 'binning-express']:
        hic.find_express_fend_corrections(iterations=args.expiter, mindistance=args.mindist,
                                          maxdistance=args.maxdist, remove_distance=args.nodist,
                                          usereads=args.expreads, mininteractions=args.minint,
                                          chroms=chroms, precorrect=precorrect)
    if rank == 0:
        hic.save()
