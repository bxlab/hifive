#!/usr/bin/env python

import sys
import subprocess

try:
    from mpi4py import MPI
except:
    pass
try:
    from pyx import *
except:
    pass
import h5py

from ..hic import HiC
from ..plotting import plot_hic_heatmap, plot_key


def run(args):
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    if not args.image is None and args.pdf and "matplotlib" not in sys.modules.keys():
        if rank == 0:
            print("-p/--pdf requires the package 'matplotlib'", file=sys.stderr)
        sys.exit(1)
    if args.chroms is None:
        chroms = []
    else:
        chroms = args.chroms.split(',')
        if len(chroms) == 1 and chroms[0] == '':
            chroms = []
    hic = HiC(args.project, 'r', silent=args.silent)
    hic.write_heatmap(args.output, binsize=args.binsize, includetrans=args.trans,
                      datatype=args.datatype, chroms=chroms, dynamically_binned=args.dynamic,
                      expansion_binsize=args.expbinsize, minobservations=args.minobs,
                      searchdistance=args.search, removefailed=args.remove, format=args.format)
    if rank > 0:
        sys.exit(0)
    if not args.image is None:
        if args.format == 'txt':
            if rank == 0:
                print("Plotting is only available for non-txt formats.", file=sys.stderr)
            return None
        kwargs = {}
        for arg in args.keywords:
            temp = arg.split("=")
            if temp[1] in ["True", "TRUE", "true"]:
                temp[1] = True
            elif temp[1] in ["False", "FALSE", "false"]:
                temp[1] = False
            elif temp[1][0] == "(":
                temp[1] = temp[1].strip('()').split(',')
                for i in range(len(temp[1])):
                    temp[1][i] = int(temp[1][i])
                temp[1] = tuple(temp[1])
            elif temp[1][0] == "[":
                temp[1] = temp[1].strip('[]').split(',')
                for i in range(len(temp[1])):
                    temp[1][i] = int(temp[1][i])
            else:
                try:
                    temp[1] = int(temp[1])
                except:
                    try:
                        temp[1] = float(temp[1])
                    except:
                        # strip off extra characters introduced by galaxy into color format
                        temp[1] = temp[1].replace('__pd__','')
            kwargs[temp[0]] = temp[1]
        if 'symmetricscaling' not in kwargs:
            if args.datatype == 'enrichment':
                kwargs['symmetricscaling'] = True
            else:
                kwargs['symmetricscaling'] = False
        minscore, maxscore = plot_hic_heatmap(args.output, returnscale=True, silent=args.silent,
                                              format=args.format, **kwargs)

