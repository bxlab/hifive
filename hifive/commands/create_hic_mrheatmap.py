#!/usr/bin/env python

import sys
import subprocess

import h5py

from ..hic import HiC


def run(args):
    if args.chroms is None:
        chroms = None
    else:
        chroms = args.chroms.split(',')
        if len(chroms) == 1 and chroms[0] == '':
            chroms = None
    hic = HiC(args.project, 'r', silent=args.silent)
    hic.write_multiresolution_heatmap(args.output, datatype=args.datatype, maxbinsize=args.maxbin,
                                      minbinsize=args.minbin, trans_maxbinsize=args.maxtransbin,
                                      trans_minbinsize=args.mintransbin, minobservations=args.minobs,
                                      chroms=chroms, includetrans=args.trans, midbinsize=args.midbin)
