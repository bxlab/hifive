#!/usr/bin/env python

from ..hic_data import HiCData


def run(args):
    data = HiCData(args.output, 'w', silent=args.silent)
    if not args.bam is None: 
        data.load_data_from_bam(args.fend, args.bam, args.insert, args.skipdups)
    elif not args.raw is None: 
        data.load_data_from_raw(args.fend, args.raw, args.insert, args.skipdups)
    elif not args.mat is None: 
        data.load_data_from_mat(args.fend, args.mat)
    elif not args.matrix is None:
        data.load_binned_data_from_matrices(args.fend, args.matrix, format=None)
    data.save()
