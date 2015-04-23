#!/usr/bin/env python

import sys

from ..fragment import Fragment
from ..fivec_data import FiveCData
from ..fivec import FiveC


def run(args):
    if args.regions is None:
        regions = []
    else:
        regions = args.regions.split(',')
        if len(regions) == 1 and regions[0] == '':
            regions = []
        for i in range(len(regions)):
            try:
                regions[i] = int(regions[i])
            except:
                print sys.stderr, ("Not all arguments in -r/--regions could be converted to integers.")
                return 1
    if args.algorithm.count('binning') > 0:
        model = args.model.split(',')
        for par in model:
            if par not in ['gc', 'len', 'distance']:
                print sys.stderr, ("Not all arguments in -v/--model are valid.")
                return 1
        modelbins = args.modelbins.split(',')
        for i in range(len(modelbins)):
            try:
                modelbins[i] = int(modelbins[i])
            except:
                print sys.stderr, ("Not all arguments in -n/--modelbins could be converted to integers.")
                return 1
        if len(model) != len(modelbins):
            print sys.stderr, ("-v/--model and -n/--modelbins not equal lengths.")
            return 1
    if args.prefix is None:
        frag_fname, data_fname, project_fname = args.output
    else:
        frag_fname = "%s.frags" % args.prefix
        data_fname = "%s.fcd" % args.prefix
        project_fname = "%s.fcp" % args.prefix
    frags = Fragment(frag_fname, mode='w', silent=args.silent)
    frags.load_fragments(args.bed)
    frags.save()
    del frags
    data = FiveCData(data_fname, 'w', silent=args.silent)
    if not args.bam is None:
        data.load_data_from_bam(frag_fname, args.bam)
    else:
        data.load_data_from_counts(frag_fname, args.count)
    data.save()
    del data
    fivec = FiveC(project_fname, 'w', silent=args.silent)
    fivec.load_data(data_fname)
    fivec.filter_fragments(mininteractions=args.minint, mindistance=args.mindist, maxdistance=args.maxdist)
    fivec.find_distance_parameters()
    precorrect = False
    if args.algorithm in ['binning', 'binning-express', 'binning-probability']:
        fivec.find_binning_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                   regions=regions, num_bins=modelbins, model=model,
                                                   usereads=args.binreads, learning_threshold=args.threshold,
                                                   max_iterations=args.biniter, precorrect=False)
        precorrect = True
    if args.algorithm in ['probability', 'binning-probability', 'probability-binning']:
        fivec.find_probability_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                        regions=regions, max_iterations=args.probiter,
                                                        minchange=args.change, learningstep=args.step,
                                                        precalculate=args.precalc, precorrect=precorrect)
    elif args.algorithm in ['express', 'binning-express', 'express-binning']:
        fivec.find_express_fragment_corrections(iterations=args.expiter, mindistance=args.mindist,
                                                maxdistance=args.maxdist, remove_distance=args.nodist,
                                                usereads=args.expreads, regions=regions,
                                                precorrect=precorrect, logged=args.logged, kr=args.kr)
    if args.algorithm in ['express-binning', 'probability-binning']:
        fivec.find_binning_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                   regions=regions, num_bins=modelbins, model=model,
                                                   usereads=args.binreads, learning_threshold=args.threshold,
                                                   max_iterations=args.biniter, precorrect=True)
    fivec.save()
