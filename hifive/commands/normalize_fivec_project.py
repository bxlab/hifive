#!/usr/bin/env python

import sys

from ..fivec import FiveC


def run(args):
    regions = args.regions.split(',')
    if len(regions) == 1 and regions[0] == '':
        regions = []
    for i in range(len(regions)):
        try:
            regions[i] = int(regions[i])
        except:
            print sys.stderr, ("Not all arguments in -r/--regions could be converted to integers.")
            return 1
    if not args.model is None:
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
    fivec = FiveC(args.project, 'r', silent=args.silent)
    precorrect = False
    if args.algorithm in ['regression', 'regression-express', 'regression-probability']:
        fivec.find_regression_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                   regions=regions, num_bins=modelbins, model=model,
                                                   usereads=args.regreads, learning_threshold=args.threshold,
                                                   max_iterations=args.regiter)
        precorrect = True
    if args.algorithm in ['probability', 'regression-probability']:
        fivec.find_probability_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                    regions=regions, burnin_iterations=args.burnin,
                                                    annealing_iterations=args.anneal, learningrate=args.rate,
                                                    precalculate=args.precalc, precorrect=precorrect)
    elif args.algorithm in ['express', 'regression-express']:
        fivec.find_express_fragment_corrections(iterations=args.expiter, mindistance=args.mindist,
                                                maxdistance=args.maxdist, remove_distance=args.nodist,
                                                usereads=args.expreads, regions=regions,
                                                precorrect=precorrect)
    fivec.save(args.outfname)
