#!/usr/bin/env python

import sys

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
        modelbins = args.modelbins.split(',')
        parameters = args.parameters.split(',')
        for i in range(len(modelbins)):
            try:
                modelbins[i] = int(modelbins[i])
            except:
                print sys.stderr, ("Not all arguments in -n/--modelbins could be converted to integers.")
                return 1
        if len(model) != len(modelbins) or len(model) != len(parameters):
            print sys.stderr, ("-v/--model, -n/--modelbins, and -u/--parameter-types are not equal lengths.")
            return 1
    fivec = FiveC(args.project, 'r', silent=args.silent)
    if args.algorithm.split('-')[0] == 'binning':
        fivec.find_binning_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                regions=regions, num_bins=modelbins, model=model,
                                                usereads=args.binreads, parameters=parameters,
                                                learning_threshold=args.threshold,
                                                max_iterations=args.biniter, precorrect=False)
    if args.algorithm.split('-')[0] == 'probability':
        fivec.find_probability_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                    regions=regions, max_iterations=args.probiter,
                                                    minchange=args.change, learningstep=args.step,
                                                    precalculate=args.precalc, precorrect=False)
    if args.algorithm.split('-')[0] == 'express':
        fivec.find_express_fragment_corrections(iterations=args.expiter, mindistance=args.mindist,
                                                maxdistance=args.maxdist, remove_distance=args.nodist,
                                                usereads=args.expreads, regions=regions,
                                                precorrect=False, logged=args.logged, kr=args.kr)
    if len(args.algorithm.split('-')) > 1:
        if args.algorithm.split('-')[1] == 'binning':
            fivec.find_binning_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                    regions=regions, num_bins=modelbins, model=model,
                                                    usereads=args.binreads, parameters=parameters,
                                                    learning_threshold=args.threshold,
                                                    max_iterations=args.biniter, precorrect=True)
        if args.algorithm.split('-')[1] == 'probability':
            fivec.find_probability_fragment_corrections(mindistance=args.mindist, maxdistance=args.maxdist,
                                                        regions=regions, max_iterations=args.probiter,
                                                        minchange=args.change, learningstep=args.step,
                                                        precalculate=args.precalc, precorrect=True)
        if args.algorithm.split('-')[1] == 'express':
            fivec.find_express_fragment_corrections(iterations=args.expiter, mindistance=args.mindist,
                                                    maxdistance=args.maxdist, remove_distance=args.nodist,
                                                    usereads=args.expreads, regions=regions,
                                                    precorrect=True, logged=args.logged, kr=args.kr)
    fivec.save(args.output)
