#!/usr/bin/env python

from ..fivec import FiveC


def run(args):
    fivec = FiveC(args.output, 'w', silent=args.silent)
    fivec.load_data(args.data)
    fivec.filter_fragments(mininteractions=args.minint, mindistance=args.mindist, maxdistance=args.maxdist)
    fivec.find_distance_parameters()
    fivec.save()
