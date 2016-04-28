#!/usr/bin/env python

import sys

from ..fend import Fend


def run(args):
    if args.binned == 0 and args.bed is None:
        print  >> sys.stderr, ("Non-uniforming binning (binned=0) must have a bed file to read bin partitions from.\n"),
        return None
    elif (args.binned is None or args.binned < 1) and args.length is not None:
        print  >> sys.stderr, ("Binning from a chromosome length file needs a positive integer value for binning.\n"),
        return None
    fends = Fend(args.output, mode='w', binned=args.binned, silent=args.silent)
    if args.bed is not None:
    	if args.binned is not None and args.binned == 0:
    		fends.load_bins(args.bed, genome_name=args.genome, format='bed')
    	else:
	        fends.load_fends(args.bed, genome_name=args.genome, re_name=args.re, format="bed")
    elif args.fend is not None:
        fends.load_fends(args.fend, genome_name=args.genome, re_name=args.re, format="fend")
    else:
		fends.load_bins(args.length, genome_name=args.genome, format='len')
    fends.save()
