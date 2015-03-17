#!/usr/bin/env python

from ..fend import Fend


def run(args):
    fends = Fend(args.output, mode='w', silent=args.silent)
    if args.bed is None:
        fends.load_fends(args.fend, genome_name=args.genome, re_name=args.re, format="fend")
    else:
        fends.load_fends(args.bed, genome_name=args.genome, re_name=args.re, format="bed")
    fends.save()
