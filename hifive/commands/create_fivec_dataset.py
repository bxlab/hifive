#!/usr/bin/env python

from ..fivec_data import FiveCData


def run(args):
    data = FiveCData(args.output, 'w', silent=args.silent)
    if args.bam is None:
        data.load_data_from_counts(args.fragment, args.count)
    else:
        data.load_data_from_bam(args.fragment, args.bam)
    data.save()
