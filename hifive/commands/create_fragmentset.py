#!/usr/bin/env python

from ..fragment import Fragment


def run(args):
    fragments = Fragment(args.output, mode='w', silent=args.silent)
    fragments.load_fragments(args.bed, genome_name=args.genome, re_name=args.re)
    fragments.save()
    return None
