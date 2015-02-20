#!/usr/bin/env python

import sys
import os
import optparse

import hifive


def main():
    usage = "usage: %prog [options] <bed_file> <out_file>\n\nArguments:"
    usage += "\n<bed_file>  bed file containined fragment bounds"
    usage += "\n<out_file>  destination for fragment file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-g", "--genome", dest="genome", default="", metavar="GENOME", type="string",
                      help="name of genome", action="store")
    parser.add_option("-r", "--re", dest="re", default="", metavar="RE", type="string",
                      help="name of restriction enzyme", action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help="silence output messages [default: %default]")
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error('incorrect number of arguments')
    fragments = hifive.Fragment(args[1], mode='w', silent=options.silent)
    fragments.load_fragments(args[0], genome_name=options.genome, re_name=options.re)
    fragments.fragments.close()
    return None


if __name__ == "__main__":
    main()
