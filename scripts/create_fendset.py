#!/usr/bin/env python

import sys
import os
import optparse

import hifive


def main():
    usage = "usage: %prog [options] <out_file>\n\nArguments:"
    usage += "\n<out_file>  destination for fend file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", "--fend", dest="fend", default=None, metavar="FEND", type="string",
                      help="fend file in HiCPipe-compatible MAT format", action="store")
    parser.add_option("-b", "--bed", dest="bed", default=None, metavar="BED", type="string",
                      help="restriction enzyme cut sites or RE fragment boundaries in BED format", action="store")
    parser.add_option("-g", "--genome", dest="genome", default="", metavar="GENOME", type="string",
                      help="name of genome", action="store")
    parser.add_option("-r", "--re", dest="re", default="", metavar="RE", type="string",
                      help="name of restriction enzyme", action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help="silence output messages [default: %default]")
    options, args = parser.parse_args()
    if len(args) < 1:
        parser.error('incorrect number of arguments')
    if ((options.bed is None and options.fend is None) or
        (not options.bed is None and not options.fend is None)):
        parser.error('requires exactly one fend file type (-f/--fend or -b/--bed)')
    fends = hifive.Fend(args[0], mode='w', silent=options.silent)
    if options.bed is None:
        fends.load_fends(options.fend, genome_name=options.genome, re_name=options.re, format="fend")
    else:
        fends.load_fends(options.bed, genome_name=options.genome, re_name=options.re, format="bed")
    fends.fends.close()


if __name__ == "__main__":
    main()
