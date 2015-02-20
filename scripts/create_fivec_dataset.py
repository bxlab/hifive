#!/usr/bin/env python

import sys
import os
import optparse

import hifive


def main():
    usage = "usage: %prog [options] <frag_file> <out_file>\n\nArguments:"
    usage += "\n<frag_file>  HiFive fragment file"
    usage += "\n<out_file>   destination for 5C data file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bam", dest="bam", default=[], metavar="BAM1 BAM2", type="string", nargs=2,
                      help="a pair of bam read end files for a single sequencing run. For multiple runs, this option can be passed multiple times",
                      action="append")
    parser.add_option("-c", "--count", dest="count", metavar="CNT", type="string", default=[],
                      help="a tab-separated text file containing a pair of fragment names and the number of observed reads for that pair (fragment1 fragment2 count), one per line. For multiple files, this option can be passed multiple times",
                      action="append")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help="silence output messages [default: %default]")
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error('incorrect number of arguments')
    if len(options.bam) == 0 and len(options.count) > 0:
        data = hifive.FiveCData(args[1], 'w', silent=options.silent)
        data.load_data_from_counts(args[0], options.count)
        data.save()
    elif len(options.bam) > 0 and len(options.count) == 0:
        data = hifive.FiveCData(args[1], 'w', silent=options.silent)
        data.load_data_from_bam(args[0], options.bam)
        data.save()
    else:
        parser.error("data must be provided by exactly one of the following options: -c/--count or -b/--bam")


if __name__ == "__main__":
    main()
