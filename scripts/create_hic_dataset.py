#!/usr/bin/env python

import sys
import os
import optparse

import hifive


def main():
    usage = "usage: %prog [options] <fend_file> <out_file>\n\nArguments:"
    usage += "\n<fend_file>    HiFive fend file"
    usage += "\n<out_file>     destination for HiC data file"
    help = {
      "bam":"a pair of comma-separated bam read end file lists (first and second end) for a single sequencing run. For multiple runs, this option can be passed multiple times",
      "mat":"a HiCPipe-style mat file containing fend pair counts",
      "raw":"a tab-separated text file containing pairs of read ends (chr1 pos1 strand1 chr2 pos2 strand2), one per line. For multiple files, this option can be passed multiple times",
      "insert":"maximum allowable distance sum between both fend ends and cutsites [default: %default]",
      "silent":"silence output messages [default: %default]",
    }
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bam", dest="bam", default=[], metavar="BAM1 BAM2", type="string", nargs=2,
                      help=help["bam"],
                      action="append")
    parser.add_option("-m", "--mat", dest="mat", metavar="MAT", type="string", default=None,
                      help=help["mat"], action="store")
    parser.add_option("-r", "--raw", dest="raw", metavar="RAW", type="string", default=[],
                      help=help["raw"], action="append")
    parser.add_option("-i", "--max-insert", dest="insert", default=1000, metavar="INSERT", type="int",
                      help=help["insert"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["silent"])
    options, args = parser.parse_args()
    if len(args) < 2:
          parser.error('incorrect number of arguments')
    if options.mat is None and len(options.raw) == 0 and len(options.bam) > 0:
        data = hifive.HiCData(args[1], 'w', silent=options.silent)
        data.load_data_from_bam(args[0], options.bam, options.insert)
        data.save()
    elif not options.mat is None and len(options.raw) == 0 and len(options.bam) == 0:
        data = hifive.HiCData(args[1], 'w', silent=options.silent)
        data.load_data_from_mat(args[0], options.mat, options.insert)
        data.save()
    elif options.mat is None and len(options.raw) > 0 and len(options.bam) == 0:
        data = hifive.HiCData(args[1], 'w', silent=options.silent)
        data.load_data_from_raw(args[0], options.raw, options.insert)
        data.save()
    else:
        parser.error("data must be provided by exactly one of the following options: -b/--bam, -m/--mat, or '-r/--raw")


if __name__ == "__main__":
    main()
