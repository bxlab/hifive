#!/usr/bin/env python

import sys
import optparse

import hifive


def main():
    usage = "usage: %prog [options] <project_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive 5C project file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-m", "--min-distance", dest="min_dist", default=0, metavar="MINDIST", type="int",
                      help="minimum interaction distance to include in learning [default: %default]",
                      action="store")
    parser.add_option("-x", "--max-distance", dest="max_dist", default=0, metavar="MAXDIST", type="int",
                      help="maximum interaction distance to include in learning (zero indicates no maximum) [default: %default]",
                      action="store")
    parser.add_option("-i", "--iterations", dest="iter", default=1000, metavar="ITER", type="int",
                      help="number of iterations to run learning algorithm for [default: %default]",
                      action="store")
    parser.add_option("-d", "--remove-distance", dest="remove_dist", default=False, metavar="NODIST",
                      help="remove the distant-dependent portion of the signal prior to learning corrections [default: %default]",
                      action="store_true")
    parser.add_option("-o", "--output-file", dest="output", default=None, metavar="OUT", type="string",
                      help="output file to save analyzed copy of HiFive 5C project to. If not specified, original file will be overwritten",
                      action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help="silence output messages [default: %default]")
    options, args = parser.parse_args()
    if len(args) < 1:
        parser.error('incorrect number of arguments')
    fivec = hifive.FiveC(args[0], 'r', silent=options.silent)
    fivec.find_express_fragment_corrections(iterations=options.iter, remove_distance=options.remove_dist,
                                            mindistance=options.min_dist, maxdistance=options.max_dist)
    fivec.save(options.output)


if __name__ == "__main__":
    main()
