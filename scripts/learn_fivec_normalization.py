#!/usr/bin/env python

import sys
import optparse

import hifive


def main():
    usage = "usage: %prog [options] <project_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive 5C project file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-l", "--learning-rate", dest="rate", default=0.1, metavar="RATE", type="float",
                      help="learning rate, defined as percentage of gradient to use for updating parameter values [default: %default]",
                      action="store")
    parser.add_option("-b", "--burnin-iterations", dest="burnin", default=1000, metavar="BURNIN", type="int",
                      help="number of iterations to run burn-in phase for with constant learning rate [default: %default]",
                      action="store")
    parser.add_option("-a", "--annealing-iterations", dest="anneal", default=1000, metavar="ANNEAL", type="int",
                      help="number of iterations to run annealing-in phase for with decreasing learning rate [default: %default]",
                      action="store")
    parser.add_option("-m", "--min-distance", dest="min_dist", default=0, metavar="MINDIST", type="int",
                      help="minimum interaction distance to include in learning [default: %default]",
                      action="store")
    parser.add_option("-x", "--max-distance", dest="max_dist", default=0, metavar="MAXDIST", type="int",
                      help="maximum interaction distance to include in learning (zero indicates no maximum) [default: %default]",
                      action="store")
    parser.add_option("-p", "--precalculate", dest="precalc", default=False, metavar="PRECALC",
                      help="indicates corrections should initially be set to fragment means [default: %default]",
                      action="store_true")
    parser.add_option("-o", "--output-file", dest="output", default=None, metavar="OUT", type="string",
                      help="output file to save analyzed copy of HiFive 5C project to. If not specified, original file will be overwritten",
                      action="store")
    parser.add_option("-d", "--display", dest="display", default=1, metavar="DISPLAY", type="int",
                      help="number of iterations to wait before explicitly calculating cost and updating display (zero indicates off) [default: %default]",
                      action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help="silence output messages [default: %default]")
    options, args = parser.parse_args()
    if len(args) < 1:
        parser.error('incorrect number of arguments')
    fivec = hifive.FiveC(args[0], 'r', silent=options.silent)
    fivec.find_fragment_corrections(display=options.display, mindistance=options.min_dist,
                                    maxdistance=options.max_dist, burnin_iterations=options.burnin,
                                    annealing_iterations=options.anneal, precalculate=options.precalc,
                                    learningrate=options.rate)
    fivec.save(options.output)


if __name__ == "__main__":
    main()
