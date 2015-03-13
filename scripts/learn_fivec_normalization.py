#!/usr/bin/env python

import sys
import os
import optparse

import hifive


def main():
    usage = "usage: %prog [options] <project_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive 5C project file"
    help = {
        "-l":"learning algorithm to use (probability, express, regression, regression-probability, or regression-express) [default: %default]",
        "-m":"minimum interaction distance to include in learning [default: %default]",
        "-x":"maximum interaction distance to include in learning [default: %default]",
        "-r":"comma-separated list of region indices to learn corrections for [default: all regions]",
        "-o":"output file to save analyzed copy of HiFive 5C project to. If not specified, original file will be overwritten",
        "-q":"silence output messages [default: %default]",
        "-b":"number of iterations to run burn-in phase for with constant learning rate [default: %default]",
        "-a":"number of iterations to run annealing-in phase for with decreasing learning rate [default: %default]",
        "-p":"precalculate correction values from fend means for probability algorithm",
        "-g":"learning rate, defined as percentage of gradient to use for updating parameter values in probability algorithm [default: %default]",
        "-e":"number of iterations to run express learning for [default: %default]",
        "-d":"remove the distant-dependent portion of the signal prior to learning corrections in express algorithm [default: %default]",
        "-w":"which set of reads, 'cis', 'trans', or 'all', to use for express learning [default: %default]",
        "-i":"maximum number of iterations to run regression modeling for [default: %default]",
        "-t":"learning threshold cutoff for regression algorithm [default %default]",
        "-y":"which set of reads, 'cis', 'trans', or 'all', to use for regression learning [default: %default]",
        "-v":"comma-separated list of parameters to include in regression model (valid parameters are 'gc', 'len', and 'distance') [default: %default]",
        "-n":"comma-separated list of numbers of bins to partition model parameter values into [default: %default]",
    }
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-l", "--algorithm", dest="algorithm", default="probability", metavar="ALGORITHM",
                      type="choice", help=help["-l"], choices=["probability", "express", "regression",
                      "regression-probability", "regression-express"])
    parser.add_option("-m", "--min-distance", dest="mindist", default=0, metavar="MINDIST", type="int",
                      help=help["-m"], action="store")
    parser.add_option("-x", "--max-distance", dest="maxdist", default=0, metavar="MAXDIST", type="int",
                      help=help["-x"], action="store")
    parser.add_option("-r", "--regions", dest="regions", default="", metavar="REGIONS", type="string",
                      help=help["-r"], action="store")
    parser.add_option("-o", "--output-file", dest="outfname", default=None, metavar="OUT", type="string",
                      help=help["-o"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    group = optparse.OptionGroup(parser, "Probability Algorithm-Specific Options")
    group.add_option("-b", "--burnin-iterations", dest="burnin", default=1000, metavar="BURNIN", type="int",
                      help=help["-b"], action="store")
    group.add_option("-a", "--annealing-iterations", dest="anneal", default=1000, metavar="ANNEAL", type="int",
                      help=help["-a"], action="store")
    group.add_option("-p", "--precalculate", dest="precalc", default=False, metavar="PRECALC",
                      help=help["-p"], action="store_true")
    group.add_option("-g", "--learning-rate", dest="rate", default=0.01, metavar="RATE", type="float",
                      help=help["-g"], action="store")
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Express Algorithm-Specific Options")
    group.add_option("-e", "--express-iterations", dest="expiter", default=1000, metavar="EXPITER", type="int",
                      help=help["-e"], action="store")
    group.add_option("-d", "--remove-dist", dest="nodist", default=False, metavar="REMOVE",
                      help=help["-d"], action="store_true")
    group.add_option("-w", "--express-reads", dest="expreads", default="cis", metavar="READS", type="choice",
                      help=help["-w"], choices=["cis", "trans", "all"])
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Regression Algorithm-Specific Options")
    group.add_option("-i", "--regression-iterations", dest="regiter", default=1000, metavar="REGITER", type="int",
                      help=help["-r"], action="store")
    group.add_option("-t", "--learning-threshold", dest="threshold", default=1.0, metavar="THRESHOLD", type="float",
                      help=help["-t"], action="store")
    group.add_option("-y", "--regression-reads", dest="regreads", default="cis", metavar="READS", type="choice",
                      help=help["-y"], choices=["cis", "trans", "all"])
    group.add_option("-v", "--model", dest="model", default="gc,len,distance", metavar="MODEL", type="string",
                      help=help["-v"], action="store")
    group.add_option("-n", "--model-bins", dest="modelbins", default="20,20,15", metavar="NUMBINS", type="string",
                      help=help["-n"], action="store")
    parser.add_option_group(group)
    options, args = parser.parse_args()
    if len(args) < 1:
        parser.error('incorrect number of arguments')
    if not os.path.exists(args[0]):
        parser.error('could not find %s' % args[0])
    regions = options.regions.split(',')
    if len(regions) == 1 and regions[0] == '':
        regions = []
    for i in range(len(regions)):
        regions[i] = int(regions[i])
    fivec = hifive.FiveC(args[0], 'r', silent=options.silent)
    precorrect = False
    if options.algorithm in ['regression', 'regression-express', 'regression-probability']:
        model = options.model.split(',')
        modelbins = options.modelbins.split(',')
        if len(model) != len(modelbins):
            parser.error('model and model-bins must be of the same length')
        for i in range(len(modelbins)):
            modelbins[i] = int(modelbins[i])
        fivec.find_regression_fragment_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                                   regions=regions, num_bins=modelbins, model=model,
                                                   usereads=options.regreads, learning_threshold=options.threshold,
                                                   max_iterations=options.regiter)
        precorrect = True
    if options.algorithm in ['probability', 'regression-probability']:
        fivec.find_probability_fragment_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                                    regions=regions, burnin_iterations=options.burnin,
                                                    annealing_iterations=options.anneal, learningrate=options.rate,
                                                    precalculate=options.precalc, precorrect=precorrect)
    elif options.algorithm in ['express', 'regression-express']:
        fivec.find_express_fragment_corrections(iterations=options.expiter, mindistance=options.mindist,
                                                maxdistance=options.maxdist, remove_distance=options.nodist,
                                                usereads=options.expreads, regions=regions,
                                                precorrect=precorrect)
    fivec.save(options.outfname)


if __name__ == "__main__":
    main()
