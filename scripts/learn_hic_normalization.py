#!/usr/bin/env python

import sys
import os
import optparse

try:
    from mpi4py import MPI
except:
    pass

import hifive


def main():
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    usage = "usage: [mpirun -np N_PROCS] %prog [options] <project_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive HiC project file"
    help = {
        "-l":"learning algorithm to use (probability, express, regression, regression-probability, or regression-express) [default: %default]",
        "-m":"minimum interaction distance to include in learning [default: %default]",
        "-x":"maximum interaction distance to include in learning [default: %default]",
        "-c":"comma-separated list of chromosome names to learn corrections for [default: all chromosomes]",
        "-o":"output file to save analyzed copy of HiFive HiC project to. If not specified, original file will be overwritten",
        "-q":"silence output messages [default: %default]",
        "-b":"number of iterations to run burn-in phase for with constant learning rate [default: %default]",
        "-a":"number of iterations to run annealing-in phase for with decreasing learning rate [default: %default]",
        "-s":"minimum mean change in fend correction parameter values needed to keep running past 'burnin_iterations' number of iterations during burn-in phase of probability learning algorithm [default: %default]",
        "-p":"precalculate correction values from fend means for probability algorithm",
        "-g":"learning rate, defined as percentage of gradient to use for updating parameter values in probability algorithm [default: %default]",
        "-i":"number of iterations to wait before explicitly calculating cost and updating display (zero indicates off) [default: %default]",
        "-e":"number of iterations to run express learning for [default: %default]",
        "-d":"remove the distant-dependent portion of the signal prior to learning corrections in express algorithm [default: %default]",
        "-w":"which set of reads, 'cis', 'trans', or 'all', to use for learning [default: %default]",
        "-f":"minimum number of interactions for fend filtering, if refiltering is required due to distance cutoff [default: %default]",
        "-r":"maximum number of iterations to run regression modeling for [default: %default]",
        "-t":"learning threshold cutoff for regression algorithm [default %default]",
        "-v":"comma-separated list of parameters to include in regression model (valid parameters are 'gc', 'len', 'distance', and 'mappability') [default: %default]",
        "-n":"comma-separated list of numbers of bins to partition model parameter values into [default: %default]",
    }
    if rank == 0:
        parser = optparse.OptionParser(usage=usage)
    else:
        parser = optparse.OptionParser(usage=optparse.SUPPRESS_USAGE, add_help_option=False)
        for key in help:
            help[key] = optparse.SUPPRESS_HELP
        parser.add_option("-h", "--help", help=optparse.SUPPRESS_HELP, dest="help", action="store_true")
    parser.add_option("-l", "--algorithm", dest="algorithm", default="probability", metavar="ALGORITHM",
                      type="choice", help=help["-l"], choices=["probability", "express", "regression",
                      "regression-probability", "regression-express"])
    parser.add_option("-m", "--min-distance", dest="mindist", default=0, metavar="MINDIST", type="int",
                      help=help["-m"], action="store")
    parser.add_option("-x", "--max-distance", dest="maxdist", default=0, metavar="MAXDIST", type="int",
                      help=help["-x"], action="store")
    parser.add_option("-c", "--chromosomes", dest="chroms", default="", metavar="CHROMS", type="string",
                      help=help["-c"], action="store")
    parser.add_option("-o", "--output-file", dest="outfname", default=None, metavar="OUT", type="string",
                      help=help["-o"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    group = optparse.OptionGroup(parser, "Probability Algorithm-Specific Options")
    group.add_option("-b", "--burnin-iterations", dest="burnin", default=1000, metavar="BURNIN", type="int",
                      help=help["-b"], action="store")
    group.add_option("-a", "--annealing-iterations", dest="anneal", default=1000, metavar="ANNEAL", type="int",
                      help=help["-a"], action="store")
    group.add_option("-s", "--min-change", dest="change", default=0.0001, metavar="CHANGE", type="float",
                      help=help["-s"], action="store")
    group.add_option("-p", "--precalculate", dest="precalc", default=False, metavar="PRECALC",
                      help=help["-p"], action="store_true")
    group.add_option("-g", "--learning-rate", dest="rate", default=0.3, metavar="RATE", type="float",
                      help=help["-g"], action="store")
    group.add_option("-i", "--display", dest="display", default=1, metavar="DISPLAY", type="int",
                      help=help["-i"], action="store")
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Express Algorithm-Specific Options")
    group.add_option("-e", "--express-iterations", dest="expiter", default=1000, metavar="EXPITER", type="int",
                      help=help["-e"], action="store")
    group.add_option("-d", "--remove-dist", dest="nodist", default=False, metavar="REMOVE",
                      help=help["-d"], action="store_true")
    group.add_option("-w", "--reads", dest="reads", default="cis", metavar="READS", type="choice",
                      help=help["-w"], choices=["cis", "trans", "all"])
    group.add_option("-f", "--min-interactions", dest="minint", default=1, metavar="MININT", type="int",
                      help=help["-f"], action="store")
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Regression Algorithm-Specific Options")
    group.add_option("-r", "--regression-iterations", dest="regiter", default=1000, metavar="REGITER", type="int",
                      help=help["-r"], action="store")
    group.add_option("-t", "--learning-threshold", dest="threshold", default=1.0, metavar="THRESHOLD", type="float",
                      help=help["-t"], action="store")
    group.add_option("-v", "--model", dest="model", default="gc,len,distance", metavar="MODEL", type="string",
                      help=help["-v"], action="store")
    group.add_option("-n", "--model-bins", dest="modelbins", default="20,20,15", metavar="NUMBINS", type="string",
                      help=help["-n"], action="store")
    parser.add_option_group(group)
    options, args = parser.parse_args()
    if len(args) < 1:
        if rank == 0:
            parser.error('incorrect number of arguments')
        else:
            sys.exit(1)
    if not os.path.exists(args[0]):
        if rank == 0:
            parser.error('could not find %s' % args[0])
        else:
            sys.exit(1)
    chroms = options.chroms.split(',')
    hic = hifive.HiC(args[0], 'r', silent=options.silent)
    precorrect = False
    if options.algorithm in ['regression', 'regression-express', 'regression-probability']:
        model = options.model.split(',')
        modelbins = options.modelbins.split(',')
        if len(model) != len(modelbins):
            if rank == 0:
                parser.error('model and model-bins must be of the same length')
            else:
                sys.exit(1)
        for i in range(len(modelbins)):
            modelbins[i] = int(modelbins[i])
        hic.find_regression_fend_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                             chroms=chroms, num_bins=modelbins, model=model,
                                             learning_threshold=options.threshold, max_iterations=options.regiter)
        precorrect = True
    if options.algorithm in ['probability', 'regression-probability']:
        hic.find_probability_fend_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                              minchange=options.change, burnin_iterations=options.burnin,
                                              annealing_iterations=options.anneal, learningrate=options.rate,
                                              display=options.display, chroms=chroms,
                                              precalculate=options.precalc, precorrect=precorrect)
    elif options.algorithm in ['express', 'regression-express']:
        hic.find_express_fend_corrections(iterations=options.expiter, mindistance=options.mindist,
                                          maxdistance=options.maxdist, remove_distance=options.nodist,
                                          usereads=options.reads, mininteractions=options.minint,
                                          chroms=chroms, precorrect=precorrect)
    if rank == 0:
        hic.save(options.outfname)


if __name__ == "__main__":
    main()
