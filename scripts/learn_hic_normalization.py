#!/usr/bin/env python

import sys
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
        "-b":"number of iterations to run burn-in phase for with constant learning rate [default: %default]",
        "-a":"number of iterations to run annealing-in phase for with decreasing learning rate [default: %default]",
        "-s":"minimum mean change in fend correction parameter values needed to keep running past 'burnin_iterations' number of iterations during burn-in phase [default: %default]",
        "-p":"precalculate correction values from fend means",
        "-m":"minimum interaction distance to include in learning [default: %default]",
        "-x":"maximum interaction distance to include in learning [default: %default]",
        "-l":"learning rate, defined as percentage of gradient to use for updating parameter values [default: %default]",
        "-d":"number of iterations to wait before explicitly calculating cost and updating display (zero indicates off) [default: %default]",
        "-c":"comma-separated list of chromosome names to learn corrections for [default: all chromosomes]",
        "-o":"output file to save analyzed copy of HiFive HiC project to. If not specified, original file will be overwritten",
        "-q":"silence output messages [default: %default]",
    }
    if rank == 0:
        parser = optparse.OptionParser(usage=usage)
    else:
        parser = optparse.OptionParser(usage=optparse.SUPPRESS_USAGE, add_help_option=False)
        for key in help:
            help[key] = optparse.SUPPRESS_HELP
        parser.add_option("-h", "--help", help=optparse.SUPPRESS_HELP, dest="help", action="store_true")
    parser.add_option("-l", "--learning-rate", dest="rate", default=0.1, metavar="RATE", type="float",
                      help=help["-l"], action="store")
    parser.add_option("-b", "--burnin-iterations", dest="burnin", default=1000, metavar="BURNIN", type="int",
                      help=help["-b"], action="store")
    parser.add_option("-a", "--annealing-iterations", dest="anneal", default=1000, metavar="ANNEAL", type="int",
                      help=help["-a"], action="store")
    parser.add_option("-m", "--min-distance", dest="min_dist", default=0, metavar="MINDIST", type="int",
                      help=help["-m"], action="store")
    parser.add_option("-x", "--max-distance", dest="max_dist", default=0, metavar="MAXDIST", type="int",
                      help=help["-x"], action="store")
    parser.add_option("-s", "--min-change", dest="change", default=0.00001, metavar="CHANGE", type="float",
                      help=help["-s"], action="store")
    parser.add_option("-p", "--precalculate", dest="precalc", default=False, metavar="PRECALC",
                      help=help["-p"], action="store_true")
    parser.add_option("-d", "--display", dest="display", default=1, metavar="DISPLAY", type="int",
                      help=help["-d"], action="store")
    parser.add_option("-c", "--chromosomes", dest="chroms", default="", metavar="CHROMS", type="string",
                      help=help["-c"], action="store")
    parser.add_option("-o", "--output-file", dest="outfname", default=None, metavar="OUT", type="string",
                      help=help["-o"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    options, args = parser.parse_args()
    if len(args) < 1:
        if rank == 0:
            parser.error('incorrect number of arguments')
        else:
            sys.exit(1)
    hic = hifive.HiC(args[0], 'r', silent=options.silent)
    hic.find_fend_corrections(display=options.display, mindistance=options.min_dist, maxdistance=options.max_dist,
                              chroms=options.chroms.split(','), learningrate=options.rate,
                              burnin_iterations=options.burnin, annealing_iterations=options.anneal,
                              minchange=options.change, precalculate=options.precalc)
    if rank == 0:
        hic.save(options.outfname)


if __name__ == "__main__":
    main()
