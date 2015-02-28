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
    usage = "usage: %prog [options] <project_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive HiC project file"
    help = {
        "-c":"number of iterations to run learning algorithm for [default: %default]",
        "-m":"minimum number of interactions for fend filtering, if refiltering is required due to distance cutoff [default: %default]",
        "-x":"maximum interaction distance to include for learning (a zero indicates no maximum) [default: %default]",
        "-r":"which set of reads, 'cis', 'trans', or 'both', to use for learning [default: %default]",
        "-i":"number of iterations to run learning algorithm for [default: %default]",
        "-d":"remove the distant-dependent portion of the signal prior to learning corrections [default: %default]",
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
    parser.add_option("-i", "--iterations", dest="iter", default=1000, metavar="ITER", type="int",
                      help=help["-i"], action="store")
    parser.add_option("-c", "--min-interactions", dest="min_int", default=0, metavar="MININT", type="int",
                      help=help["-c"], action="store")
    parser.add_option("-m", "--min-distance", dest="min_dist", default=0, metavar="MINDIST", type="int",
                      help=help["-m"], action="store")
    parser.add_option("-x", "--max-distance", dest="max_dist", default=0, metavar="MAXDIST", type="int",
                      help=help["-x"], action="store")
    parser.add_option("-r", "--reads", dest="reads", default="cis", metavar="READS", type="string",
                      help=help["-r"], action="store")
    parser.add_option("-d", "--remove-distance", dest="remove_dist", default=False, metavar="NODIST",
                      help=help["-d"], action="store_true")
    parser.add_option("-o", "--output-file", dest="output", default=None, metavar="OUT", type="string",
                      help=help["-o"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    options, args = parser.parse_args()
    if len(args) < 1:
        if rank == 0:
            parser.error('incorrect number of arguments')
        else:
            sys.exit(1)
    if options.reads not in ['cis', 'trans', 'both']:
        if rank == 0:
            parser.error("-r/--reads accepts 'cis', 'trans', or 'both'")
        else:
            sys.exit(1)
    hic = hifive.HiC(args[0], 'r', silent=options.silent)
    hic.find_express_fend_corrections(iterations=options.iter, mindistance=options.min_dist,
                                      maxdistance=options.max_dist, mininteractions=options.min_int,
                                      usereads=options.reads, remove_distance=options.remove_dist)
    if rank == 0:
        hic.save(options.output)


if __name__ == "__main__":
    main()
