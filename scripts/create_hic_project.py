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
    usage = "usage: [mpirun -np N_PROCS] %prog [options] <data_file> <out_file>\n\nArguments:"
    usage += "\n<data_file>  HiFive HiC data file"
    usage += "\n<out_file>   destination for HiC project file"
    help = {
        "-i":"minimum number of interactions needed for valid fends [default: %default]",
        "-d":"smallest interaction distance to be included for filtering fends [default: %default]",
        "-m":"largest interaction distance to be included for filtering fends (zero indicates no maximum) [default: %default]",
        "-s":"smallest interaction distance bin size for distance function [default: %default]",
        "-b":"number of bins to partion interaction distance range into for distance function [default: %default]",
        "-q":"silence output messages [default: %default]",
    }
    if rank == 0:
        parser = optparse.OptionParser(usage=usage)
    else:
        parser = optparse.OptionParser(usage=optparse.SUPPRESS_USAGE, add_help_option=False)
        for key in help:
            help[key] = optparse.SUPPRESS_HELP
        parser.add_option("-h", "--help", help=optparse.SUPPRESS_HELP, dest="help", action="store_true")
    parser.add_option("-i", "--min-interactions", dest="min_int", default=10, metavar="MININT", type="int",
                      help=help["-i"], action="store")
    parser.add_option("-d", "--min-distance", dest="min_dist", default=0, metavar="MINDIST", type="int",
                      help=help["-d"], action="store")
    parser.add_option("-m", "--max-distance", dest="max_dist", default=0, metavar="MAXDIST", type="int",
                      help=help["-m"], action="store")
    parser.add_option("-s", "--min-size", dest="min_size", default=1000, metavar="MINBIN", type="int",
                      help=help["-s"], action="store")
    parser.add_option("-b", "--num-bins", dest="num_bins", default=1000, metavar="NUMBINS", type="int",
                      help=help["-b"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    options, args = parser.parse_args()
    if len(args) < 2:
        if rank == 0:
            parser.error('incorrect number of arguments')
        else:
            sys.exit(1)
    if rank == 0:
        hic = hifive.HiC(args[1], 'w', silent=options.silent)
        hic.load_data(args[0])
        hic.filter_fends(mininteractions=options.min_int, mindistance=options.min_dist, maxdistance=options.max_dist)
        hic.save()
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
    else:
        comm.recv(source=0, tag=11)
        hic = hifive.HiC(args[1], 'r', silent=True)
    hic.find_distance_parameters(minsize=options.min_size, numbins=options.num_bins)
    if rank == 0:
        hic.save()
    return None

    
if __name__ == "__main__":
    main()
