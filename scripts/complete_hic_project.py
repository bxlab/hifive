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
    parser = setup_parser()
    options, args = parser.parse_args()
    if rank > 0:
        options.silent = True
    check_options(parser, options, args)
    if not options.silent:
        display_options(options, args)
    fend_fname = args[0]
    data_fname = args[1]
    hic_fname = args[2]
    if rank == 0:
        fends = hifive.Fend(fend_fname, mode='w', silent=options.silent)
        if not options.fend is None:
            fends.load_fends(options.fend, format='fend')
            format = 'fend'
        else:
            fends.load_fends(options.bed, format='bed')
        fends.fends.close()
        del fends
        data = hifive.HiCData(data_fname, 'w', silent=options.silent)
        if len(options.bam) > 0:
            data.load_data_from_bam(fend_fname, options.bam, options.insert)
        elif len(options.raw) > 0:
            data.load_data_from_raw(fend_fname, options.raw, options.insert)
        else:
            data.load_data_from_mat(fend_fname, options.mat, options.insert)
        data.save()
        del data
        for i in range(1, num_procs):
            comm.send(1, dest=i, tag=11)
    else:
        comm.recv(source=0, tag=11)
    hic = hifive.HiC(hic_fname, 'w', silent=options.silent)
    hic.load_data(data_fname)
    hic.filter_fends(mininteractions=options.minint, mindistance=options.mindist, maxdistance=options.maxdist)
    hic.find_distance_parameters(minsize=options.minbin, numbins=options.distbins)
    precorrect = False
    if options.algorithm in ['regression', 'regression-express', 'regression-probability']:
        hic.find_regression_fend_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                             chroms=options.chroms, num_bins=options.modelbins, model=options.model,
                                             learning_threshold=options.threshold, max_iterations=options.regiter,
                                             usereads=options.regreads)
        precorrect = True
    if options.algorithm in ['probability', 'regression-probability']:
        hic.find_probability_fend_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                              minchange=options.change, burnin_iterations=options.burnin,
                                              annealing_iterations=options.anneal, learningrate=options.rate,
                                              display=options.display, chroms=options.chroms,
                                              precalculate=options.precalc, precorrect=precorrect)
    elif options.algorithm in ['express', 'regression-express']:
        hic.find_express_fend_corrections(iterations=options.expiter, mindistance=options.mindist,
                                          maxdistance=options.maxdist, remove_distance=options.nodist,
                                          usereads=options.expreads, mininteractions=options.minint,
                                          chroms=options.chroms, precorrect=precorrect)
    if rank == 0:
        hic.save()
    return 0

def setup_parser():
    usage = "usage: [mpirun -np N_PROCS] %prog [options] <out_fend> <out_data> <out_project>\n\nArguments:"
    usage += "\n<out_fend>      a file to write the HiFive HiC fend file to"
    usage += "\n<out_data>      a file to write the HiFive HiC data file to"
    usage += "\n<out_project>   a file to write the HiFive HiC project file to"
    help = {
        "-j":"maximum allowable distance sum between both fend ends and cutsites [default: %default]",
        "-f":"minimum number of interactions for fend filtering [default: %default]",
        "-k":"smallest interaction distance bin size for distance function [default: %default]",
        "-u":"number of bins to partion interaction distance range into for distance function [default: %default]",
        "-m":"minimum interaction distance to include in filtering and learning [default: %default]",
        "-x":"maximum interaction distance to include in filtering and learning [default: %default]",
        "-l":"learning algorithm to use (probability, express, regression, regression-probability, or regression-express) [default: %default]",
        "-c":"comma-separated list of chromosome names to learn corrections for [default: all chromosomes]",
        "-q":"silence output messages [default: %default]",
        "-F":"fend file in HiCPipe-compatible MAT format (one-indexed)",
        "-B":"restriction enzyme cut sites or RE fragment boundaries in BED format",
        "-S":"a pair of comma-separated bam read end file lists (first and second end) for a single sequencing run. For multiple runs, this option can be passed multiple times",
        "-M":"a HiCPipe-style MAT file containing fend pair counts",
        "-R":"a tab-separated text file containing pairs of read ends (chr1 pos1 strand1 chr2 pos2 strand2), one per line. For multiple files, this option can be passed multiple times",
        "-b":"number of iterations to run burn-in phase for with constant learning rate [default: %default]",
        "-a":"number of iterations to run annealing-in phase for with decreasing learning rate [default: %default]",
        "-s":"minimum mean change in fend correction parameter values needed to keep running past 'burnin_iterations' number of iterations during burn-in phase of probability learning algorithm [default: %default]",
        "-p":"precalculate correction values from fend means for probability algorithm",
        "-g":"learning rate, defined as percentage of gradient to use for updating parameter values in probability algorithm [default: %default]",
        "-i":"number of iterations to wait before explicitly calculating cost and updating display (zero indicates off) [default: %default]",
        "-e":"number of iterations to run express learning for [default: %default]",
        "-d":"remove the distant-dependent portion of the signal prior to learning corrections in express algorithm [default: %default]",
        "-w":"which set of reads, 'cis', 'trans', or 'all', to use for express learning [default: %default]",
        "-r":"maximum number of iterations to run regression modeling for [default: %default]",
        "-t":"learning threshold cutoff for regression algorithm [default %default]",
        "-y":"which set of reads, 'cis', 'trans', or 'all', to use for regression learning [default: %default]",
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
    parser.add_option("-j", "--insert", dest="insert", default=500, metavar="INSERT", type="int", action="store",
                      help=help["-j"])
    parser.add_option("-f", "--min-interactions", dest="minint", default=10, metavar="MININT", type="int",
                      help=help["-f"], action="store")
    parser.add_option("-k", "--min-binsize", dest="minbin", default=1000, metavar="MINBIN", type="int",
                      action="store", help=help["-k"])
    parser.add_option("-u", "--distance-bins", dest="distbins", default=100, metavar="DISTBINS", type="int",
                      action="store", help=help["-u"])
    parser.add_option("-m", "--min-distance", dest="mindist", default=0, metavar="MINDIST", type="int",
                      help=help["-m"], action="store")
    parser.add_option("-x", "--max-distance", dest="maxdist", default=0, metavar="MAXDIST", type="int",
                      help=help["-x"], action="store")
    parser.add_option("-l", "--algorithm", dest="algorithm", default="probability", metavar="ALGORITHM",
                      type="choice", help=help["-l"], choices=["probability", "express", "regression",
                      "regression-probability", "regression-express"])
    parser.add_option("-c", "--chromosomes", dest="chroms", default="", metavar="CHROMS", type="string",
                      help=help["-c"], action="store")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    group = optparse.OptionGroup(parser, "Fend Feature-Specific Options", "Exactly one of these options needs to be specified.")
    group.add_option("-F", "--fend", dest="fend", default=None, metavar="FEND", type="string", action="store",
                      help=help["-F"])
    group.add_option("-B", "--bed", dest="bed", default=None, metavar="BED", type="string", action="store",
                      help=help["-B"])
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Data-Specific Options", "Exactly one of these option types needs to be specified. -S/--bam and -R/--raw can be used multiple times to specify multiple input data files.")
    group.add_option("-S", "--bam", dest="bam", default=[], metavar="BAM1 BAM2", type="string", action="append",
                      help=help["-S"], nargs=2)
    group.add_option("-R", "--raw", dest="raw", default=[], metavar="RAW", type="string", action="append",
                      help=help["-R"])
    group.add_option("-M", "--mat", dest="mat", default=None, metavar="MAT", type="string", action="store",
                     help=help["-M"])
    parser.add_option_group(group)
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
    group.add_option("-w", "--express-reads", dest="expreads", default="cis", metavar="E-READS", type="choice",
                      help=help["-w"], choices=["cis", "trans", "all"])
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Regression Algorithm-Specific Options")
    group.add_option("-r", "--regression-iterations", dest="regiter", default=1000, metavar="REGITER", type="int",
                      help=help["-r"], action="store")
    group.add_option("-t", "--learning-threshold", dest="threshold", default=1.0, metavar="THRESHOLD", type="float",
                      help=help["-t"], action="store")
    group.add_option("-y", "--regression-reads", dest="regreads", default="cis", metavar="R-READS", type="choice",
                      help=help["-w"], choices=["cis", "trans", "all"])
    group.add_option("-v", "--model", dest="model", default="gc,len,distance", metavar="MODEL", type="string",
                      help=help["-v"], action="store")
    group.add_option("-n", "--model-bins", dest="modelbins", default="20,20,15", metavar="NUMBINS", type="string",
                      help=help["-n"], action="store")
    parser.add_option_group(group)
    return parser

def check_options(parser, options, args):
    if len(args) < 3:
        if rank == 0:
            parser.error('incorrect number of arguments')
        else:
            sys.exit(1)
    if ((options.bed is None and options.fend is None) or
        (not options.bed is None and not options.fend is None)):
        if rank == 0:
            parser.error('requires exactly one fend file type (-F/--fend or -B/--bed)')
        else:
            sys.exit(1)
    if ((len(options.bam) > 0 and (len(options.raw) > 0 or not options.mat is None)) or
        (len(options.raw) > 0 and (len(options.bam) > 0 or not options.mat is None)) or
        (not options.mat is None and (len(options.raw) > 0 or len(options.bam) > 0)) or
        (len(options.bam) == 0 and len(options.raw) == 0 and options.mat is None)):
        if rank == 0:
            parser.error('requires exactly one data argument type (-S/--bam, -R/--raw, or -M/--mat)')
        else:
            sys.exit(1)
    options.model = options.model.split(',')
    options.modelbins = options.modelbins.split(',')
    for i in range(len(options.modelbins)):
        options.modelbins[i] = int(options.modelbins[i])
    if (len(options.model) != len(options.modelbins) and
        options.algorithm in ['regression', 'regression-express', 'regression-probability']):
        if rank == 0:
            parser.error("'model' and 'model-bins' options must be the same length")
        else:
            sys.exit(1)
    options.chroms = options.chroms.split(',')
    return None

def display_options(options, args):
    settings =      "Running with the following settings:\n"
    settings +=     "Fend file destination:          %s\n" % args[0]
    settings +=     "Data file destination:          %s\n" % args[1]
    settings +=     "Project file destination:       %s\n" % args[2]
    if not options.fend is None:
        settings += "Fend data read from:            %s\n" % options.fend
    else:
        settings += "Fend data read from:            %s\n" % options.bed
    if len(options.bam) > 0:
        datafiles = "%s & %s" % (options.bam[0][0], options.bam[0][1])
        if len(options.bam) == 2:
            datafiles += " and %s & %s" % (options.bam[1][0], options.bam[1][1])
        elif len(options.bam) > 2:
            for i in range(1, len(options.bam) - 1):
                datafiles += ", %s & %s" % (options.bam[i][0], options.bam[i][1])
            datafiles += ", and %s & %s" % (options.bam[-1][0], options.bam[-1][1])
    elif len(options.raw) > 0:
        datafiles = "%s" % (options.raw[0])
        if len(options.raw) == 2:
            datafiles += " and %s" % (options.raw[1])
        elif len(options.raw) > 2:
            for i in range(1, len(options.raw) - 1):
                datafiles += ", %s" % (options.raw[i])
            datafiles += ", and %s" % (options.raw[-1])
    else:
        datafiles = options.mat
    settings +=     "Interaction data read from:     %s\n" % datafiles
    settings +=     "Maximum insert size:            %i\n" % options.insert
    settings +=     "Minimum interactions per fend:  %i\n" % options.minint
    settings +=     "Minimum interaction distance:   %i\n" % options.mindist
    if options.maxdist == 0:
        maxdist = "no limit"
    else:
        maxdist = str(options.maxdist)
    settings +=     "Maximum interaction distance:   %s\n" % maxdist
    settings +=     "Learning algorithm:             %s\n" % options.algorithm
    if options.chroms == "":
        chroms = "all"
    else:
        temp = options.chroms
        chroms = "%s" % (temp[0])
        if len(temp) == 2:
            chroms += " and %s" % (temp[1])
        elif len(temp) > 2:
            for i in range(1, len(temp) - 1):
                chroms += ", %s" % (temp[i])
            chroms += ", and %s" % (temp[-1])
    settings +=     "Chromosomes to learn:           %s\n" % chroms
    if options.algorithm in ['regression', 'regression-probability', 'regression-express']:
        settings += "Regression Algorithm Settings\n"
        settings += "Regression model parameters:    %s\n" % ','.join(options.model)
        temp = '%i' % options.modelbins[0]
        for i in range(1, len(options.modelbins)):
            temp += ',%i' % options.modelbins[i]
        settings += "Regression model bin sizes:     %s\n" % temp
        settings += "Reads to include in learning:   %s\n" % options.regreads
        settings += "Number of regression iterations:%i\n" % options.regiter
        settings += "Regression learning cutoff:     %i\n" % options.threshold
    if options.algorithm in ['express', 'regression-express']:
        settings += "Express Algorithm Settings\n"
        settings += "Remove distance-dependence:     %s\n" % str(options.nodist)
        settings += "Reads to include in learning:   %s\n" % options.expreads
        settings += "Number of learning iterations:  %i\n" % options.expiter
    if options.algorithm in ['probability', 'regression-probability']:
        settings += "Probability Algorithm Settings\n"
        settings += "Precalculate corrections:       %s\n" % str(options.precalc)
        settings += "Learning rate:                  %f\n" % options.rate
        settings += "Correction change cutoff:       %f\n" % options.change
        if options.display == 0:
            display = "off"
        else:
            display = str(options.display)
        settings += "Iterations between display:     %s\n" % display
        settings += "Number of burnin iterations:    %i\n" % options.burnin
        settings += "Number of annealing iterations: %i\n" % options.anneal
    print settings
    return None


if 'mpi4py' in sys.modules.keys():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()
else:
    comm = None
    rank = 0
    num_procs = 1

if __name__ == "__main__":
    main()
