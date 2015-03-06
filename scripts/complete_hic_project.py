#!/usr/bin/env python

import sys
import optparse

import hifive


def main():
    parser = setup_parser()
    options, args = parser.parse_args()
    check_options(parser, options, args)
    if not options.silent:
        display_options(options, args)
    fend_fname = args[0]
    data_fname = args[1]
    hic_fname = args[2]
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
    hic = hifive.HiC(hic_fname, 'w', silent=options.silent)
    hic.load_data(data_fname)
    hic.filter_fends(mininteractions=options.minint, mindistance=options.mindist, maxdistance=options.maxdist)
    hic.find_distance_parameters(minsize=options.minbin, numbins=options.numbin)
    if options.algorithm in ['regression', 'regression-probability', 'regression-express']:
        hic.find_regression_fend_corrections(max_iterations=options.iter, chroms=options.chroms.split(','),
                                             mindistance=options.mindist, maxdistance=options.maxdist,
                                             model=options.model, num_bins=options.modelbins)
    if options.algorithm == 'express':
        hic.find_express_fend_corrections(iterations=options.iter, mindistance=options.mindist,
                                          maxdistance=options.maxdist, mininteractions=options.minint,
                                          usereads=options.reads, remove_distance=options.nodist,
                                          chroms=options.chroms.split(','))
    elif options.algorithm == 'regression-express':
        hic.find_express_fend_corrections(iterations=options.iter, mindistance=options.mindist,
                                          maxdistance=options.maxdist, mininteractions=options.minint,
                                          usereads=options.reads, remove_distance=options.nodist,
                                          chroms=options.chroms.split(','), precorrect=True)
    elif options.algorithm == 'probability':
        hic.find_fend_corrections(display=options.display, mindistance=options.mindist, maxdistance=options.maxdist,
                                  chroms=options.chroms.split(','), learningrate=options.rate,
                                  burnin_iterations=options.burnin, annealing_iterations=options.anneal,
                                  minchange=options.change, precalculate=options.precalc)
    elif options.algorithm == 'regression-probability':
        hic.find_fend_corrections(display=options.display, mindistance=options.mindist, maxdistance=options.maxdist,
                                  chroms=options.chroms.split(','), learningrate=options.rate,
                                  burnin_iterations=options.burnin, annealing_iterations=options.anneal,
                                  minchange=options.change, precalculate=options.precalc, precorrect=True)
    hic.save()
    return None


def setup_parser():
    usage = "usage: [mpirun -np N_PROCS] %prog [options] <out_fend> <out_data> <out_project>\n\nArguments:"
    usage += "\n<out_fend>   a file to write the HiFive HiC fend file to"
    usage += "\n<out_data>   a file to write the HiFive HiC data file to"
    usage += "\n<out_project>   a file to write the HiFive HiC project file to"
    help = {
        "FEND":"fend file in HiCPipe-compatible MAT format",
        "BED":"restriction enzyme cut sites or RE fragment boundaries in BED format",
        "BAM":"a pair of comma-separated bam read end file lists (first and second end) for a single sequencing run. For multiple runs, this option can be passed multiple times",
        "MAT":"a HiCPipe-style MAT file containing fend pair counts",
        "RAW":"a tab-separated text file containing pairs of read ends (chr1 pos1 strand1 chr2 pos2 strand2), one per line. For multiple files, this option can be passed multiple times",
        "INSERT":"maximum allowable distance sum between both fend ends and cutsites [default: %default]",
        "MININT":"minimum number of interactions needed for valid fends [default: %default]",
        "MINDIST":"smallest interaction distance to be included for filtering and learning fend corrections [default: %default]",
        "MAXDIST":"largest interaction distance to be included for filtering and, if not using Express, learning fend corrections (zero indicates no maximum) [default: %default]",
        "MINBIN":"smallest interaction distance bin size for distance function [default: %default]",
        "NUMBIN":"number of bins to partion interaction distance range into for distance function [default: %default]",
        "ALGORITHM":"select the learning algorithm that should be used (probability, express, regression, regression-probability, or regression-express) [default: %default]",
        "BURNIN":"number of iterations to run burn-in phase for under the standard learning algorithm with constant learning rate [default: %default]",
        "ANNEAL":"number of iterations to run annealing-in phase for under the standard learning algorithm with decreasing learning rate [default: %default]",
        "ITER":"number of iterations to run learning for under Express or Regression learning algorithm [default %default]",
        "CHANGE":"minimum mean change in the fend correction parameter values needed to keep running past 'burnin_iterations' number of iterations during burn-in phase [default: %default]",
        "PRECALC":"precalculate correction values from fend means",
        "RATE":"learning rate, defined as percentage of gradient to use for updating parameter values [default: %default]",
        "READS":"which set of reads, 'cis', 'trans', or 'both', to use for learning [default: %default]",
        "REMOVE":"remove the distant-dependent portion of the signal prior to learning corrections [default: %default]",
        "DISPLAY":"number of iterations to wait before explicitly calculating cost and updating display (zero indicates off) [default: %default]",
        "MODEL":"comma-separated list of parameters to include in regression model (gc, len, mappability, and distance) [default: %default]",
        "MODELBINS":"comma-separated list of number of bins, one for each model parameter [default: %default]",
        "CHROMS":"comma-separated list of chromosome names to learn corrections for [default: all chromosomes]",
        "QUIET":"silence output messages [default: %default]",
    }
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-F", "--fend", dest="fend", default=None, metavar="FEND", type="string", action="store",
                      help=help["FEND"])
    parser.add_option("-B", "--bed", dest="bed", default=None, metavar="BED", type="string", action="store",
                      help=help["BED"])
    parser.add_option("-S", "--bam", dest="bam", default=[], metavar="BAM1 BAM2", type="string", action="append",
                      help=help["BAM"], nargs=2)
    parser.add_option("-R", "--raw", dest="raw", default=[], metavar="RAW", type="string", action="append",
                      help=help["RAW"])
    parser.add_option("-M", "--mat", dest="mat", default=None, metavar="MAT", type="string", action="store",
                      help=help["MAT"])
    parser.add_option("-i", "--insert", dest="insert", default=500, metavar="INSERT", type="int", action="store",
                      help=help["INSERT"])
    parser.add_option("-t", "--min-interactions", dest="minint", default=20, metavar="MININT", type="int",
                      action="store", help=help["MININT"])
    parser.add_option("-m", "--min-distance", dest="mindist", default=0, metavar="MINDIST", type="int", action="store",
                      help=help["MINDIST"])
    parser.add_option("-x", "--max-distance", dest="maxdist", default=0, metavar="MAXDIST", type="int", action="store",
                      help=help["MAXDIST"])
    parser.add_option("-s", "--min-binsize", dest="minbin", default=1000, metavar="MINBIN", type="int",
                      action="store", help=help["MINBIN"])
    parser.add_option("-n", "--num-bins", dest="numbin", default=100, metavar="NUMBIN", type="int", action="store",
                      help=help["NUMBIN"])
    parser.add_option("-a", "--algorithm", dest="algorithm", action="choice", default='probability',
                      choices=["probability", "express", "regression", "regression-probability", "regression-express"],
                      help=help["ALGORITHM"], metavar="ALGORITHM")
    parser.add_option("-p", "--precalculate", dest="precalc", action="store_true", default=False, help=help["PRECALC"])
    parser.add_option("-b", "--burnin", dest="burnin", default=1000, metavar="BURNIN", type="int", action="store",
                      help=help["BURNIN"])
    parser.add_option("-a", "--anneal", dest="anneal", default=1000, metavar="ANNEAL", type="int", action="store",
                      help=help["ANNEAL"])
    parser.add_option("-e", "--iterations", dest="iter", default=1000, metavar="ITER", type="int", action="store",
                      help=help["ITER"])
    parser.add_option("-g", "--min-change", dest="change", default=0.0001, metavar="CHANGE", type="float",
                      action="store", help=help["CHANGE"])
    parser.add_option("-r", "--learning-rate", dest="rate", default=0.1, metavar="RATE", type="float", action="store",
                      help=help["RATE"])
    parser.add_option("-d", "--display", dest="display", default=100, metavar="DISPLAY", type="int", action="store",
                      help=help["DISPLAY"])
    parser.add_option("-v", "--remove-distance", dest="nodist", action="store_true", default=False,
                      help=help["REMOVE"])
    parser.add_option("-f", "--reads", dest="reads", default="cis", metavar="READS", type="choice",
                      action="store", help=help["READS"], choices=["cis", "trans", "both"])
    parser.add_option("-h", "--model", dest="model", default="gc,len,distance", metavar="MODEL", type="store",
                      help=help["MODEL"])
    parser.add_option("-j", "--modelbins", dest="modelbins", default="20,20,15", metavar="MODELBINS", type="store",
                      help=help["MODELBINS"])
    parser.add_option("-c", "--chromosomes", dest="chroms", default="", metavar="CHROMS", type="string",
                      action="store", help=help["CHROMS"])
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False, help=help["QUIET"])
    return parser

def check_options(parser, options, args):
    if len(args) < 3:
        parser.error('incorrect number of arguments')
    if ((options.bed is None and options.fend is None) or
        (not options.bed is None and not options.fend is None)):
        parser.error('requires exactly one fend file type (-F/--fend or -B/--bed)')
    if ((len(options.bam) > 0 and (len(options.raw) > 0 or not options.mat is None)) or
        (len(options.raw) > 0 and (len(options.bam) > 0 or not options.mat is None)) or
        (not options.mat is None and (len(options.raw) > 0 or len(options.bam) > 0)) or
        (len(options.bam) == 0 and len(options.raw) == 0 and options.mat is None)):
        parser.error('requires exactly one data argument type (-S/--bam, -R/--raw, or -M/--mat)')
    options.model = options.model.split(',')
    options.modelbins = options.modelbins.split(',')
    for i in range(len(options.modelbins)):
        options.modelbins[i] = int(options.modelbins[i])
    if len(options.model) != len(options.modelbins):
        parser.error("'model' and 'modelbins' options must be the same length")
    return None

def display_options(options, args):
    settings =      "Running with the following settings:\n"
    settings +=     "Output file prefix:             %s\n" % args[0]
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
    settings +=     "Minimum interactions per fend:  %i\n" % options.minint
    settings +=     "Minimum interaction distance:   %i\n" % options.mindist
    if options.maxdist == 0:
        maxdist = "no limit"
    else:
        maxdist = str(options.maxdist)
    settings +=     "Maximum interaction distance:   %s\n" % maxdist
    settings +=     "Learning algorithm:             %s\n" % options.algorithm
    if options.algorithm in ['express', 'regression-express']:
        settings += "Remove distance-dependence:     %s\n" % str(options.nodist)
        settings += "Reads to include in learning:   %s\n" % options.reads
        settings += "Number of learning iterations:  %i\n" % options.iter
    elif options.algorithm == 'regression':
        settings += "Number of learning iterations:  %i\n" % options.iter
    else:
        settings += "Learning algorithm:             probabilistic (Standard)\n"
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
    if options.algorithm == 'regression-probability':
        settings += "Number of regression iterations:%i\n" % options.iter
    if options.algorithm in ['regression', 'regression-probability', 'regression-express']:
        settings += "Regression model parameters:   %s\n" % ','.join(options.model)
        temp = '%i' % options.modelbins[0]
        for i in range(1, len(options.modelbins)):
            temp += ',%i' % options.modelbins[i]
        settings += "Regression model bin sizes:    %s\n" % temp
    if options.chroms == "":
        chroms = "all"
    else:
        temp = options.chroms.split(',')
        chroms = "%s" % (temp[0])
        if len(temp) == 2:
            chroms += " and %s" % (temp[1])
        elif len(temp) > 2:
            for i in range(1, len(temp) - 1):
                chroms += ", %s" % (temp[i])
            chroms += ", and %s" % (temp[-1])
    settings +=     "Chromosomes to learn:           %s" % chroms
    print settings
    return None

if __name__ == "__main__":
    main()
