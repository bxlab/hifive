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
    frag_fname = args[1]
    data_fname = args[2]
    project_fname = args[3]
    frags = hifive.Fragment(frag_fname, mode='w', silent=options.silent)
    frags.load_fragments(args[0], fastafile=options.fasta)
    frags.fragments.close()
    del frags
    data = hifive.FiveCData(data_fname, 'w', silent=options.silent)
    if len(options.bam) > 0:
        data.load_data_from_bam(frag_fname, options.bam)
    else:
        data.load_data_from_counts(frag_fname, options.counts)
    data.save()
    del data
    fivec = hifive.FiveC(project_fname, 'w', silent=options.silent)
    fivec.load_data(data_fname)
    fivec.filter_fragments(mininteractions=options.minint, mindistance=options.mindist, maxdistance=options.maxdist)
    fivec.find_distance_parameters()
    precorrect = False
    if options.algorithm in ['regression', 'regression-express', 'regression-probability']:
        fivec.find_regression_fragment_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                                   regions=options.regions, num_bins=options.modelbins,
                                                   model=options.model, usereads=options.regreads,
                                                   learning_threshold=options.threshold,
                                                   max_iterations=options.regiter)
        precorrect = True
    if options.algorithm in ['probability', 'regression-probability']:
        fivec.find_probability_fragment_corrections(mindistance=options.mindist, maxdistance=options.maxdist,
                                                    regions=options.regions, burnin_iterations=options.burnin,
                                                    annealing_iterations=options.anneal, learningrate=options.rate,
                                                    precalculate=options.precalc, precorrect=precorrect)
    elif options.algorithm in ['express', 'regression-express']:
        fivec.find_express_fragment_corrections(iterations=options.expiter, mindistance=options.mindist,
                                                maxdistance=options.maxdist, remove_distance=options.nodist,
                                                usereads=options.expreads, regions=options.regions,
                                                precorrect=precorrect)
    fivec.save(options.outfname)
    return None


def setup_parser():
    usage = "usage: %prog [options] <frag_file> <out_frag> <out_data> <out_project>\n\nArguments:"
    usage += "\n<frag_file>  a file containing probed restriction enzyme fragment boundaries in BED format"
    usage += "\n<out_frag>   a file to write the HiFive 5C fragment file to"
    usage += "\n<out_data>   a file to write the HiFive 5C data file to"
    usage += "\n<out_project>   a file to write the HiFive 5C project file to"
    help = {
        "-f":"FASTA file containing sequence of 5C primers (for finding GC content for regression modeling)"
        "-n":"minimum number of interactions needed for valid fragment [default: %default]",
        "-l":"learning algorithm to use (probability, express, regression, regression-probability, or regression-express) [default: %default]",
        "-m":"minimum interaction distance to include in learning [default: %default]",
        "-x":"maximum interaction distance to include in learning [default: %default]",
        "-r":"comma-separated list of region indices to learn corrections for [default: all regions]",
        "-o":"output file to save analyzed copy of HiFive 5C project to. If not specified, original file will be overwritten",
        "-q":"silence output messages [default: %default]",
        "-B":"a pair of comma-separated bam read end file lists (first and second end) for a single sequencing run. For multiple runs, this option can be passed multiple times",
        "-C":"a tab-separated text file containing pairs of fragments and their associated number of observed reads (fragment1 fragment2 count), one per line. For multiple files, this option can be passed multiple times",
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
    parser.add_option("-f", "--fasta-file", dest="fasta", default=None, metavar="FASTA", type="string",
                      help=help["-f"], action="store")
    parser.add_option("-n", "--min-interactions", dest="minint", default=20, metavar="MININT", type="int",
                      action="store", help=help["-n"])
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
    group = optparse.OptionGroup(parser, "Data-Specific Options", "Exactly one of these option types needs to be specified. Both options can be used multiple times to specify multiple input data files.")
    group.add_option("-B", "--bam", dest="bam", default=[], metavar="BAM1 BAM2", type="string", action="append",
                      help=help["-B"], nargs=2)
    group.add_option("-C", "--count", dest="counts", default=[], metavar="COUNT", type="string", action="append",
                      help=help["-C"])
    parser.add_option_group(group)
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
    return parser

def check_options(parser, options, args):
    if len(args) < 4:
        parser.error('incorrect number of arguments')
    return None
    if len(args) < 4:
        parser.error('incorrect number of arguments')
    if ((len(options.bam) > 0 and len(options.counts) > 0) or
        (len(options.bam) == 0 and len(options.counts) == 0)):
        parser.error('requires exactly one data argument type (-B/--bam or -C/--count)')
    options.model = options.model.split(',')
    options.modelbins = options.modelbins.split(',')
    for i in range(len(options.modelbins)):
        options.modelbins[i] = int(options.modelbins[i])
    if (len(options.model) != len(options.modelbins) and
        options.algorithm in ['regression', 'regression-express', 'regression-probability']):
        parser.error("'model' and 'model-bins' options must be the same length")
    options.regions = options.regions.split(',')
    if len(options.regions) == 1 and options.regions[0] == '':
        options.regions = []
    for i in range(len(options.regions)):
        options.regions[i] = int(options.regions[i])
    return None

def display_options(options, args):
    settings =      "Running with the following settings:\n"
    settings +=     "Fragment data read from:         %s\n" % args[0]
    settings +=     "Frag file destination:           %s\n" % args[1]
    settings +=     "Data file destination:           %s\n" % args[2]
    settings +=     "Project file destination:        %s\n" % args[3]
    if not options.fasta is None:
        settings += "Primer Fasta file:               %s\n" % options.fasta
    if len(options.bam) > 0:
        datafiles = "%s & %s" % (options.bam[0][0], options.bam[0][1])
        if len(options.bam) == 2:
            datafiles += " and %s & %s" % (options.bam[1][0], options.bam[1][1])
        elif len(options.bam) > 2:
            for i in range(1, len(options.bam) - 1):
                datafiles += ", %s & %s" % (options.bam[i][0], options.bam[i][1])
            datafiles += ", and %s & %s" % (options.bam[-1][0], options.bam[-1][1])
    else:
        datafiles = "%s" % (options.counts[0])
        if len(options.counts) == 2:
            datafiles += " and %s" % (options.counts[1])
        elif len(options.counts) > 2:
            for i in range(1, len(options.counts) - 1):
                datafiles += ", %s" % (options.counts[i])
            datafiles += ", and %s" % (options.counts[-1])
    settings +=     "Interaction data read from:      %s\n" % datafiles
    settings +=     "Minimum interactions per frag:   %i\n" % options.minint
    settings +=     "Minimum interaction distance:    %i\n" % options.mindist
    if options.maxdist == 0:
        maxdist = "no limit"
    else:
        maxdist = str(options.maxdist)
    settings +=     "Maximum interaction distance:    %s\n" % maxdist
    settings +=     "Learning algorithm:              %s\n" % options.algorithm
    if options.regions == "":
        regions = "all"
    else:
        temp = []
        for i in range(len(options.regions)):
            temp.append(str(options.regions[i]))
        regions = "%s" % (temp[0])
        if len(temp) == 2:
            regions += " and %s" % (temp[1])
        elif len(temp) > 2:
            for i in range(1, len(temp) - 1):
                regions += ", %s" % (temp[i])
            regions += ", and %s" % (temp[-1])
    settings +=     "Regions to learn:                %s\n" % regions
    if options.algorithm in ['regression', 'regression-probability', 'regression-express']:
        settings += "Regression Algorithm Settings\n"
        settings += "Regression model parameters:     %s\n" % ','.join(options.model)
        temp = '%i' % options.modelbins[0]
        for i in range(1, len(options.modelbins)):
            temp += ',%i' % options.modelbins[i]
        settings += "Regression model bin sizes:      %s\n" % temp
        settings += "Reads to include in learning:    %s\n" % options.regreads
        settings += "Number of regression iterations: %i\n" % options.regiter
        settings += "Regression learning cutoff:      %i\n" % options.threshold
    if options.algorithm in ['express', 'regression-express']:
        settings += "Express Algorithm Settings\n"
        settings += "Remove distance-dependence:      %s\n" % str(options.nodist)
        settings += "Reads to include in learning:    %s\n" % options.expreads
        settings += "Number of learning iterations:   %i\n" % options.expiter
    if options.algorithm in ['probability', 'regression-probability']:
        settings += "Probability Algorithm Settings\n"
        settings += "Precalculate corrections:        %s\n" % str(options.precalc)
        settings += "Learning rate:                   %f\n" % options.rate
        settings += "Number of burnin iterations:     %i\n" % options.burnin
        settings += "Number of annealing iterations:  %i\n" % options.anneal
    print settings
    return None

if __name__ == "__main__":
    main()
