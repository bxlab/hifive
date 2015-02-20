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
    frags.load_fragments(args[0])
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
    fivec.filter_fragments(mininteractions=options.minint)
    fivec.find_distance_parameters()
    if options.express:
        fivec.find_express_fragment_corrections(iterations=options.iter, mindistance=options.mindist,
                                                maxdistance=options.maxdist, remove_distance=options.nodist)
    else:
        fivec.find_fragment_corrections(display=options.display, mindistance=options.mindist,
                                        maxdistance=options.maxdist, learningrate=options.rate,
                                        burnin_iterations=options.burnin, annealing_iterations=options.anneal,
                                        precalculate=options.precalc)
    fivec.save()
    return None


def setup_parser():
    usage = "usage: %prog [options] <frag_file> <out_frag> <out_data> <out_project>\n\nArguments:"
    usage += "\n<frag_file>  a file containing probed restriction enzyme fragment boundaries in BED format"
    usage += "\n<out_frag>   a file to write the HiFive 5C fragment file to"
    usage += "\n<out_data>   a file to write the HiFive 5C data file to"
    usage += "\n<out_project>   a file to write the HiFive 5C project file to"
    help = {
        "BAM":"a pair of comma-separated bam read end file lists (first and second end) for a single sequencing run. For multiple runs, this option can be passed multiple times",
        "COUNT":"a tab-separated text file containing pairs of fragments and their associated number of observed reads (fragment1 fragment2 count), one per line. For multiple files, this option can be passed multiple times",
        "MININT":"minimum number of interactions needed for valid fragments [default: %default]",
        "MINDIST":"smallest interaction distance to be included for filtering and learning fragment corrections [default: %default]",
        "MAXDIST":"largest interaction distance to be included for filtering and, if not using Express, learning fragment corrections (zero indicates no maximum) [default: %default]",
        "EXPRESS":"indicates that the express learning algorithm should be used",
        "BURNIN":"number of iterations to run burn-in phase for under the standard learning algorithm with constant learning rate [default: %default]",
        "ANNEAL":"number of iterations to run annealing-in phase for under the standard learning algorithm with decreasing learning rate [default: %default]",
        "ITER":"number of iterations to run learning for under Express learning algorithm [default %default]",
        "PRECALC":"precalculate correction values from fend means",
        "RATE":"learning rate, defined as percentage of gradient to use for updating parameter values [default: %default]",
        "REMOVE":"remove the distant-dependent portion of the signal prior to learning corrections [default: %default]",
        "DISPLAY":"number of iterations to wait before explicitly calculating cost and updating display (zero indicates off) [default: %default]",
        "QUIET":"silence output messages [default: %default]",
    }
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-B", "--bam", dest="bam", default=[], metavar="BAM1 BAM2", type="string", action="append",
                      help=help["BAM"], nargs=2)
    parser.add_option("-C", "--count", dest="counts", default=[], metavar="COUNT", type="string", action="append",
                      help=help["COUNT"])
    parser.add_option("-t", "--min-interactions", dest="minint", default=20, metavar="MININT", type="int",
                      action="store", help=help["MININT"])
    parser.add_option("-m", "--min-distance", dest="mindist", default=0, metavar="MINDIST", type="int", action="store",
                      help=help["MINDIST"])
    parser.add_option("-x", "--max-distance", dest="maxdist", default=0, metavar="MAXDIST", type="int", action="store",
                      help=help["MAXDIST"])
    parser.add_option("-l", "--express", dest="express", action="store_true", default=False, help=help["EXPRESS"])
    parser.add_option("-p", "--precalculate", dest="precalc", action="store_true", default=False, help=help["PRECALC"])
    parser.add_option("-b", "--burnin", dest="burnin", default=1000, metavar="BURNIN", type="int", action="store",
                      help=help["BURNIN"])
    parser.add_option("-a", "--anneal", dest="anneal", default=1000, metavar="ANNEAL", type="int", action="store",
                      help=help["ANNEAL"])
    parser.add_option("-e", "--iterations", dest="iter", default=1000, metavar="ITER", type="int", action="store",
                      help=help["ITER"])
    parser.add_option("-r", "--learning-rate", dest="rate", default=0.1, metavar="RATE", type="float", action="store",
                      help=help["RATE"])
    parser.add_option("-d", "--display", dest="display", default=100, metavar="DISPLAY", type="int", action="store",
                      help=help["DISPLAY"])
    parser.add_option("-v", "--remove-distance", dest="nodist", action="store_true", default=False,
                      help=help["REMOVE"])
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False, help=help["QUIET"])
    return parser

def check_options(parser, options, args):
    if len(args) < 4:
        parser.error('incorrect number of arguments')
    if ((len(options.bam) > 0 and len(options.counts) > 0) or
        (len(options.bam) == 0 and len(options.counts) == 0)):
        parser.error('requires exactly one data argument type (-B/--bam or -C/--count)')
    return None

def display_options(options, args):
    settings =      "Running with the following settings:\n"
    settings +=     "Output file prefix:              %s\n" % args[1]
    settings +=     "Fragment data read from:         %s\n" % args[0]
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
    if options.express:
        settings += "Learning algorithm:              approximate (Express)\n"
        settings += "Remove distance-dependence:      %s\n" % str(options.nodist)
        settings += "Number of learning iterations:   %i\n" % options.iter
    else:
        settings += "Learning algorithm:              probabilistic (Standard)\n"
        settings += "Precalculate corrections:        %s\n" % str(options.precalc)
        settings += "Learning rate:                   %f\n" % options.rate
        settings += "Number of burnin iterations:     %i\n" % options.burnin
        settings += "Number of annealing iterations:  %i\n" % options.anneal
        if options.display == 0:
            display = "off"
        else:
            display = str(options.display)
        settings += "Iterations between display:      %s\n" % display
    print settings
    return None

if __name__ == "__main__":
    main()
