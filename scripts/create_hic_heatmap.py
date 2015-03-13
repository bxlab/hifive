#!/usr/bin/env python

import sys
import optparse
import subprocess

try:
    from mpi4py import MPI
except:
    pass
try:
    from pyx import *
except:
    pass
import h5py

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
    usage = "usage: [mpirun -np N_PROCS] %prog [options] <project_file> <out_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive HiC project file"
    usage += "\n<out_file>      destination for HiC heatmaps file"
    help = {
        "-b":"size of bins, in base pairs, to group data into (a size of zero indicates unbinned) [default: %default]",
        "-t":"calculate and include trans interactions in heatmaps [default: %default]",
        "-d":"which corrections (if any) to apply to counts (raw, distance, fend, or enrichment) [default: %default]",
        "-c":"comma-separated list of chromosomes to include in heatmaps [default: all chromosomes]",
        "-i":"save the data as an image to this file",
        "-p":"format image in PDF file [default %default]",
        "-l":"add color scale to plot (pdf format only) [default %default]",
        "-n":"add chromosome labels to plot (pdf format only) [default %default]",
        "-k":"additional keyword arguments to pass function",
        "-q":"silence output messages [default: %default]",
    }
    if rank == 0:
        parser = optparse.OptionParser(usage=usage)
    else:
        parser = optparse.OptionParser(usage=optparse.SUPPRESS_USAGE, add_help_option=False)
        for key in help:
            help[key] = optparse.SUPPRESS_HELP
        parser.add_option("-h", "--help", help=optparse.SUPPRESS_HELP, dest="help", action="store_true")
    parser.add_option("-b", "--bin-size", dest="binsize", default=10000, metavar="BIN", type="int",
                      help=help["-b"], action="store")
    parser.add_option("-t", "--trans", dest="trans", default=False,
                      help=help["-t"], action="store_true")
    parser.add_option("-d", "--datatype", dest="datatype", default='enrichment', action="store",
                      help=help["-d"], type="choice", choices=['raw', 'distance', 'fend', 'enrichment'])
    parser.add_option("-c", "--chromosomes", dest="chroms", default="", metavar="CHROMS", type="string",
                      help=help["-c"], action="store")
    parser.add_option("-i", "--image-file", dest="imagefile", default=None, metavar="IMAGE", type="string",
                      help=help["-i"], action="store")
    parser.add_option("-p", "--pdf", dest="pdf", default=False,
                      help=help["-p"], action="store_true")
    parser.add_option("-l", "--add-legend", dest="legend", default=False,
                      help=help["-l"], action="store_true")
    parser.add_option("-n", "--names", dest="names", default=False,
                      help=help["-n"], action="store_true")
    parser.add_option("-k", "--keyword", dest="keywords", default=[], metavar="KEYWORD", type="string",
                      help=help["-k"], action="append")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help=help["-q"])
    options, args = parser.parse_args()
    if len(args) < 2:
        if rank == 0:
            parser.error('incorrect number of arguments')
        else:
            sys.exit(1)
    if not options.imagefile is None and options.pdf and "pyx" not in sys.modules.keys():
        if rank == 0:
            parser.error("-p/--pdf requires the package 'pyx'")
        else:
            sys.exit(1)
    options.regions = options.chroms.split(',')
    if len(options.chroms) == 1 and options.chroms[0] == '':
        options.chroms = []
    fivec = hifive.FiveC(args[0], 'r', silent=options.silent)
    hic = hifive.HiC(args[0], 'r', silent=options.silent)
    hic.write_heatmap(args[1], binsize=options.binsize, includetrans=options.trans,
                      datatype=options.datatype, chroms=options.chroms)
    if rank > 0:
        sys.exit(0)
    if not options.imagefile is None:
        kwargs = {}
        for arg in options.keywords:
            temp = arg.split("=")
            if temp[1] in ["True", "TRUE", "true"]:
                temp[1] = True
            elif temp[1] in ["False", "FALSE", "false"]:
                temp[1] = False
            elif temp[1][0] == "(":
                temp[1] = temp[1].strip('()').split(',')
                for i in range(len(temp[1])):
                    temp[1][i] = int(temp[1][i])
                temp[1] = tuple(temp[1])
            elif temp[1][0] == "[":
                temp[1] = temp[1].strip('[]').split(',')
                for i in range(len(temp[1])):
                    temp[1][i] = int(temp[1][i])
            else:
                try:
                    temp[1] = int(temp[1])
                except:
                    try:
                        temp[1] = float(temp[1])
                    except:
                        temp[1] = temp[1].replace('__pd__','')
            kwargs[temp[0]] = temp[1]
        if 'symmetricscaling' not in kwargs:
            if options.datatype == 'enrichment':
                kwargs['symmetricscaling'] = True
            else:
                kwargs['symmetricscaling'] = False
        img, minscore, maxscore = hifive.plotting.plot_hic_heatmap(args[1], returnscale=True,
                                                                   silent=options.silent, **kwargs)
        if not options.pdf:
            img.save(options.imagefile, format='png')
        else:
            unit.set(defaultunit="cm")
            text.set(mode="latex")
            text.preamble(r"\usepackage{times}")
            text.preamble(r"\usepackage{sansmath}")
            text.preamble(r"\sansmath")
            text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
            hm = h5py.File(args[1], 'r')
            chroms = hm['chromosomes'][...]
            sizes = [0]
            minsize = 999999999
            for chrom in chroms:
                sizes.append(hm['%s.positions' % chrom].shape[0])
                minsize = min(minsize, sizes[-1])
            if len(sizes) > 2:
                sizes[1] += 0.5
                sizes[-1] += 0.5
                if len(sizes) > 3:
                    for i in range(2, len(sizes) - 1):
                        sizes[i] += 1.0
            for i in range(1, len(sizes)):
                sizes[i] += sizes[i - 1]
            height = width = max(5.0, sizes[-1] * 0.5 / minsize)
            for i in range(len(sizes)):
                sizes[i] = sizes[i] / sizes[-1] * height
            c = canvas.canvas()
            c.insert(bitmap.bitmap(0, 0, img, width=width))
            if options.legend:
                if 'min_color' in kwargs:
                    min_color = kwargs['min_color']
                else:
                    min_color = "0000ff"
                if 'mid_color' in kwargs:
                    mid_color = kwargs['mid_color']
                else:
                    mid_color = "ffffff"
                if 'max_color' in kwargs:
                    max_color = kwargs['max_color']
                else:
                    max_color = "ff0000"
                if 'logged' in kwargs:
                    logged = kwargs['logged']
                else:
                    logged = True
                c.insert(hifive.plotting.plot_key(min_score=minscore, max_score=maxscore, height=(height * 0.05),
                                                  width=width, orientation='top', num_ticks=5, min_color=min_color,
                                                  mid_color=mid_color, max_color=max_color, log_display=False),
                                                  [trafo.translate(0, height * 1.05)])
                if logged:
                    label = "Log2 "
                else:
                    label = ""
                if options.datatype == 'enrichment':
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sEnrichment" % label,
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                elif options.datatype == 'raw':
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sCounts" % label,
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                else:
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sNormalized Counts" % label,
                           [text.halign.center, text.valign.bottom, text.size(-2)])
            if options.names:
                for i, chrom in enumerate(chroms):
                    c.text(width + 0.25, height - (sizes[i] + sizes[i + 1]) / 2, 'chr%s' % chrom,
                           [text.halign.left, text.valign.middle, text.size(-2)])
            c.writePDFfile(options.imagefile)
            if len(options.imagefile.split('.')) <= 1 or options.imagefile.split('.')[-1] != 'pdf':
                subprocess.call('mv %s.pdf %s' % (options.imagefile, options.imagefile), shell=True)

if __name__ == "__main__":
    main()
