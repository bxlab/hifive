#!/usr/bin/env python

import sys
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

from ..hic import HiC
from ..plotting import plot_hic_heatmap, plot_key


def run(args):
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    if not args.image is None and args.pdf and "pyx" not in sys.modules.keys():
        if rank == 0:
            print >> sys.stderr, ("-p/--pdf requires the package 'pyx'"),
        sys.exit(1)
    if args.chroms is None:
        chroms = []
    else:
        chroms = args.chroms.split(',')
        if len(chroms) == 1 and chroms[0] == '':
            chroms = []
    hic = HiC(args.project, 'r', silent=args.silent)
    hic.write_heatmap(args.output, binsize=args.binsize, includetrans=args.trans,
                      datatype=args.datatype, chroms=chroms, dynamically_binned=args.dynamic,
                      expansion_binsize=args.expbinsize, minobservations=args.minobs,
                      searchdistance=args.search, removefailed=args.remove, format=args.format)
    if rank > 0:
        sys.exit(0)
    if not args.image is None:
        if args.format == 'txt':
            if rank == 0:
                print >> sys.stderr, ("Plotting is only available for non-txt formats.\n"),
            return None
        kwargs = {}
        for arg in args.keywords:
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
                        # strip off extra characters introduced by galaxy into color format
                        temp[1] = temp[1].replace('__pd__','')
            kwargs[temp[0]] = temp[1]
        if 'symmetricscaling' not in kwargs:
            if args.datatype == 'enrichment':
                kwargs['symmetricscaling'] = True
            else:
                kwargs['symmetricscaling'] = False
        img, minscore, maxscore = plot_hic_heatmap(args.output, returnscale=True, silent=args.silent,
                                                   format=args.format, **kwargs)
        if not args.pdf:
            img_format = args.image.split('.')[-1].upper()
            if img_format not in ['PNG', 'TIF', 'JPG', 'JPEG']:
                img_format = 'PNG'
            img.save(args.image, img_format)
        else:
            unit.set(defaultunit="cm")
            text.set(mode="latex")
            text.preamble(r"\usepackage{times}")
            text.preamble(r"\usepackage{sansmath}")
            text.preamble(r"\sansmath")
            text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
            hm = h5py.File(args.output, 'r')
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
            if args.legend:
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
                c.insert(plot_key(min_score=minscore, max_score=maxscore, height=(height * 0.05),
                                  width=width, orientation='top', num_ticks=5, min_color=min_color,
                                  mid_color=mid_color, max_color=max_color, log_display=False),
                                  [trafo.translate(0, height * 1.05)])
                if logged:
                    label = "Log2 "
                else:
                    label = ""
                if args.datatype == 'enrichment':
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sEnrichment" % label,
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                elif args.datatype == 'raw':
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sCounts" % label,
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                else:
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sNormalized Counts" % label,
                           [text.halign.center, text.valign.bottom, text.size(-2)])
            if args.names:
                for i, chrom in enumerate(chroms):
                    c.text(width + 0.25, height - (sizes[i] + sizes[i + 1]) / 2, 'chr%s' % chrom,
                           [text.halign.left, text.valign.middle, text.size(-2)])
            c.writePDFfile(args.image)
            if len(args.image.split('.')) <= 1 or args.image.split('.')[-1] != 'pdf':
                subprocess.call('mv %s.pdf %s' % (args.image, args.image), shell=True)
