#!/usr/bin/env python

import sys
import subprocess
from math import floor, ceil, log10

import numpy
try:
    from pyx import *
except:
    pass

from ..hic import HiC
from ..plotting import plot_compact_array, plot_upper_array, plot_diagonal_from_compact_array
from ..plotting import plot_full_array, plot_diagonal_from_upper_array, plot_key


def run(args):
    if not args.image is None and args.pdf and "pyx" not in sys.modules.keys():
        parser.error("-p/--pdf requires the package 'pyx'")
    hic = HiC(args.project, 'r', silent=args.silent)
    if 'binned' in hic.fends['/'].attrs and hic.fends['/'].attrs['binned'] is not None:
        binned = True
        chr_indices = 'bin_indices'
        fends = 'bins'
    else:
        binned = False
        chr_indices = 'chr_indices'
        fends = 'fends'
    if args.stop == 0 or args.stop is None:
        maxstop = hic.fends[fends]['stop'][hic.fends[chr_indices][hic.chr2int[args.chrom] + 1] - 1]
    else:
        maxstop = args.stop
    if args.stop is None:
        args.stop = maxstop
    if args.start is None:
        args.start = hic.fends[fends]['start'][hic.fends[chr_indices][hic.chr2int[args.chrom]]]
    if not args.chrom2 is None:
        if args.stop2 == 0 or args.stop2 is None:
            maxstop2 = hic.fends[fends]['stop'][hic.fends[chr_indices][hic.chr2int[args.chrom2] + 1] - 1]
        else:
            maxstop2 = args.stop2
        if args.stop2 is None:
            args.stop2 = maxstop2
        if args.start2 is None:
            args.start2 = hic.fends[fends]['start'][hic.fends[chr_indices][hic.chr2int[args.chrom2]]]
    else:
        if args.maxdist is None:
            args.maxdist = 0
        if args.maxdist == 0 or args.maxdist >= (maxstop - args.start) / 2:
            arraytype = 'upper'
        else:
            arraytype = 'compact'
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
                # strip off extra characters introduced by galaxy into color format
                temp[1] = temp[1].replace('__pd__','')
        kwargs[temp[0]] = temp[1]
    if not args.chrom2 is None:
        data, mapping1, mapping2 = hic.trans_heatmap(chrom1=args.chrom, chrom2=args.chrom2, binsize=args.binsize,
                                                     start1=args.start, stop1=args.stop, start2=args.start2,
                                                     stop2=args.stop2, datatype=args.datatype,
                                                     maxdistance=args.maxdist, returnmapping=True, skipfiltered=True,
                                                     dynamically_binned=args.dynamic,
                                                     expansion_binsize=args.expbinsize, minobservations=args.minobs,
                                                     searchdistance=args.search, removefailed=args.remove)
    else:
        data, mapping = hic.cis_heatmap(chrom=args.chrom, binsize=args.binsize, start=args.start,
                                        stop=args.stop, datatype=args.datatype, arraytype=arraytype,
                                        maxdistance=args.maxdist, returnmapping=True, skipfiltered=True,
                                        dynamically_binned=args.dynamic, expansion_binsize=args.expbinsize,
                                        minobservations=args.minobs, searchdistance=args.search,
                                        removefailed=args.remove)
    output = open(args.output, 'w')
    if args.matrix:
        if args.chrom2 is None:
            diag = int(hic.binned is not None)
            if arraytype == 'upper':
                temp = numpy.zeros((mapping.shape[0], mapping.shape[0]), dtype=numpy.float64)
                indices = numpy.triu_indices(mapping.shape[0], 1 - diag)
                where = numpy.where(data[:, 1] > 0)[0]
                temp[indices[0][where], indices[1][where]] = data[where, 0] / data[where, 1]
                temp[indices[1][where], indices[0][where]] += data[where, 0] / data[where, 1]
            else:
                temp = numpy.zeros((mapping.shape[0], mapping.shape[0]), dtype=numpy.float64)
                for i in range(mapping.shape[0] - 1 + diag):
                    where = numpy.where(data[i, :, 1] > 0)[0]
                    temp[i, where + i + 1 - diag] = data[i, where, 0] / data[i, where, 1]
                indices = numpy.triu_indices(mapping.shape[0], 1 - diag)
                temp[indices[1], indices[0]] += temp[indices]
        else:
            temp = numpy.zeros((data.shape[0], data.shape[1]), dtype=numpy.float64)
            where = numpy.where(data[:, :, 1] > 0)
            temp[where] = data[where[0], where[1], 0] / data[where[0], where[1], 1]
        if args.datatype == 'raw':
            for i in range(temp.shape[0]):
                tempout = []
                for j in range(temp.shape[1]):
                    tempout.append("%i" % temp[i, j])
                print >> output, '\t'.join(tempout)
        else:
            for i in range(temp.shape[0]):
                tempout = []
                for j in range(temp.shape[1]):
                    tempout.append("%0.6f" % temp[i, j])
                print >> output, '\t'.join(tempout)
    else:
        if args.chrom2 is None:
            diag = binned
            if arraytype == 'upper':
                pos = 0
                for i in range(mapping.shape[0] - 1 + diag):
                    for j in range(i + 1 - diag, mapping.shape[0]):
                        if data[pos, 0] > 0.0 and data[pos, 1] > 0.0:
                            print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (args.chrom, mapping[i, 0],
                                                                           mapping[i, 1], args.chrom,
                                                                           mapping[j, 0], mapping[j, 1],
                                                                           numpy.log2(data[pos, 0] / data[pos, 1]))
                        pos += 1
            else:
                for i in range(mapping.shape[0] - 1 + diag):
                    for pos in range(min(mapping.shape[0] - i - 1 + diag, data.shape[1])):
                        j = i + pos + 1 - diag
                        if data[i, pos, 0] > 0.0 and data[i, pos, 1] > 0.0:
                            print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (args.chrom, mapping[i, 0],
                                                                            mapping[i, 1], args.chrom,
                                                                            mapping[j, 0], mapping[j, 1],
                                                                            numpy.log2(data[i, pos, 0] /
                                                                            data[i, pos, 1]))
        else:
            for i in range(mapping1.shape[0]):
                for j in range(mapping2.shape[0]):
                    if data[i, j, 0] > 0.0 and data[i, j, 1] > 0.0:
                        print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (args.chrom,
                                                                        mapping1[i, 0], mapping1[i, 1],
                                                                        args.chrom2, mapping2[j, 0], mapping2[j, 1],
                                                                        numpy.log2(data[i, j, 0] / data[i, j, 1]))
    output.close()
    if not args.image is None:
        width = max(5.0, (args.stop - args.start) / 1000000.)
        if args.datatype == 'enrichment':
            symmetricscaling = True
        else:
            symmetricscaling = False
        if 'symmetricscaling' in kwargs:
            symmetricscaling = kwargs['symmetricscaling']
        if not args.chrom2 is None:
            img, minscore, maxscore = plot_full_array(data, returnscale=True, symmetricscaling=symmetricscaling,
                                                      silent=args.silent, **kwargs)
            offset = 0.0
            height = (width / data.shape[0]) * data.shape[1]
        elif arraytype == 'compact':
            if args.rotate:
                img, minscore, maxscore = plot_diagonal_from_compact_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=args.silent,
                                          diagonal_included=diag, **kwargs)
                offset = width / 2. / (data.shape[0] * 2 - 1 + diag)
                height = width / (data.shape[0] * 2.0 - 2) * data.shape[1]
            else:
                img, minscore, maxscore = plot_compact_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=args.silent,
                                          diagonal_included=diag, **kwargs)
                offset = 0.0
                height = width
        else:
            if args.rotate:
                img, minscore, maxscore = plot_diagonal_from_upper_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=args.silent,
                                          diagonal_included=diag, **kwargs)
                offset = width / 2. / (mapping.shape[0] * 2 - 1 + diag)
                height = width / 2.
            else:
                img, minscore, maxscore = plot_upper_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=args.silent,
                                          diagonal_included=diag, **kwargs)
                offset = 0.0
                height = width
        if args.pdf:
            c = canvas.canvas()
            if args.chrom2 is None:
                c1 = canvas.canvas([canvas.clip(path.rect(0, 0, width, height))])
                c1.insert(bitmap.bitmap(-offset, -offset, img, width=width))
            else:
                c1 = canvas.canvas([canvas.clip(path.rect(0, 0, width, height))])
                c1.insert(bitmap.bitmap(-offset, -offset, img, width=width))
            c.insert(c1)
            if args.ticks and args.binsize > 0:
                if args.chrom2 is None:
                    c.stroke(path.line(0, 0, width, 0))
                    xmin = (mapping[0, 0] + mapping[0, 1]) / 2
                    xmax = (mapping[-1, 0] + mapping[-1, 1]) / 2
                    #order = int(floor(log10(xmax - xmin))) - 1
                    #step = int(floor((xmax - xmin) / (10.0 ** order))) * 10 ** order
                    
                    order = int(floor(log10((xmax - xmin) / (width * 2.0))))
                    step = int(floor((xmax - xmin) / (width * 2.0) / (10.0 ** order))) * 10 ** order
                    values = numpy.arange(((xmin - 1) / step + 1) * step, (xmax / step) * step + 1, step)
                    ticks = (values - float(mapping[0, 0] + mapping[0, 1]) / 2) / (mapping[-1, 0] -
                                                                                   mapping[0, 0]) * width
                    for i in range(values.shape[0]):
                        c.stroke(path.line(ticks[i], 0, ticks[i], -0.25), [style.linewidth.Thin])
                        c.text(ticks[i], -0.3, "%0.2e" % values[i],
                               [text.valign.middle, text.halign.left, text.size(-2), trafo.rotate(-90)])
                    if not args.rotate:
                        c.stroke(path.line(width, 0, width, height))
                        for i in range(values.shape[0]):
                            c.stroke(path.line(width, height - ticks[i], width + 0.25, height - ticks[i]), [style.linewidth.Thin])
                            c.text(width + 0.3, height - ticks[i], "%0.2e" % values[i], [text.valign.middle, text.halign.left,
                                                                              text.size(-2)])
                else:
                    c.stroke(path.line(0, 0, width, 0))
                    xmin = (mapping1[0, 0] + mapping1[0, 1]) / 2
                    xmax = (mapping1[-1, 0] + mapping1[-1, 1]) / 2
                    order = int(floor(log10((xmax - xmin) / (width * 2.0))))
                    step = int(floor((xmax - xmin) / (width * 2.0) / (10.0 ** order))) * 10 ** order
                    values = numpy.arange(((xmin - 1) / step + 1) * step, (xmax / step) * step + 1, step)
                    ticks = (values - float(mapping1[0, 0] + mapping1[0, 1]) / 2) / (mapping1[-1, 0] -
                                                                                     mapping1[0, 0]) * width
                    for i in range(values.shape[0]):
                        c.stroke(path.line(ticks[i], 0, ticks[i], -0.25), [style.linewidth.Thin])
                        c.text(ticks[i], -0.3, "%0.2e" % values[i],
                               [text.valign.middle, text.halign.left, text.size(-2), trafo.rotate(-90)])
                    c.stroke(path.line(0, 0, width, 0))
                    xmin = (mapping2[0, 0] + mapping2[0, 1]) / 2
                    xmax = (mapping2[-1, 0] + mapping2[-1, 1]) / 2
                    order = int(floor(log10((xmax - xmin) / (width * 2.0))))
                    step = int(floor((xmax - xmin) / (width * 2.0) / (10.0 ** order))) * 10 ** order
                    values = numpy.arange(((xmin - 1) / step + 1) * step, (xmax / step) * step + 1, step)
                    ticks = (values - float(mapping2[0, 0] + mapping2[0, 1]) / 2) / (mapping2[-1, 0] -
                                                                                     mapping2[0, 0]) * height
                    for i in range(values.shape[0]):
                        c.stroke(path.line(width, height - ticks[i], width + 0.25, height - ticks[i]), [style.linewidth.Thin])
                        c.text(width + 0.3, height - ticks[i], "%0.2e" % values[i],
                               [text.valign.middle, text.halign.left, text.size(-2)])
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
                c.insert(plot_key(min_score=minscore, max_score=maxscore, height=0.25, width=min(5., width),
                                  orientation='top', num_ticks=5, min_color=min_color,
                                  mid_color=mid_color, max_color=max_color,
                                  log_display=False), [trafo.translate(width * 0.5 - min(2.5, width * 0.5), height + 0.25)])
                if logged:
                    label = "Log2 "
                else:
                    label = ""
                if args.datatype == 'enrichment':
                    c.text(width * 0.5, height + 0.8, "%sEnrichment" % label, [text.halign.center, text.valign.bottom,
                                                                       text.size(-2)])
                elif args.datatype == 'raw':
                    c.text(width * 0.5, height + 0.8, "%sCounts" % label, [text.halign.center, text.valign.bottom,
                                                                   text.size(-2)])
                else:
                    c.text(width * 0.5, height + 0.8, "%sNormalized Counts" % label,
                           [text.halign.center, text.valign.bottom, text.size(-2)])

            c.writePDFfile(args.image)
            if len(args.image.split('.')) <= 1 or args.image.split('.')[-1] != 'pdf':
                subprocess.call('mv %s.pdf %s' % (args.image, args.image), shell=True)

        else:
            img_format = args.image.split('.')[-1].upper()
            if img_format not in ['PNG', 'TIF', 'JPG', 'JPEG']:
                img_format = 'PNG'
            img.save(args.image, img_format)
