#!/usr/bin/env python

import sys
import subprocess
from math import floor, ceil, log10

import numpy
try:
    from pyx import *
except:
    pass

from ..fivec import FiveC
from ..plotting import plot_full_array, plot_upper_array, plot_diagonal_from_upper_array, plot_key


def run(args):
    if not args.image is None and args.pdf and "pyx" not in sys.modules.keys():
        print sys.stderr, ("-p/--pdf requires the package 'pyx'"),
        return 1
    fivec = FiveC(args.project, 'r', silent=args.silent)
    if args.binsize == 0:
        arraytype = 'compact'
    else:
        arraytype = 'upper'
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
    chrom = fivec.frags['regions']['chromosome'][args.region]
    if args.region2 is None:
        temp = fivec.cis_heatmap(region=args.region, binsize=args.binsize, start=args.start,
                                 stop=args.stop, datatype=args.datatype, arraytype=arraytype,
                                 returnmapping=True, skipfiltered=True, dynamically_binned=args.dynamic,
                                 expansion_binsize=args.expbinsize, minobservations=args.minobs,
                                 searchdistance=args.search, removefailed=args.remove)
    else:
        chrom2 = fivec.frags['regions']['chromosome'][args.region2]
        temp = fivec.trans_heatmap(region1=args.region, region2=args.region2, binsize=args.binsize, start1=args.start,
                                 stop1=args.stop, start2=args.start2, stop2=args.stop2, datatype=args.datatype,
                                 arraytype='full', returnmapping=True, skipfiltered=True,
                                 dynamically_binned=args.dynamic, expansion_binsize=args.expbinsize,
                                 minobservations=args.minobs, searchdistance=args.search, removefailed=args.remove)
    output = open(args.output, 'w')
    if args.region2 is None:
        if args.binsize == 0:
            data = temp[0]
            xmapping = temp[1][:, :2]
            ymapping = temp[2][:, :2]
            all_data = []
            for i in range(xmapping.shape[0]):
                for j in range(ymapping.shape[0]):
                    if data[i, j, 0] <= 0.0:
                        continue
                    if xmapping[i, 0] < ymapping[j, 0]:
                        all_data.append((xmapping[i, 0], xmapping[i, 1], ymapping[j, 0], ymapping[j, 1],
                                         numpy.log2(data[i, j, 0] / data[i, j, 1])))
                    else:
                        all_data.append((ymapping[j, 0], ymapping[j, 1], xmapping[i, 0], xmapping[i, 1],
                                         numpy.log2(data[i, j, 0] / data[i, j, 1])))
            all_data = numpy.array(all_data, dtype=numpy.dtype([('start1', numpy.int32), ('stop1', numpy.int32),
                                                                ('start2', numpy.int32), ('stop2', numpy.int32),
                                                                ('value', numpy.float32)]))
            order = numpy.lexsort((all_data['start2'], all_data['start1']))
            for i in order:
                print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (chrom, all_data['start1'][i],
                    all_data['stop1'][i], chrom, all_data['start2'][i], all_data['stop2'][i], all_data['value'][i])
        else:
            data = temp[0]
            mapping = temp[1][:, :2]
            pos = 0
            for i in range(mapping.shape[0] - 1):
                for j in range(i + 1, mapping.shape[0]):
                    if data[pos, 0] > 0.0 and data[pos, 1] > 0.0:
                        print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (chrom, mapping[i, 0], mapping[i, 1],
                            chrom, mapping[j, 0], mapping[j, 1], numpy.log2(data[pos, 0] / data[pos, 1]))
                    pos += 1
    else:
        data, mapping1, mapping2 = temp
        for i in range(mapping1.shape[0]):
            for j in range(mapping2.shape[0]):
                if data[i, j, 0] <= 0.0:
                    continue
                print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (chrom, mapping1[i, 0], mapping1[i, 1], chrom2,
                                                                       mapping2[j, 0], mapping2[j, 1],
                                                                       numpy.log2(data[i, j, 0] / data[i, j, 1]))
    output.close()
    if not args.image is None:
        if args.datatype == 'enrichment':
            symmetricscaling = True
        else:
            symmetricscaling = False
        if 'symmetricscaling' in kwargs:
            symmetricscaling = kwargs['symmetricscaling']
        if not args.region2 is None:
            img, minscore, maxscore = plot_full_array(data, returnscale=True, symmetricscaling=symmetricscaling,
                                                      silent=args.silent, **kwargs)
            offset = 0.0
            width = 5.0
            height = width / mapping1.shape[0] * mapping2.shape[0]
        elif arraytype == 'compact':
            img, minscore, maxscore = plot_full_array(data, returnscale=True, symmetricscaling=symmetricscaling,
                                                      silent=args.silent, **kwargs)
            offset = 0.0
            width = 5.0
            height = width / xmapping.shape[0] * ymapping.shape[0]
        else:
            if args.rotate:
                img, minscore, maxscore = plot_diagonal_from_upper_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=args.silent, **kwargs)
                offset = 2.5 / (mapping.shape[0] * 2 - 2)
                height = 2.5
                width = 5.0
            else:
                img, minscore, maxscore = plot_upper_array(data, returnscale=True, symmetricscaling=symmetricscaling,
                                                           silent=args.silent, **kwargs)
                offset = 0.0
                height = width = 5.0
        if args.pdf:
            c = canvas.canvas()
            c1 = canvas.canvas([canvas.clip(path.rect(0, 0, width, height))])
            c1.insert(bitmap.bitmap(-offset, -offset, img, width=width))
            c.insert(c1)
            if args.region2 is None:
                if args.ticks and args.binsize > 0:
                    c.stroke(path.line(0, 0, width, 0))
                    xmin = (mapping[0, 0] + mapping[0, 1]) / 2
                    xmax = (mapping[-1, 0] + mapping[-1, 1]) / 2
                    order = int(floor(log10(xmax - xmin))) - 1
                    step = int(floor((xmax - xmin) / (10.0 ** order * width))) * 10 ** order
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
                            c.stroke(path.line(width, height - ticks[i], width + 0.25, height - ticks[i]),
                                     [style.linewidth.Thin])
                            c.text(width + 0.3, height - ticks[i], "%0.2e" % values[i],
                                   [text.valign.middle, text.halign.left, text.size(-2)])
            elif args.ticks:
                c.stroke(path.line(0, 0, width, 0))
                xmin = (mapping1[0, 0] + mapping1[0, 1]) / 2
                xmax = (mapping1[-1, 0] + mapping1[-1, 1]) / 2
                order = int(floor(log10(xmax - xmin))) - 1
                step = int(floor((xmax - xmin) / (10.0 ** order * width))) * 10 ** order
                values = numpy.arange(((xmin - 1) / step + 1) * step, (xmax / step) * step + 1, step)
                ticks = (values - float(mapping1[0, 0] + mapping1[0, 1]) / 2) / (mapping1[-1, 0] -
                                                                               mapping1[0, 0]) * width
                for i in range(values.shape[0]):
                    c.stroke(path.line(ticks[i], 0, ticks[i], -0.25), [style.linewidth.Thin])
                    c.text(ticks[i], -0.3, "%0.2e" % values[i],
                           [text.valign.middle, text.halign.left, text.size(-2), trafo.rotate(-90)])

                c.stroke(path.line(width, 0, width, height))
                xmin = (mapping2[0, 0] + mapping2[0, 1]) / 2
                xmax = (mapping2[-1, 0] + mapping2[-1, 1]) / 2
                order = int(floor(log10(xmax - xmin))) - 1
                step = int(floor((xmax - xmin) / (10.0 ** order * width))) * 10 ** order
                values = numpy.arange(((xmin - 1) / step + 1) * step, (xmax / step) * step + 1, step)
                ticks = (values - float(mapping2[0, 0] + mapping2[0, 1]) / 2) / (mapping2[-1, 0] -
                                                                               mapping2[0, 0]) * height
                for i in range(values.shape[0]):
                    c.stroke(path.line(width, height - ticks[i], width + 0.25, height - ticks[i]),
                             [style.linewidth.Thin])
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
                c.insert(plot_key(min_score=minscore, max_score=maxscore, height=0.25, width=width,
                                  orientation='top', num_ticks=5, min_color=min_color,
                                  mid_color=mid_color, max_color=max_color,
                                  log_display=False), [trafo.translate(0, height + 0.25)])
                if logged:
                    label = "Log2 "
                else:
                    label = ""
                if args.datatype == 'enrichment':
                    c.text(width * 0.5, height + 0.8, "%sEnrichment" % label, [text.halign.center, text.valign.bottom,
                                                                       text.size(-2)])
                elif args.datatype == 'raw':
                    c.text(width + 0.5, height + 0.8, "%sCounts" % label, [text.halign.center, text.valign.bottom,
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
