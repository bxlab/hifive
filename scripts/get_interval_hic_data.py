#!/usr/bin/env python

import sys
import optparse
import subprocess
from math import floor, ceil, log10

import numpy
try:
    from pyx import *
except:
    pass

import hifive


def main():
    usage = "usage: %prog [options] <project_file> <out_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive HiC project file"
    usage += "\n<out_file>      destination for interval-formatted text file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bin-size", dest="binsize", default=10000, metavar="BIN", type="int",
                      help="size of bins, in base pairs, to group data into (zero indicates unbinned) [default: %default]", action="store")
    parser.add_option("-c", "--chrom", dest="chrom", default="1", metavar="CHROM", type="string",
                      help="chromosome to pull data from", action="store")
    parser.add_option("-s", "--start", dest="start", default=None, type="int", metavar="START", action="store",
                      help="start coordinate of region to return [default: %default]")
    parser.add_option("-e", "--stop", dest="stop", default=None, type="int", metavar="STOP", action="store",
                      help="stop coordinate of region to return (zero indicates end of chromosome) [default: %default]")
    parser.add_option("-m", "--max-distance", dest="maxdist", default=0, type="int", metavar="MAXDIST",
                      help="maximum distance interactions to return (zero indicates a maximum of stop - start) [default: %default]", action="store")
    parser.add_option("-d", "--datatype", dest="datatype", default="enrichment", metavar="TYPE", type="string",
                      help="which corrections (if any) to apply to counts (raw, distance, fend, or enrichment) [default: %default]", action="store")
    parser.add_option("-i", "--image-file", dest="imagefile", default=None, metavar="IMAGE", type="string",
                      help="save the data as an image to this file", action="store")
    parser.add_option("-r", "--rotate", dest="rotate", default=False,
                      help="rotate plot 45 degrees [default %default]", action="store_true")
    parser.add_option("-p", "--pdf", dest="pdf", default=False,
                      help="format image in PDF file [default %default]", action="store_true")
    parser.add_option("-t", "--tick-marks", dest="ticks", default=False,
                      help="add tick marks and labels to plot (pdf format only) [default %default]", action="store_true")
    parser.add_option("-l", "--add-legend", dest="legend", default=False,
                      help="add color scale to plot (pdf format only) [default %default]", action="store_true")
    parser.add_option("-k", "--keyword", dest="keywords", default=[], metavar="KEYWORD", type="string",
                      help="additional keyword arguments to pass function", action="append")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help="silence output messages [default: %default]")
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error('incorrect number of arguments')
    if options.datatype not in ['raw', 'distance', 'fend', 'enrichment']:
        parser.error("-d/--datatype accepts only 'raw', 'distance', 'fend', or 'enrichment'")
    if not options.imagefile is None and options.pdf and "pyx" not in sys.modules.keys():
        parser.error("-p/--pdf requires the package 'pyx'")
    hic = hifive.HiC(args[0], 'r', silent=options.silent)
    if options.stop == 0:
        maxstop = hic.fends['fends']['stop'][hic.fends['chr_indices'][hic.chr2int[options.chrom] + 1] - 1]
    else:
        maxstop = options.stop
    if options.maxdist == 0 or options.maxdist >= (maxstop - options.start) / 2:
        arraytype = 'upper'
    else:
        arraytype = 'compact'
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
                temp[1] = temp[1].replace('__pd__','')
        kwargs[temp[0]] = temp[1]
    data, mapping = hic.cis_heatmap(chrom=options.chrom, binsize=options.binsize, start=options.start,
                                    stop=options.stop, datatype=options.datatype, arraytype=arraytype,
                                    maxdistance=options.maxdist, returnmapping=True, skipfiltered=True)
    if options.binsize == 0:
        mapping = numpy.hstack((hic.fends['fends']['start'][mapping].reshape(-1, 1),
                               hic.fends['fends']['stop'][mapping].reshape(-1, 1)))
    else:
        mapping = mapping[:, 2:]
    output = open(args[1], 'w')
    if arraytype == 'upper':
        pos = 0
        for i in range(mapping.shape[0] - 1):
            for j in range(i + 1, mapping.shape[0]):
                if data[pos, 0] > 0.0 and data[pos, 1] > 0.0:
                    print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (options.chrom, mapping[i, 0], mapping[i, 1],
                                                                        options.chrom, mapping[j, 0], mapping[j, 1],
                                                                        numpy.log2(data[pos, 0] / data[pos, 1]))
                pos += 1
    else:
        for i in range(mapping.shape[0] - 1):
            for pos in range(min(mapping.shape[0] - i - 1, data.shape[1])):
                j = i + pos + 1
                if data[i, pos, 0] > 0.0 and data[i, pos, 1] > 0.0:
                    print >> output, "chr%s\t%i\t%i\tchr%s\t%i\t%i\t%f" % (options.chrom, mapping[i, 0], mapping[i, 1],
                                                                        options.chrom, mapping[j, 0], mapping[j, 1],
                                                                        numpy.log2(data[i, pos, 0] / data[i, pos, 1]))
    output.close()
    if not options.imagefile is None:
        if options.datatype == 'enrichment':
            symmetricscaling = True
        else:
            symmetricscaling = False
        if 'symmetricscaling' in kwargs:
            symmetricscaling = kwargs['symmetricscaling']
        if arraytype == 'compact':
            if options.rotate:
                img, minscore, maxscore = hifive.plotting.plot_diagonal_from_compact_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=options.silent, **kwargs)
                offset = 2.5 / (data.shape[0] * 2 - 2)
                height = 2.5
            else:
                img, minscore, maxscore = hifive.plotting.plot_compact_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=options.silent, **kwargs)
                offset = 0.0
                height = 5.0
        else:
            if options.rotate:
                img, minscore, maxscore = hifive.plotting.plot_diagonal_from_upper_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=options.silent, **kwargs)
                offset = 2.5 / (mapping.shape[0] * 2 - 2)
                height = 2.5
            else:
                img, minscore, maxscore = hifive.plotting.plot_upper_array(data, returnscale=True,
                                          symmetricscaling=symmetricscaling, silent=options.silent, **kwargs)
                offset = 0.0
                height = 5.0
        if options.pdf:
            c = canvas.canvas()
            c1 = canvas.canvas([canvas.clip(path.rect(0, 0, 5, 5))])
            c1.insert(bitmap.bitmap(-offset, -offset, img, width=5.0))
            c.insert(c1)
            if options.ticks and options.binsize > 0:
                c.stroke(path.line(0, 0, 5.0, 0))
                xmin = (mapping[0, 0] + mapping[0, 1]) / 2
                xmax = (mapping[-1, 0] + mapping[-1, 1]) / 2
                order = int(floor(log10(xmax - xmin))) - 1
                step = int(floor((xmax - xmin) / (10.0 ** order * 5.0))) * 10 ** order
                values = numpy.arange(((xmin - 1) / step + 1) * step, (xmax / step) * step + 1, step)
                ticks = (values - float(mapping[0, 0] + mapping[0, 1]) / 2) / (mapping[-1, 0] - mapping[0, 0]) * 5.0
                for i in range(values.shape[0]):
                    c.stroke(path.line(ticks[i], 0, ticks[i], -0.25), [style.linewidth.Thin])
                    c.text(ticks[i], -0.3, "%0.2e" % values[i], [text.valign.middle, text.halign.left, text.size(-2),
                                                              trafo.rotate(-90)])
                if not options.rotate:
                    c.stroke(path.line(5.0, 0, 5.0, 5.0))
                    for i in range(values.shape[0]):
                        c.stroke(path.line(5.0, 5.0 - ticks[i], 5.25, 5.0 - ticks[i]), [style.linewidth.Thin])
                        c.text(5.3, 5.0 - ticks[i], "%0.2e" % values[i], [text.valign.middle, text.halign.left,
                                                                          text.size(-2)])
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
                c.insert(hifive.plotting.plot_key(min_score=minscore, max_score=maxscore, height=0.25, width=5.0,
                                                  orientation='top', num_ticks=5, min_color=min_color,
                                                  mid_color=mid_color, max_color=max_color,
                                                  log_display=False), [trafo.translate(0, height + 0.25)])
                if logged:
                    label = "Log2 "
                else:
                    label = ""
                if options.datatype == 'enrichment':
                    c.text(2.5, height + 0.8, "%sEnrichment" % label, [text.halign.center, text.valign.bottom,
                                                                       text.size(-2)])
                elif options.datatype == 'raw':
                    c.text(2.5, height + 0.8, "%sCounts" % label, [text.halign.center, text.valign.bottom,
                                                                   text.size(-2)])
                else:
                    c.text(2.5, height + 0.8, "%sNormalized Counts" % label, [text.halign.center, text.valign.bottom,
                                                                              text.size(-2)])

            c.writePDFfile(options.imagefile)
            if len(options.imagefile.split('.')) <= 1 or options.imagefile.split('.')[-1] != 'pdf':
                subprocess.call('mv %s.pdf %s' % (options.imagefile, options.imagefile), shell=True)

        else:
            img.save(options.imagefile, format='png')


if __name__ == "__main__":
    main()