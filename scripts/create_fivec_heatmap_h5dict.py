#!/usr/bin/env python

import sys
import optparse
import subprocess

import h5py
try:
    from pyx import *
except:
    pass

import hifive


def main():
    usage = "usage: %prog [options] <project_file> <out_file>\n\nArguments:"
    usage += "\n<project_file>  HiFive 5C project file"
    usage += "\n<out_file>      destination for 5C heatmaps file"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b", "--bin-size", dest="binsize", default=10000, metavar="BIN", type="int",
                      help="size of bins, in base pairs, to group data into [default: %default]", action="store")
    parser.add_option("-t", "--trans", dest="trans", default=False,
                      help="calculate and include trans interactions in heatmaps [default: %default]", action="store_true")
    parser.add_option("-d", "--datatype", dest="datatype", default=False, action="store",
                      help="which corrections (if any) to apply to counts (raw, distance, fragment, or enrichment) [default: %default]",
                      type="choice", choices=['raw', 'distance', 'fragment', 'enrichment'])
    parser.add_option("-r", "--regions", dest="regions", default="", metavar="REGIONS", type="string",
                      help="comma-separated list of regions to include in heatmaps [default: all regions]", action="store")
    parser.add_option("-a", "--array-type", dest="atype", action="store", default="full",
                      help="type of array layout ('full' or 'compact') to store data in (only used for unbinned data) [default: %default]",
                      type="choice", choices=["full", "compact"])
    parser.add_option("-i", "--image-file", dest="imagefile", default=None, metavar="IMAGE", type="string",
                      help="save the data as an image to this file", action="store")
    parser.add_option("-p", "--pdf", dest="pdf", default=False,
                      help="format image in PDF file [default %default]", action="store_true")
    parser.add_option("-l", "--add-legend", dest="legend", default=False,
                      help="add color scale to plot (pdf format only) [default %default]", action="store_true")
    parser.add_option("-n", "--names", dest="names", default=False,
                      help="add region labels to plot (pdf format only) [default %default]", action="store_true")
    parser.add_option("-k", "--keyword", dest="keywords", default=[], metavar="KEYWORD", type="string",
                      help="additional keyword arguments to pass function", action="append")
    parser.add_option("-q", "--quiet", dest="silent", action="store_true", default=False,
                      help="silence output messages [default: %default]")
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error('incorrect number of arguments')
    if not options.imagefile is None and options.pdf and "pyx" not in sys.modules.keys():
        parser.error("-p/--pdf requires the package 'pyx'")
    if options.binsize > 0:
        options.atype = 'full'
    fivec = hifive.FiveC(args[0], 'r', silent=options.silent)
    fivec.write_heatmap_dict(args[1], binsize=options.binsize, includetrans=options.trans, arraytype=options.atype,
                             datatype=options.datatype, regions=options.regions.split(','))
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
        if options.atype == 'full':
            img, minscore, maxscore = hifive.plotting.plot_fivec_full_heatmap_dict(args[1], returnscale=True,
                                                                                   silent=options.silent, **kwargs)
        else:
            img, minscore, maxscore = hifive.plotting.plot_fivec_compact_heatmap_dict(args[1], returnscale=True,
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
            totalx = 0
            totaly = 0
            minsize = 999999999
            sizes = [0]
            names = []
            all_regions = hm['regions'][...]
            for i, region in enumerate(all_regions['index']):
                names.append("chr%s:%i-%i" % (all_regions['chromosome'][i], all_regions['start'][i],
                                              all_regions['stop'][i]))
                if "%i.fragments" % region in hm:
                    name = "%i.fragments" % int(region)
                    sizes.append(hm[name].shape[0])
                    totalx += hm[name].shape[0]
                    totaly += hm[name].shape[0]
                    minsize = min(minsize, hm[name].shape[0])
                elif "%i.reverse_fragments" % int(region) in hm:
                    name = "%i.reverse_fragments" % int(region)
                    sizes.append(hm[name].shape[0])
                    totaly += hm[name].shape[0]
                    totalx += hm["%i.forward_fragments" % int(region)].shape[0]
                    minsize = min(minsize, hm[name].shape[0])
                elif "%i.positions" % int(region) in hm:
                    name = "%i.positions" % int(region)
                    sizes.append(hm[name].shape[0])
                    totalx += hm[name].shape[0]
                    totaly += hm[name].shape[0]
                    minsize = min(minsize, hm[name].shape[0])
            totaly += all_regions.shape[0] - 1
            totalx += all_regions.shape[0] - 1
            height = max(5.0, totaly * 0.5 / minsize)
            width = height / totaly * totalx
            if len(sizes) > 2:
                sizes[1] += 0.5
                sizes[-1] += 0.5
                if len(sizes) > 3:
                    for i in range(2, len(sizes) - 1):
                        sizes[i] += 1.0
            for i in range(1, len(sizes)):
                sizes[i] += sizes[i - 1]
            for i in range(1, len(sizes)):
                sizes[i] = sizes[i] / totaly * height
            c = canvas.canvas()
            c.insert(bitmap.bitmap(0, 0, img, width=width))
            if options.legend:
                if 'min_color' in kwargs:
                    min_color = kwargs['min_color']
                else:
                    min_color = (0, 0, 1)
                if 'mid_color' in kwargs:
                    mid_color = kwargs['mid_color']
                else:
                    mid_color = (1, 1, 1)
                if 'max_color' in kwargs:
                    max_color = kwargs['max_color']
                else:
                    max_color = (1, 0, 0)
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
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sEnrichment" % (label),
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                elif options.datatype == 'raw':
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sCounts" % (label),
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                else:
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sNormalized Counts" % (label),
                           [text.halign.center, text.valign.bottom, text.size(-2)])
            if options.names:
                for i, name in enumerate(names):
                    c.text(width + 0.25, height - (sizes[i] + sizes[i + 1]) / 2, name,
                           [text.halign.left, text.valign.middle, text.size(-2)])
            c.writePDFfile(options.imagefile)
            if len(options.imagefile.split('.')) <= 1 or options.imagefile.split('.')[-1] != 'pdf':
                subprocess.call('mv %s.pdf %s' % (options.imagefile, options.imagefile), shell=True)



if __name__ == "__main__":
    main()
