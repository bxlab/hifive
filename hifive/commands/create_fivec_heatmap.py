#!/usr/bin/env python

import sys
import subprocess

import h5py
try:
    from pyx import *
except:
    pass

from ..fivec import FiveC
from ..plotting import plot_fivec_heatmap, plot_key


def run(args):
    if not args.image is None and args.pdf and "pyx" not in sys.modules.keys():
        print >> sys.stderr, ("-p/--pdf requires the package 'pyx'"),
        return 1
    if args.binsize > 0:
        args.arraytype = 'full'
    if args.regions is None:
        regions = []
    else:
        regions = args.regions.split(',')
        if len(regions) == 1 and regions[0] == '':
            regions = []
    for i in range(len(regions)):
        try:
            regions[i] = int(regions[i])
        except:
            print sys.stderr, ("Not all arguments in -r/--regions could be converted to integers.")
            return 1
    fivec = FiveC(args.project, 'r', silent=args.silent)
    fivec.write_heatmap(args.output, binsize=args.binsize, includetrans=args.trans, arraytype=args.arraytype,
                        datatype=args.datatype, regions=regions, dynamically_binned=args.dynamic,
                        expansion_binsize=args.expbinsize, minobservations=args.minobs, searchdistance=args.search,
                        removefailed=args.remove)
    if not args.image is None:
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
        img, minscore, maxscore = plot_fivec_heatmap(args.output, returnscale=True, silent=args.silent, **kwargs)
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
            if args.legend:
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
                c.insert(plot_key(min_score=minscore, max_score=maxscore, height=(height * 0.05),
                                  width=width, orientation='top', num_ticks=5, min_color=min_color,
                                  mid_color=mid_color, max_color=max_color, log_display=False),
                                  [trafo.translate(0, height * 1.05)])
                if logged:
                    label = "Log2 "
                else:
                    label = ""
                if args.datatype == 'enrichment':
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sEnrichment" % (label),
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                elif args.datatype == 'raw':
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sCounts" % (label),
                           [text.halign.center, text.valign.bottom, text.size(-2)])
                else:
                    c.text(width * 0.5, height * 1.1 + 0.75, "%sNormalized Counts" % (label),
                           [text.halign.center, text.valign.bottom, text.size(-2)])
            if args.names:
                for i, name in enumerate(names):
                    c.text(width + 0.25, height - (sizes[i] + sizes[i + 1]) / 2, name,
                           [text.halign.left, text.valign.middle, text.size(-2)])
            c.writePDFfile(args.image)
            if len(args.image.split('.')) <= 1 or args.image.split('.')[-1] != 'pdf':
                subprocess.call('mv %s.pdf %s' % (args.image, args.image), shell=True)
