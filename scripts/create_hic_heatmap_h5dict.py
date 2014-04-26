#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys

import hifive


def create_hic_heatmap_h5dict(hic_fname, out_fname, binsize, includetrans=True, removedistance=False, chroms=[]):
    hic = hifive.analysis.HiC(hic_fname, 'r')
    hifive.hic.binning.write_heatmap_dict(hic, out_fname, binsize, includetrans=includetrans,
                                          removedistance=removedistance, chroms=chroms)


if __name__ == "__main__":
    hic_fname, out_fname, binsize, includetrans, removedistance = sys.argv[1:6]
    if len(sys.argv) > 6:
        chroms = sys.argv[6].split(',')
    else:
        chroms = []
    binsize = int(binsize)
    if includetrans in ['True','true','1']:
        includetrans = True
    else:
        includetrans = False
    if removedistance in ['True','true','1']:
        removedistance = True
    else:
        removedistance = False
    create_hic_heatmap_h5dict(hic_fname, out_fname, binsize, includetrans, removedistance, chroms)
