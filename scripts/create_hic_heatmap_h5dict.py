#!/usr/bin/env python

import sys

import hifive


def main():
    if len(sys.argv) < 6 and rank == 0:
        print "Usage python create_hic_heatmap_h5dict.py HIC_FILE OUT_FILE BINSIZE INCLUDE_TRANS REMOVE_DISTANCE CHROMS"
        print "HIC_FILE         File name of a HiC h5dict to pull data from."
        print "OUT_FILE         File name of heatmap h5dict to write data to."
        print "BINSIZE          Size of bins, in base pairs, to group data into."
        print "INCLUDE_TRANS    Specifies whether to find inter-chromosome interactions."
        print "REMOVE_DISTANCE  Specifies whether to remove distance-dependent portion of signal."
        print "CHROMS           Comma-separated list of chromosomes to find heatmaps for."
        print "This script is MPI compatible."
        return None
    elif len(sys.argv) < 6:
        return None
    hic_fname, out_fname, binsize, includetrans, remove_distance, chroms = sys.argv[1:7]
    chroms = sys.argv[6].split(',')
    binsize = int(binsize)
    if includetrans in ['True','true','1']:
        includetrans = True
    else:
        includetrans = False
    if remove_distance in ['True','true','1']:
        remove_distance = True
    else:
        remove_distance = False
    hic = hifive.HiC(hic_fname, 'r')
    hifive.hic_binning.write_heatmap_dict(hic, out_fname, binsize, includetrans=includetrans,
                                          remove_distance=remove_distance, chroms=chroms)


if __name__ == "__main__":
    main()
