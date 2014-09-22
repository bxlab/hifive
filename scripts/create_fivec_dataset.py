#!/usr/bin/env python

import sys
import os

import hifive


def main():
    if len(sys.argv) < 4:
        print "Usage: python create_fivec_dataset.py FRAG_FILE DATA_FILE_1[,...,DATA_FILE_N] OUT_FILE"
        print "FRAG_FILE                     File name of Fragment h5dict to link with data."
        print "DATA_FILE_1[,...DATA_FILE_N]  A comma-separated list of either count file names or BAM file prefices."
        print "OUT_FILE                      File name to write FiveCData h5dict to."
        return None
    frag_fname, data_fnames, out_fname = sys.argv[1:4]
    data = hifive.FiveCData(out_fname, 'w')
    if data_fnames.lower().endswith('counts'):
        data_fnames = data_fnames.split(',')
        data.load_data_from_counts(frag_fname, data_fnames)
    else:
        data_fnames = data_fnames.split(',')
        data.load_data_from_bam(frag_fname, data_fnames)


if __name__ == "__main__":
    main()
