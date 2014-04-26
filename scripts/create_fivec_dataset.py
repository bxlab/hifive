#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys
import os

import hifive


def create_fivec_dataset(frag_fname, data_fnames, out_fname):
    data = hifive.fivec.data.FiveCData(out_fname, 'w')
    if data_fnames.lower().endswith('counts'):
        data_fnames = data_fnames.split(',')
        data.load_data_from_counts(frag_fname, data_fnames)
    else:
        data_fnames = data_fnames.split(',')
        data.load_data_from_bam(frag_fname, data_fnames)


if __name__ == "__main__":
    frag_fname, data_fnames, out_fname = sys.argv[1:4]
    create_fivec_dataset(frag_fname, data_fnames, out_fname)
