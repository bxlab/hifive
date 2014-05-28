#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

import sys
import os

import hifive


def create_dataset(fend_fname, data_fnames, out_fname, maxinsert):
    data = hifive.hic.data.HiCData(out_fname, 'w')
    if data_fnames.endswith('mat'):
        data.load_data_from_mat(fend_fname, data_fnames, maxinsert)
    elif data_fnames.endswith('raw'):
        data_fnames = data_fnames.split(',')
        data.load_data_from_raw(fend_fname, data_fnames, maxinsert)
    else:
        data_fnames = data_fnames.split(',')
        data.load_data_from_bam(fend_fname, data_fnames, maxinsert)


if __name__ == "__main__":
    fend_fname, data_fnames, out_fname, maxinsert = sys.argv[1:5]
    create_dataset(fend_fname, data_fnames, out_fname, int(maxinsert))
