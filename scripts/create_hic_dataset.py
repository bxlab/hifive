#!/usr/bin/env python

import sys
import os

import hifive


def main():
    if len(sys.argv) < 5:
        print "Usage: python create_hic_dataset.py FEND_FILE DATA_FILE_1[,...,DATA_FILE_N] OUT_FILE MAX_INSERT"
        print "FEND_FILE                     File name of Fend h5dict to link with data."
        print "DATA_FILE_1[,...DATA_FILE_N]  A comma-separated list of either BAM file prefices, raw coordinate read pairs or HiCPipe-compatible MAT files."
        print "OUT_FILE                      File name to write HiCData h5dict to."
        print "MAX_INSERT                    Integer specifying the maximum distance sum from each mapped end to restriction site."
        return None
    fend_fname, data_fnames, out_fname, maxinsert = sys.argv[1:5]
    maxinsert = int(maxinsert)
    data = hifive.HiCData(out_fname, 'w')
    if data_fnames.endswith('mat'):
        data.load_data_from_mat(fend_fname, data_fnames, maxinsert)
    elif data_fnames.endswith('raw'):
        data_fnames = data_fnames.split(',')
        data.load_data_from_raw(fend_fname, data_fnames, maxinsert)
    else:
        data_fnames = data_fnames.split(',')
        data.load_data_from_bam(fend_fname, data_fnames, maxinsert)


if __name__ == "__main__":
    main()
