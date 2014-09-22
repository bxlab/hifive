#!/usr/bin/env python

import sys

import numpy

import hifive


def main():
    if len(sys.argv) < 3:
        print "Usage python data2mat.py DATA_FILE OUT_FILE"
        print "DATA_FILE  File name of HiCData h5dict."
        print "OUT_FILE   File name to write HiCPipe-compatible MAT-formatted data to."
        return None
    data_fname, mat_fname = sys.argv[1:3]
    data = hifive.HiCData(data_fname)
    data.export_to_mat(mat_fname)


if __name__ == "__main__":
    main()
