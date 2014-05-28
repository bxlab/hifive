#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

import sys

import numpy

import hifive


def export_data_to_mat(data_fname, mat_fname):
    data = hifive.hic.data.HiCData(data_fname)
    data.export_to_mat(mat_fname)


if __name__ == "__main__":
    data_fname, mat_fname = sys.argv[1:3]
    export_data_to_mat(data_fname, mat_fname)
