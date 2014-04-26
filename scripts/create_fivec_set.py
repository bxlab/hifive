#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys

import hifive


def create_fivec_set(data_fname, fivec_fname, mininteractions=10):
    fivec = hifive.analysis.FiveC(fivec_fname, 'w')
    fivec.load_data(data_fname)
    fivec.filter_fragments(mininteractions=mininteractions)
    fivec.find_distance_parameters()
    fivec.save()


if __name__ == "__main__":
    data_fname, fivec_fname = sys.argv[1:3]
    if len(sys.argv) > 3:
        mininteractions = int(sys.argv[3])
        create_fivec_set(data_fname, fivec_fname, mininteractions)
    else:
        create_fivec_set(data_fname, fivec_fname)
