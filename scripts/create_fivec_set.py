#!/usr/bin/env python

import sys

import hifive


def main():
    if len(sys.argv) < 4:
        print "Usage python create_fivec_set.py DATA_FILE OUT_FILE MIN_INTERACTIONS"
        print "DATA_FILE         File name of FiveCData h5dict to link with analysis."
        print "OUT_FILE          File name to write FiveC h5dict to."
        print "MIN_INTERACTIONS  Minimum number of interactions needed for valid fragment."
        return None
    data_fname, fivec_fname, mininteractions = sys.argv[1:4]
    mininteractions = int(mininteractions)
    fivec = hifive.FiveC(fivec_fname, 'w')
    fivec.load_data(data_fname)
    fivec.filter_fragments(mininteractions=mininteractions)
    fivec.find_distance_parameters()
    fivec.save()


if __name__ == "__main__":
    main()
