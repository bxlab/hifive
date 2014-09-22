#!/usr/bin/env python

import sys

import hifive


def main():
    if len(sys.argv) < 5:
        print "Usage python learn_fivec_normalization_express.py FIVEC_FILE ITERATIONS REMOVE_DIST RECALC"
        print "FIVEC_FILE    File name of FiveC h5dict to analyze."
        print "ITERATIONS    Number of iterations to run learning for."
        print "REMOVE_DIST   Specifies whether to remove distance-dependent portion of the signal prior to learning."
        print "RECALC        Number of iterations to wait between recalculating distance function parameters."
        return None
    fivec_fname, iterations, removedistance, recalculatedistance = sys.argv[1:5]
    fivec = hifive.FiveC(fivec_fname, 'r')
    if removedistance in ['1', 'true', 'True', 'TRUE']:
        removedistance = True
    else:
        removedistance = False
    fivec.find_express_fragment_corrections(iterations=int(iterations), remove_distance=removedistance,
                                             recalculatedistance=int(recalculatedistance))
    fivec.save()


if __name__ == "__main__":
    main()
