#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys

import hifive


fivec_fname = sys.argv[1]
fivec = hifive.analysis.FiveC(fivec_fname, 'r')
if len(sys.argv) >= 5:
    iterations = int(sys.argv[2])
    ignoredistance = sys.argv[3]
    if ignoredistance in ['1', 'true', 'True', 'TRUE']:
        ignoredistance = True
    else:
        ignoredistance = False
    recalculatedistance = int(sys.argv[4])
    fivec.find_express_fragment_corrections(iterations=iterations, ignoredistance=ignoredistance,
                                             recalculatedistance=recalculatedistance)
else:
    fivec.find_express_fragment_corrections()
fivec.save()
