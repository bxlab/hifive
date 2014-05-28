#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

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
