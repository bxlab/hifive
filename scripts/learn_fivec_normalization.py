#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

import sys

import hifive


fivec_fname = sys.argv[1]
fivec = hifive.analysis.FiveC(fivec_fname, 'r')
if len(sys.argv) >= 7:
    burnin_iterations = int(sys.argv[2])
    annealing_iterations = int(sys.argv[3])
    maxdistance = int(sys.argv[4])
    recalculatedistance = int(sys.argv[5])
    display = int(sys.argv[6])
    fivec.find_fragment_corrections(display=display, maxdistance=maxdistance,
                                    burnin_iterations=burnin_iterations,
                                    annealing_iterations=annealing_iterations,
                                    recalculatedistance=recalculatedistance)
else:
    fivec.find_fragment_corrections()
fivec.save()
