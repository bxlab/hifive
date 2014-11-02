#!/usr/bin/env python

import sys

import hifive

def main():
    if len(sys.argv) < 8:
        print "Usage python learn_fivec_normalization.py FIVEC_FILE RATE BURNIN ANNEALING MAX_DIST RECALC DISPLAY"
        print "FIVEC_FILE  File name of FiveC h5dict to analyze."
        print "RATE        Percent of gradient to use for updating parameter values."
        print "BURNIN      Number of iterations to run burn-in phase for."
        print "ANNEALING   Number of iterations to run annealing phase for."
        print "MAX_DIST    Maximum interaction distance to include in learning."
        print "RECALC      Number of iterations to wait between recalculating distance function parameters."
        print "DISPLAY     Number of iterations to wait before explicitly calculating cost and updating display."
        return None
    fivec_fname, learningrate, burnin_iterations, annealing_iterations, maxdistance, recalculate_distance, display = (
        sys.argv[1:8])
    fivec = hifive.FiveC(fivec_fname, 'r')
    fivec.find_fragment_corrections(display=int(display), maxdistance=int(maxdistance),
                                    burnin_iterations=int(burnin_iterations),
                                    annealing_iterations=int(annealing_iterations),
                                    recalculate_distance=int(recalculate_distance),
                                    learningrate=float(learningrate))
    fivec.save()


if __name__ == "__main__":
    main()
