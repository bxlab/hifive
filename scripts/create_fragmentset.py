#!/usr/bin/env python

import sys
import os

import hifive


def main():
    if len(sys.argv) < 3:
        print "Usage: python create_fragmentset.py BED_FILE OUT_FILE [GENOME RE]"
        print "BED_FILE  File containing 5C targeted fragments and primer names in BED format."
        print "OUT_FILE  File name to write Fragment h5dict to."
        print "GENOME    Name of genome."
        print "RE        Name of restriction enzyme."
        return None
    fragment_fname, out_fname = sys.argv[1:3]
    if len(sys.argv) > 4:
        genome_name = sys.argv[3]
        re_name = sys.argv[4]
    else:
        genome_name = None
        re_name = None
    fragments = hifive.Fragment(out_fname, mode='w')
    fragments.load_fragments(fragment_fname, genome_name=genome_name, re_name=re_name)
    fragments.fragments.close()
    return None


if __name__ == "__main__":
    main()
