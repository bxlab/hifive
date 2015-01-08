#!/usr/bin/env python

import sys
import os

import hifive


def main():
    if len(sys.argv) < 3:
        print "Usage: python create_fendset.py FEND_FILE OUT_FILE [GENOME RE]"
        print "FEND_FILE  File containing restriction fragment data in HiCPipe-compatible or BED format."
        print "OUT_FILE   File name to write Fend h5dict to."
        print "GENOME     Name of genome."
        print "RE         Name of restriction enzyme."
        return None
    fend_fname, out_fname = sys.argv[1:3]
    if len(sys.argv) > 4:
        genome_name = sys.argv[3]
        re_name = sys.argv[4]
    else:
        genome_name = None
        re_name = None
    fends = hifive.Fend(out_fname, mode='w')
    fends.load_fends(fend_fname, genome_name=genome_name, re_name=re_name)
    fends.fends.close()


if __name__ == "__main__":
    main()
