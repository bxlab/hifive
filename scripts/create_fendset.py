#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

import sys
import os

import hifive


def create_fendset(fend_fname, out_fname, genome_name=None, RE_name=None):
    fends = hifive.fend.Fend(out_fname, mode='w')
    fends.load_fends(fend_fname, genome_name=genome_name, re_name=RE_name)
    fends.fends.close()


if __name__ == "__main__":
    fend_fname, out_fname = sys.argv[1:3]
    if len(sys.argv) > 3:
        genome_name = sys.argv[3]
    else:
        genome_name = None
    if len(sys.argv) > 4:
        re_name = sys.argv[4]
    else:
        re_name = None
    create_fendset(fend_fname, out_fname, genome_name, re_name)
