#!/usr/bin/env python

import sys
import os

import h5py


def run(args):
    target_fname = os.path.abspath(args.target)
    source_fname = os.path.abspath(args.source)
    infile = h5py.File(target_fname, 'r')
    if args.type in ['hic-data', '5c-data']:
        if args.type == 'hic-data':
            if 'fendfilename' in infile['/'].attrs:
                secondary_fname = infile['/'].attrs['fendfilename']
            else:
                print >> sys.stderr, ('Unable to find fend filename in data file. The fend file link will not be updated.\n'),
                secondary_fname = None
        elif args.type == '5c-data':
            if 'fragfilename' in infile['/'].attrs:
                secondary_fname = infile['/'].attrs['fragfilename']
            else:
                print >> sys.stderr, ('Unable to find fragment filename in data file. The fragment file link will not be updated.\n'),
                secondary_fname = None
    else:
        secondary_fname = None
    print 'source: %s' % source_fname
    print 'target: %s' % target_fname
    print 'secondary:', secondary_fname
    if secondary_fname is not None:
        if secondary_fname[:2] == './':
            secondary_fname = secondary_fname[2:]
        parent_count = secondary_fname.count('../')
        secondary_fname = '/'.join(target_fname.split('/')[:-(1 + parent_count)] +
                            secondary_fname.lstrip('/').split('/')[parent_count:])
        secondary_fname = "%s/%s" % (os.path.relpath(os.path.dirname(secondary_fname),
                                  os.path.dirname(source_fname)), os.path.basename(secondary_fname))
    target_fname = "%s/%s" % (os.path.relpath(os.path.dirname(target_fname),
                              os.path.dirname(source_fname)), os.path.basename(target_fname))
    print 'target: %s' % target_fname
    print 'secondary:', secondary_fname
    outfile = h5py.File(source_fname, 'a')
    if args.type == 'fends':
        outfile['/'].attrs['fendfilename'] = target_fname
    elif args.type == 'fragments':
        outfile['/'].attrs['fragfilename'] = target_fname
    if args.type == 'hic-data':
        if secondary_fname is not None:
            outfile['/'].attrs['fendfilename'] = secondary_fname
        outfile['/'].attrs['datafilename'] = target_fname
    if args.type == '5c-data':
        if secondary_fname is not None:
            outfile['/'].attrs['fragfilename'] = secondary_fname
        outfile['/'].attrs['datafilename'] = target_fname
    infile.close()
    outfile.close()
