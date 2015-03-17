#!/usr/bin/env python

import sys
import os

import numpy
import h5py


def run(args):
    in_fname1 = args.replicate1
    in_fname2 = args.replicate2
    out_fname = args.output
    silent = args.silent
    history = ""
    infile1 = h5py.File(in_fname1, 'r')
    history += infile1['/'].attrs['history']
    infile2 = h5py.File(in_fname2, 'r')
    history += infile2['/'].attrs['history']
    fendfilename = infile1['/'].attrs['fendfilename']
    if fendfilename[:2] == './':
        fendfilename = fendfilename[2:]
    parent_count = fendfilename.count('../')
    fendfilename = '/'.join(os.path.abspath(in_fname1).split('/')[:-(1 + parent_count)] +
                            fendfilename.lstrip('/').split('/')[parent_count:])
    maxinsert = infile1['/'].attrs['maxinsert']
    output = h5py.File(out_fname, 'w')
    output.attrs['maxinsert'] = maxinsert
    output.attrs['fendfilename'] = "%s/%s" % (os.path.relpath(os.path.dirname(fendfilename),
                                              os.path.dirname(out_fname)), os.path.basename(fendfilename))
    fends = h5py.File(fendfilename, 'r')
    num_fends = fends['fends'].shape[0]
    fends.close()
    data1 = infile1['cis_data'][...]
    data2 = infile2['cis_data'][...]
    indices1 = infile1['cis_indices'][...]
    indices2 = infile2['cis_indices'][...]
    cindices = numpy.zeros(indices1.shape[0], dtype=numpy.int64)
    fend2 = numpy.zeros(data1.shape[0] + data2.shape[0], dtype=numpy.int32)
    pos = 0
    for i in range(num_fends):
        if not silent:
            print >> sys.stderr, ("\rCis phase 1, %08i of %08i") % (i, num_fends),
        f1 = indices1[i + 1] - indices1[i]
        f2 = indices2[i + 1] - indices2[i]
        if f1 > 0 and f2 > 0:
            temp = numpy.unique(numpy.hstack((data1[indices1[i]:indices1[i + 1], 1],
                                              data2[indices2[i]:indices2[i + 1], 1])))
            fend2[pos:(pos + temp.shape[0])] = temp
            pos += temp.shape[0]
        elif f1 > 0:
            fend2[pos:(pos + f1)] = data1[indices1[i]:indices1[i + 1], 1]
            pos += f1
        elif f2 > 0:
            fend2[pos:(pos + f2)] = data2[indices2[i]:indices2[i + 1], 1]
            pos += f2
        cindices[i + 1] = pos
    cis_data = numpy.zeros((cindices[-1], 3), dtype=numpy.int32)
    cis_data[:, 1] = fend2[:cindices[-1]]
    del fend2
    for i in range(num_fends):
        if not silent:
            print >> sys.stderr, ("\rCis phase 2, %08i of %08i") % (i, num_fends),
        cis_data[cindices[i]:cindices[i + 1], 0] = i
        f1 = indices1[i + 1] - indices1[i]
        f2 = indices2[i + 1] - indices2[i]
        if f1 > 0:
            pos = numpy.searchsorted(cis_data[cindices[i]:cindices[i + 1], 1], data1[indices1[i]:indices1[i + 1], 1])
            if pos.shape[0] > 0: 
                cis_data[pos + cindices[i], 2] += data1[indices1[i]:indices1[i + 1], 2]
        if f2 > 0:
            pos = numpy.searchsorted(cis_data[cindices[i]:cindices[i + 1], 1], data2[indices2[i]:indices2[i + 1], 1])
            if pos.shape[0] > 0: 
                cis_data[pos + cindices[i], 2] += data2[indices2[i]:indices2[i + 1], 2]
    del data1
    del data2
    del indices1
    del indices2
    output.create_dataset(name='cis_data', data=cis_data)
    output.create_dataset(name='cis_indices', data=cindices)
    count = 0
    valid = 0
    for i in range(cindices.shape[0] - 1):
        if cindices[i + 1] > cindices[i]:
            count += 1
            if cis_data[cindices[i], 0] == cis_data[cindices[i + 1] - 1, 0]:
                valid += 1
    if not silent:
        print >> sys.stderr, ("\r%s\r%i out of %i cis fends validated\n") % (" " * 80, valid, count),
    del cis_data
    del cindices
    data1 = infile1['trans_data'][...]
    data2 = infile2['trans_data'][...]
    indices1 = infile1['trans_indices'][...]
    indices2 = infile2['trans_indices'][...]
    tindices = numpy.zeros(indices1.shape[0], dtype=numpy.int64)
    fend2 = numpy.zeros(data1.shape[0] + data2.shape[0], dtype=numpy.int32)
    pos = 0
    for i in range(num_fends):
        if not silent:
            print >> sys.stderr, ("\rTrans phase 1, %08i of %08i") % (i, num_fends),
        f1 = indices1[i + 1] - indices1[i]
        f2 = indices2[i + 1] - indices2[i]
        if f1 > 0 and f2 > 0:
            temp = numpy.unique(numpy.hstack((data1[indices1[i]:indices1[i + 1], 1],
                                              data2[indices2[i]:indices2[i + 1], 1])))
            fend2[pos:(pos + temp.shape[0])] = temp
            pos += temp.shape[0]
        elif f1 > 0:
            fend2[pos:(pos + f1)] = data1[indices1[i]:indices1[i + 1], 1]
            pos += f1
        elif f2 > 0:
            fend2[pos:(pos + f2)] = data2[indices2[i]:indices2[i + 1], 1]
            pos += f2
        tindices[i + 1] = pos
    trans_data = numpy.zeros((tindices[-1], 3), dtype=numpy.int32)
    trans_data[:, 1] = fend2[:tindices[-1]]
    del fend2
    for i in range(num_fends):
        if not silent:
            print >> sys.stderr, ("\rTrans phase 2, %08i of %08i") % (i, num_fends),
        trans_data[tindices[i]:tindices[i + 1], 0] = i
        f1 = indices1[i + 1] - indices1[i]
        f2 = indices2[i + 1] - indices2[i]
        if f1 > 0:
            pos = numpy.searchsorted(trans_data[tindices[i]:tindices[i + 1], 1], data1[indices1[i]:indices1[i + 1], 1])
            if pos.shape[0] > 0: 
                trans_data[pos + tindices[i], 2] += data1[indices1[i]:indices1[i + 1], 2]
        if f2 > 0:
            pos = numpy.searchsorted(trans_data[tindices[i]:tindices[i + 1], 1], data2[indices2[i]:indices2[i + 1], 1])
            if pos.shape[0] > 0: 
                trans_data[pos + tindices[i], 2] += data2[indices2[i]:indices2[i + 1], 2]
    del data1
    del data2
    del indices1
    del indices2
    output.create_dataset(name='trans_data', data=trans_data)
    output.create_dataset(name='trans_indices', data=tindices)
    output.attrs['history'] = history
    count = 0
    valid = 0
    for i in range(tindices.shape[0] - 1):
        if tindices[i + 1] > tindices[i]:
            count += 1
            if trans_data[tindices[i], 0] == trans_data[tindices[i + 1] - 1, 0]:
                valid += 1
    if not silent:
        print >> sys.stderr, ("\r%s\r%i out of %i trans fends validate\n") % (" " * 80, valid, count),
    del trans_data
    del tindices
    output.close()
    if not silent:
        print >> sys.stderr, ("\r%s\rCombining HiC replicates... Done\n") % (" " * 80),
    return None
