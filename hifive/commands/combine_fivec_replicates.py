#!/usr/bin/env python

import sys
import os

import numpy
import h5py


def run(args):
    in_fnames = args.replicates
    out_fname = os.path.abspath(args.output)
    silent = args.silent
    data = []
    history = ""
    for fname in in_fnames:
        data.append(h5py.File(fname, 'r'))
        history += data[-1]['/'].attrs['history']
    fragfilename = data[0]['/'].attrs['fragfilename']
    if fragfilename[:2] == './':
        fragfilename = fragfilename[2:]
    parent_count = fragfilename.count('../')
    fragfilename = '/'.join(os.path.abspath(in_fnames[0]).split('/')[:-(1 + parent_count)] +
                            fragfilename.lstrip('/').split('/')[parent_count:])
    frags = h5py.File(fragfilename, 'r')
    output = h5py.File(out_fname, 'w')
    output.attrs['fragfilename'] = "%s/%s" % (os.path.relpath(os.path.dirname(fragfilename),
                                              os.path.dirname(out_fname)), os.path.basename(fragfilename))
    num_frags = frags['fragments'].shape[0]
    temp_data = []
    all_cis_pairs = numpy.zeros(0, dtype=numpy.int64)
    for i in range(len(data)):
        if not silent:
            print >> sys.stderr, ('\r%s\rLoading %s cis data...') % (' ' * 80, in_fnames[i].split('/')[-1]),
        temp = data[i]['cis_data'][...].astype(numpy.int64)
        temp_data.append([temp[:, 0] * num_frags + temp[:, 1], temp[:, 2]])
        del temp
        all_cis_pairs = numpy.hstack((all_cis_pairs, temp_data[i][0]))
    if not silent:
        print >> sys.stderr, ('\r%s\rCompiling cis data...') % (' ' * 80),
    all_cis_pairs = numpy.unique(all_cis_pairs)
    cis_data = numpy.zeros((all_cis_pairs.shape[0], 3), dtype=numpy.int32)
    cis_data[:, 0] = all_cis_pairs / num_frags
    cis_data[:, 1] = all_cis_pairs % num_frags
    for i in range(len(data)):
        cis_data[numpy.searchsorted(all_cis_pairs, temp_data[i][0]), 2] += temp_data[i][1]
    del all_cis_pairs
    del temp_data
    temp_data = []
    all_trans_pairs = numpy.zeros(0, dtype=numpy.int64)
    for i in range(len(data)):
        if not silent:
            print >> sys.stderr, ('\r%s\rLoading %s trans data...') % (' ' * 80, in_fnames[i].split('/')[-1]),
        temp = data[i]['trans_data'][...].astype(numpy.int64)
        temp_data.append([temp[:, 0] * num_frags + temp[:, 1], temp[:, 2]])
        all_trans_pairs = numpy.hstack((all_trans_pairs, temp_data[i][0]))
    if not silent:
        print >> sys.stderr, ('\r%s\rCompiling trans data...') % (' ' * 80),
    all_trans_pairs = numpy.unique(all_trans_pairs)
    trans_data = numpy.zeros((all_trans_pairs.shape[0], 3), dtype=numpy.int32)
    trans_data[:, 0] = all_trans_pairs / num_frags
    trans_data[:, 1] = all_trans_pairs % num_frags
    for i in range(len(data)):
        trans_data[numpy.searchsorted(all_trans_pairs, temp_data[i][0]), 2] += temp_data[i][1]
    del all_trans_pairs
    del temp_data
    if not silent:
        print >> sys.stderr, ('\r%s\rWriting data to file...') % (' ' * 80),
    cis_indices = numpy.r_[0, numpy.bincount(cis_data[:, 0], minlength=num_frags)].astype(numpy.int32)
    trans_indices = numpy.r_[0, numpy.bincount(trans_data[:, 0], minlength=num_frags)].astype(numpy.int32)
    for i in range(1, cis_indices.shape[0]):
        cis_indices[i] += cis_indices[i - 1]
        trans_indices[i] += trans_indices[i - 1]
    if cis_data.shape[0] > 0:
        output.create_dataset('cis_data', data=cis_data)
        output.create_dataset('cis_indices', data=cis_indices)
    if trans_data.shape[0] > 0:
        output.create_dataset('trans_data', data=trans_data)
        output.create_dataset('trans_indices', data=trans_indices)
    history += "Combined replicates - Succes\n"
    output.attrs['history'] = history
    output.close()
    for i in range(len(data)):
        data[i].close()
    if not silent:
        print >> sys.stderr, ('\r%s\rCombine 5C replicates... Done\n') % (' ' * 80),
    return None
