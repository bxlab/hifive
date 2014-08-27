#!/usr/bin/env python
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys
import os

import numpy
import h5py


def main():
    in_fnames, out_fname = sys.argv[1:3]
    data = []
    in_fnames = in_fnames.split(',')
    for fname in in_fnames:
        data.append(h5py.File(fname, 'r'))
    basedir = '/' + '/'.join(in_fnames[0].split('/')[:-1]).strip('/') + '/'
    fendfilename = data[0]['/'].attrs['fendfilename']
    if fendfilename[:2] == './':
        fendfilename = fendfilename[2:]
    fendfilename = "%s/%s" % (basedir, '/'.join(fendfilename.split('/')[fendfilename.count('../'):]))
    fends = h5py.File(fendfilename, 'r')
    maxinsert = data[0]['/'].attrs['maxinsert']
    output = h5py.File(out_fname, 'w')
    output.attrs['maxinsert'] = maxinsert
    output.attrs['fendfilename'] = os.path.relpath(fendfilename, out_fname)
    num_fends = fends['fends'].shape[0]
    temp_data = []
    all_cis_pairs = numpy.zeros(0, dtype=numpy.int64)
    for i in range(len(data)):
        temp_data.append(data[i]['cis_data'][...])
        all_cis_pairs = numpy.hstack((all_cis_pairs, temp_data[i][:, 0].astype(numpy.int64) * num_fends +
                                                     temp_data[i][:, 1].astype(numpy.int64)))
    all_cis_pairs = numpy.unique(all_cis_pairs)
    cis_data = numpy.zeros((all_cis_pairs.shape[0], 3), dtype=numpy.int32)
    cis_data[:, 0] = all_cis_pairs / num_fends
    cis_data[:, 1] = all_cis_pairs % num_fends
    for i in range(len(data)):
        indices = temp_data[i][:, 0].astype(numpy.int64) * num_fends + temp_data[i][:, 1].astype(numpy.int64)
        cis_data[numpy.searchsorted(all_cis_pairs, indices), 2] += temp_data[i][:, 2]
    del all_cis_pairs
    del temp_data
    temp_data = []
    all_trans_pairs = numpy.zeros(0, dtype=numpy.int64)
    for i in range(len(data)):
        temp_data.append(data[i]['trans_data'][...])
        all_trans_pairs = numpy.hstack((all_trans_pairs, temp_data[i][:, 0].astype(numpy.int64) * num_fends +
                                                     temp_data[i][:, 1].astype(numpy.int64)))
    all_trans_pairs = numpy.unique(all_trans_pairs)
    trans_data = numpy.zeros((all_trans_pairs.shape[0], 3), dtype=numpy.int32)
    trans_data[:, 0] = all_trans_pairs / num_fends
    trans_data[:, 1] = all_trans_pairs % num_fends
    for i in range(len(data)):
        indices = temp_data[i][:, 0].astype(numpy.int64) * num_fends + temp_data[i][:, 1].astype(numpy.int64)
        trans_data[numpy.searchsorted(all_trans_pairs, indices), 2] += temp_data[i][:, 2]
    del all_trans_pairs
    del temp_data
    cis_indices = numpy.r_[0, numpy.bincount(cis_data[:, 0], minlength=num_fends)].astype(numpy.int32)
    trans_indices = numpy.r_[0, numpy.bincount(trans_data[:, 0], minlength=num_fends)].astype(numpy.int32)
    for i in range(1, cis_indices.shape[0]):
        cis_indices[i] += cis_indices[i - 1]
        trans_indices[i] += trans_indices[i - 1]
    if cis_data.shape[0] > 0:
        output.create_dataset('cis_data', data=cis_data)
        output.create_dataset('cis_indices', data=cis_indices)
    if trans_data.shape[0] > 0:
        output.create_dataset('trans_data', data=trans_data)
        output.create_dataset('trans_indices', data=trans_indices)
    output.close()
    for i in range(len(data)):
        data[i].close()
    return None


if __name__ == "__main__":
    main()