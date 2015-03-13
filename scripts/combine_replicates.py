#!/usr/bin/env python

import sys
import os

import numpy
import h5py


def main():
    if len(sys.argv) < 3:
        print "Usage python combine_replicates.py REP_FILE_1,REP_FILE_2[,...,REP_FILE_N] OUT_FILE"
        print "REP_FILE_1,REP_FILE_2  A comma-separated list of HiCData h5dict files."
        print "OUT_FILE               File name to write new HiCData h5dict to."
        return None
    in_fnames, out_fname = sys.argv[1:3]
    out_fname = os.path.abspath(out_fname)
    data = []
    in_fnames = in_fnames.split(',')
    history = ""
    for fname in in_fnames:
        data.append(h5py.File(fname, 'r'))
        history += data[-1]['/'].attrs['history']
    fendfilename = data[0]['/'].attrs['fendfilename']
    if fendfilename[:2] == './':
        fendfilename = fendfilename[2:]
    parent_count = fendfilename.count('../')
    fendfilename = '/'.join(os.path.abspath(in_fnames[0]).split('/')[:-(1 + parent_count)] +
                            fendfilename.lstrip('/').split('/')[parent_count:])
    fends = h5py.File(fendfilename, 'r')
    maxinsert = data[0]['/'].attrs['maxinsert']
    output = h5py.File(out_fname, 'w')
    output.attrs['maxinsert'] = maxinsert
    output.attrs['fendfilename'] = "%s/%s" % (os.path.relpath(os.path.dirname(fendfilename),
                                              os.path.dirname(out_fname)), os.path.basename(fendfilename))
    num_fends = fends['fends'].shape[0]
    temp_data = []
    all_cis_pairs = numpy.zeros(0, dtype=numpy.int64)
    for i in range(len(data)):
        temp = data[i]['cis_data'][...].astype(numpy.int64)
        temp_data.append([temp[:, 0] * num_fends + temp[:, 1], temp[:, 2]])
        del temp
        all_cis_pairs = numpy.hstack((all_cis_pairs, temp_data[i][0]))
    all_cis_pairs = numpy.unique(all_cis_pairs)
    cis_data = numpy.zeros((all_cis_pairs.shape[0], 3), dtype=numpy.int32)
    cis_data[:, 0] = all_cis_pairs / num_fends
    cis_data[:, 1] = all_cis_pairs % num_fends
    for i in range(len(data)):
        cis_data[numpy.searchsorted(all_cis_pairs, temp_data[i][0]), 2] += temp_data[i][1]
    del all_cis_pairs
    del temp_data
    temp_data = []
    all_trans_pairs = numpy.zeros(0, dtype=numpy.int64)
    for i in range(len(data)):
        temp = data[i]['trans_data'][...].astype(numpy.int64)
        temp_data.append([temp[:, 0] * num_fends + temp[:, 1], temp[:, 2]])
        all_trans_pairs = numpy.hstack((all_trans_pairs, temp_data[i][0]))
    all_trans_pairs = numpy.unique(all_trans_pairs)
    trans_data = numpy.zeros((all_trans_pairs.shape[0], 3), dtype=numpy.int32)
    trans_data[:, 0] = all_trans_pairs / num_fends
    trans_data[:, 1] = all_trans_pairs % num_fends
    for i in range(len(data)):
        trans_data[numpy.searchsorted(all_trans_pairs, temp_data[i][0]), 2] += temp_data[i][1]
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
    history += "Combined replicates - Succes\n"
    output.attrs['history'] = history
    output.close()
    for i in range(len(data)):
        data[i].close()
    return None


if __name__ == "__main__":
    main()