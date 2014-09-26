#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys

import numpy
import h5py
import mlpy

import hifive


def main():
    if len(sys.argv) < 6:
        print "Usage python model_single_chr_binned.py HIC_FILE OUT_PREFIX BIN_SIZE MIN_OBS EXP_BINSIZE CHROM"
        print "HIC_FILE     File name of HiC h5dict to pull data from."
        print "OUT_PREFIX   File prefix for all output files of script."
        print "BIN_SIZE     Size of bins, in base pairs, to group data into."
        print "MIN_OBS      Minimum number of observations for valid dynamic bins."
        print "EXP_BINSIZE  Size of bins for additional data used for dynamic bin expansion. This may be set to zero for unbinned data."
        print "CHROM        Name of chromosome to model."
        return None
    hic_fname, out_prefix, binsize, minobservations, exp_binsize, chrom = sys.argv[1:7]
    binsize, minobservations, exp_binsize = int(binsize), int(minobservations), int(exp_binsize)
    hic = hifive.HiC(hic_fname, 'r')
    dynamic, bounds = find_dynamic_signal(hic, chrom, binsize, minobservations, exp_binsize)
    coordinates = learn_model(dynamic)
    write_coordinates(coordinates, bounds, "%s_%s_coordinates.txt" % (out_prefix, chrom))
    write_model(dynamic, coordinates, bounds, chrom, minobservations, "%s_%s_model.hdf5" % (out_prefix, chrom))
    find_model_fits(hic, dynamic, coordinates, binsize, chrom, out_prefix)
    return None


def find_model_fits(hic, dynamic, coordinates, binsize, chrom, out_prefix):
    estimated = hifive.hic_binning.bin_cis_signal(hic, chrom, binsize=binsize, datatype='distance',
                            arraytype='upper')
    estimated[:, 0] = estimated[:, 1]
    estimated[:, 1] = 1.0
    distances = numpy.sum((coordinates.reshape(-1, 1, 3) - coordinates.reshape(1, -1, 3)) ** 2.0, axis=2) ** 0.5
    distances = 1.0 / distances[numpy.triu_indices(distances.shape[0], 1)]
    where = numpy.where(estimated[:, 0] > 0)[0]
    log_dynamic = numpy.log(dynamic[:, 0] / dynamic[:, 1])
    log_distances = numpy.log(distances)
    log_estimated = numpy.log(estimated[where, 0])
    model_corr = numpy.corrcoef(log_dynamic, log_distances)[0, 1]
    estimated_corr = numpy.corrcoef(log_dynamic[where], log_estimated)[0, 1]
    output = open("%s_%s_correlations.txt" % (out_prefix, chrom), 'w')
    print >> output, "chr\tmodel_pearson\testimated_pearson"
    print >> output, "%s\t%f\t%f" % (chrom, model_corr, estimated_corr)
    output.close()
    return None


def write_model(dynamic, coordinates, bounds, chrom, minobservations, out_fname):
    print >> sys.stderr, ("Writing model..."),
    output = h5py.File(out_fname, 'w')
    output.create_dataset(name='bin_coordinates', data=bounds)
    output.create_dataset(name='model_coordinates', data=coordinates)
    output.create_dataset(name='enrichments', data=dynamic)
    output['chromosome'] = chrom
    output['minobservations'] = minobservations
    print >> sys.stderr, ("Done\n"),
    output.close()


def write_coordinates(coordinates, bounds, out_fname):
    print >> sys.stderr, ("Writing coordinates..."),
    output = open(out_fname, 'w')
    print >> output, "Start\tStop\tX\tY\tZ"
    for i in range(bounds.shape[0]):
        print >> output, "%i\t%i\t%f\t%f\t%f" % (bounds[i, 0], bounds[i, 1], coordinates[i, 0], coordinates[i, 1],
                                                 coordinates[i, 2])
    output.close()
    print >> sys.stderr, ("Done\n"),
    return None


def learn_model(binned):
    print >> sys.stderr, ("Learning model..."),
    dim = int(1 + (1 + 8 * binned.shape[0])**0.5) / 2
    data_mat = numpy.zeros((dim, dim), dtype=numpy.float32)
    indices = numpy.triu_indices(dim, 1)
    data_mat[indices] = numpy.log(binned[:, 0] / binned[:, 1])
    data_mat[indices[1], indices[0]] = data_mat[indices]
    data_mat[numpy.arange(dim), numpy.arange(dim)] = numpy.amax(data_mat[indices])
    pca_fast = mlpy.PCAFast(k=3)
    pca_fast.learn(data_mat)
    coordinates = pca_fast.transform(data_mat)
    print >> sys.stderr, ("Done\n"),
    return coordinates


def find_dynamic_signal(hic, chrom, binsize, minobservations, exp_binsize):
    dynamic, mapping2 = hifive.hic_binning.bin_cis_signal(hic, chrom, binsize=binsize, datatype='fend',
                        arraytype='upper', returnmapping=True)
    if exp_binsize == 0:
        unbinned = numpy.copy(dynamic)
        mids = (mapping2[:, 2] + mapping2[:, 3]) / 2
    else:
        unbinned, mapping = hifive.hic_binning.bin_cis_signal(hic, chrom, binsize=binsize, datatype='fend',
                            arraytype='upper', returnmapping=True)
        mids = (mapping[:, 2] + mapping[:, 3]) / 2
    hifive.hic_binning.dynamically_bin_cis_array(unbinned, mids, dynamic, mapping2[:, 2:],
                                                 minobservations=int(minobservations))
    return [dynamic, mapping2[:, 2:]]


if __name__ == "__main__":
    main()
