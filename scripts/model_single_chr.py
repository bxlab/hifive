#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

import sys

import numpy
import h5py
import mlpy

import hifive


def main():
    hic_fname, bi_fname, out_prefix, chrom, cutoff, minobservations = sys.argv[1:7]
    cutoff, minobservations = float(cutoff), int(minobservations)
    BI = hifive.bi.BI()
    BI.load(bi_fname)
    hic = hifive.analysis.HiC(hic_fname, 'r')
    bounds = find_bounds(BI, hic, chrom, cutoff)
    dynamic = find_dynamic_signal(hic, chrom, bounds, minobservations)
    coordinates = learn_model(dynamic)
    write_coordinates(coordinates, bounds, "%s_coordinates.txt" % out_prefix)
    write_model(dynamic, coordinates, bounds, chrom, cutoff, minobservations, "%s_model.hdf5" % out_prefix)
    find_model_fits(hic, dynamic, coordinates, bounds, chrom, out_prefix)
    return None


def find_model_fits(hic, dynamic, coordinates, bounds, chrom, out_prefix):
    unbinned_est, mapping = hifive.hic.binning.unbinned_cis_signal(hic, chrom, datatype='distance', arraytype='upper',
                                                          skipfiltered=True, returnmapping=True)
    unbinned_est[:, 0] = unbinned_est[:, 1]
    unbinned_est[:, 1] = 1.0
    estimated = hifive.hic.binning.bin_cis_array(hic, unbinned_est, mapping, binbounds=bounds, arraytype='upper')
    distances = numpy.sum((coordinates.reshape(-1, 1, 3) - coordinates.reshape(1, -1, 3)) ** 2.0, axis=2) ** 0.5
    distances = 1.0 / distances[numpy.triu_indices(distances.shape[0], 1)]
    where = numpy.where(estimated[:, 0] > 0)[0]
    log_dynamic = numpy.log(dynamic[:, 0] / dynamic[:, 1])
    log_distances = numpy.log(distances)
    log_estimated = numpy.log(estimated[where, 0])
    model_corr = numpy.corrcoef(log_dynamic, log_distances)[0, 1]
    estimated_corr = numpy.corrcoef(log_dynamic[where], log_estimated)[0, 1]
    output = open("%s_correlations.txt" % out_prefix, 'w')
    print >> output, "chr\tmodel_pearson\testimated_pearson"
    print >> output, "%s\t%f\t%f" % (chrom, model_corr, estimated_corr)
    output.close()
    return None


def write_model(dynamic, coordinates, bounds, chrom, cutoff, minobservations, out_fname):
    print >> sys.stderr, ("Writing model..."),
    output = h5py.File(out_fname, 'w')
    output.create_dataset(name='bin_coordinates', data=bounds)
    output.create_dataset(name='model_coordinates', data=coordinates)
    output.create_dataset(name='enrichments', data=dynamic)
    output['chromosome'] = chrom
    output['cutoff'] = cutoff
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


def find_dynamic_signal(hic, chrom, bounds, minobservations):
    unbinned, mapping = hifive.hic.binning.unbinned_cis_signal(hic, chrom, datatype='fend', arraytype='upper',
                                                               skipfiltered=True, returnmapping=True)
    mids = hic.fends['fends']['mid'][mapping]
    dynamic = hifive.hic.binning.bin_cis_signal(hic, chrom, binbounds=bounds, datatype='fend', arraytype='upper')
    hifive.hic.binning.dynamically_bin_cis_array(unbinned, mids, dynamic, bounds,
                                                 minobservations=int(minobservations))
    return dynamic


def find_bounds(bi, hic, chrom, cutoff):
    bounds = bi.find_bi_bounds(cutoff, chroms=[])
    bounds = bounds[numpy.where(bounds['chr'] == bi.chr2int[chrom])]
    binbounds = numpy.zeros((bounds.shape[0] + 1, 2), dtype=numpy.int32)
    binbounds[1:, 0] = bounds['coord'][:]
    binbounds[:-1, 1] = binbounds[1:, 0]
    chrint = hic.chr2int[chrom]
    start_fend = hic.fends['chr_indices'][chrint]
    while hic.filter[start_fend] == 0:
        start_fend += 1
    stop_fend = hic.fends['chr_indices'][chrint + 1] - 1
    while hic.filter[stop_fend] == 0:
        stop_fend -= 1
    binbounds[0, 0] = hic.fends['fends']['mid'][start_fend] - 1
    binbounds[-1, 1] = hic.fends['fends']['mid'][stop_fend] + 1
    return binbounds


if __name__ == "__main__":
    main()
