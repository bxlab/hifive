#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys

import numpy
import h5py
import mlpy
from math import ceil
from random import shuffle
try:
    from mpi4py import MPI
except:
    pass
from PIL import Image

import hifive


def main():
    if 'mpi4py' in sys.modules.keys():
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        num_procs = comm.Get_size()
    else:
        comm = None
        rank = 0
        num_procs = 1
    if len(sys.argv) < 7 and rank == 0:
        print "Usage python model_whole_genome_binned.py HIC_FILE BI_FILE OUT_PREFIX BIN_SIZE MIN_OBS CIS_SCALING CHROMS"
        print "HIC_FILE     File name of HiC h5dict to pull data from."
        print "OUT_PREFIX   File prefix for all output files of script."
        print "BIN_SIZE    Size of bins, in base pairs, to group data into."
        print "MIN_OBS      Minimum number of observations for valid dynamic bins."
        print "CIS_SCALING  Scaling factor to adjust cis interactions by prior to modeling."
        print "CHROMS       Comma-separated list of names of chromosomes to model."
        print "This script is MPI compatible."
        return None
    elif len(sys.argv) < 7:
        return None
    hic_fname, out_prefix, binsize, minobservations, cis_scaling, chroms = sys.argv[1:7]
    binsize, minobservations, cis_scaling = int(binsize), int(minobservations), float(cis_scaling)
    hic = hifive.analysis.HiC(hic_fname, 'r')
    chroms = chroms.split(',')
    if rank == 0:
        needed_pairs = []
        for i in range(len(chroms)):
            needed_pairs.append((i,))
            for j in range(i + 1, len(chroms)):
                needed_pairs.append((i, j))
        shuffle(needed_pairs)
        starts = numpy.zeros((len(chroms), 3), dtype=numpy.int32)
        for i, chrom in enumerate(chroms):
            chrint = hic.chr2int[chrom]
            start_fend = hic.fends['chr_indices'][chrint]
            stop_fend = hic.fends['chr_indices'][chrint + 1] - 1
            while start_fend < stop_fend and hic.filter[start_fend] == 0:
                start_fend += 1
            while stop_fend > start_fend and hic.filter[stop_fend] == 0:
                stop_fend -= 1
            start = (hic.fends['fends']['mid'][start_fend] / binsize) * binsize
            stop = ((hic.fends['fends']['mid'][stop_fend] - 1) / binsize + 1) * binsize
            starts[i, 0] = start
            starts[i, 1] = stop
            starts[i, 2] = (stop - start) / binsize
        bound_indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int32)
        for i in range(bound_indices.shape[0] - 1):
            bound_indices[i + 1] = starts[i, 2] + bound_indices[i]
        worker_size = int(ceil(len(needed_pairs) / float(num_procs)))
        for i in range(1,min(num_procs, len(needed_pairs))):
            comm.send(needed_pairs[(worker_size * (i - 1)):(worker_size * i)], dest=i, tag=11)
            comm.send(starts, dest=i, tag=11)
        for i in range(min(num_procs, len(needed_pairs)), num_procs):
            comm.send([], dest=i, tag=11)
            comm.send(starts, dest=i, tag=11)
        node_pairs = needed_pairs[(min(num_procs - 1, len(needed_pairs) - 1) * worker_size):]
        dynamics = numpy.zeros((bound_indices[-1], bound_indices[-1], 2), dtype=numpy.float32)
    else:
        node_pairs = comm.recv(source=0, tag=11)
        starts = comm.recv(source=0, tag=11)
        dynamics = {}
    for pair in node_pairs:
        if rank == 0:
            if len(pair) == 1:
                dynamics[bound_indices[pair[0]]:bound_indices[pair[0] + 1],
                         bound_indices[pair[0]]:bound_indices[pair[0] + 1],
                         :] = find_dynamic_signal(hic, chroms, starts, pair, binsize, minobservations)
            else:
                dynamics[bound_indices[pair[0]]:bound_indices[pair[0] + 1],
                         bound_indices[pair[1]]:bound_indices[pair[1] + 1],
                         :] = find_dynamic_signal(hic, chroms, starts, pair, binsize, minobservations)
        else:
            dynamics[pair] = find_dynamic_signal(hic, chroms, starts, pair, binsize, minobservations)
    if rank == 0:
        for i in range(1, num_procs):
            keys = comm.recv(source=i, tag=11)
            for j in range(len(keys)):
                if len(keys[j]) == 1:
                    dynamics[bound_indices[keys[j][0]]:bound_indices[keys[j][0] + 1],
                             bound_indices[keys[j][0]]:bound_indices[keys[j][0] + 1], :] = comm.recv(source=i, tag=11)
                else:
                    dynamics[bound_indices[keys[j][0]]:bound_indices[keys[j][0] + 1],
                             bound_indices[keys[j][1]]:bound_indices[keys[j][1] + 1], :] = comm.recv(source=i, tag=11)
        coordinates = learn_model(dynamics, bound_indices, cis_scaling)
        write_coordinates(coordinates, starts, bound_indices, binsize, chroms, "%s_coordinates.txt" % out_prefix)
        write_model(dynamics, coordinates, starts, bound_indices, chroms, binsize, minobservations,
                    cis_scaling, "%s_model.hdf5" % out_prefix)
        find_model_fits(dynamics, coordinates, bound_indices, out_prefix)
    else:
        keys = dynamics.keys()
        comm.send(keys, dest=0, tag=11)
        for key in keys:
            comm.send(dynamics[key], dest=0, tag=11)
            del dynamics[key]
    return None


def find_model_fits(dynamic, coordinates, indices, out_prefix):
    distances = numpy.sum((coordinates.reshape(-1, 1, 3) - coordinates.reshape(1, -1, 3)) ** 2.0, axis=2) ** 0.5
    triu_indices = numpy.triu_indices(distances.shape[0], 1)
    distances[triu_indices] = 1.0 / distances[triu_indices]
    for i in range(indices.shape[0] - 1):
        dynamic[indices[i]:indices[i + 1], indices[i]:indices[i + 1], 0] = 0
    where = numpy.where(dynamic[:, :, 0] > 0)
    log_dynamic = numpy.log(dynamic[where[0], where[1], 0] / dynamic[where[0], where[1], 1])
    log_distances = numpy.log(distances[where])
    model_corr = numpy.corrcoef(log_dynamic, log_distances)
    hifive.plotting.plot_full_array(dynamic, symmetricscaling=False).save('dynamic.png')
    new_distances = numpy.ones(dynamic.shape, dtype=numpy.float32)
    new_distances[:, :, 0] = distances
    new_distances[where[0], where[1], 0] = 0
    new_distances[where[0], where[1], 1] = 0
    hifive.plotting.plot_full_array(new_distances, symmetricscaling=False).save('distance.png')
    print model_corr
    model_corr = model_corr[0, 1]
    output = open("%s_correlations.txt" % out_prefix, 'w')
    print >> output, "model_pearson"
    print >> output, "%f" % (model_corr)
    output.close()
    return None


def write_model(dynamic, coordinates, starts, bound_indices, chroms, binsize, minobservations,
                cis_scaling, out_fname):
    print >> sys.stderr, ("Writing model..."),
    output = h5py.File(out_fname, 'w')
    bounds = numpy.zeros((dynamic.shape[0], 2), dtype=numpy.int32)
    for i in range(bound_indices.shape[0] - 1):
        bounds[bound_indices[i]:bound_indices[i + 1], 0] = (numpy.arange(bound_indices[i + 1] - bound_indices[i]) *
                                                            binsize + starts[i, 0])
        bounds[bound_indices[i]:bound_indices[i + 1], 1] = (numpy.arange(1, bound_indices[i + 1] - bound_indices[i] +
                                                            1) * binsize + starts[i, 0])
    output.create_dataset(name='bin_coordinates', data=bounds)
    output.create_dataset(name='model_coordinates', data=coordinates)
    indices = numpy.triu_indices(dynamic.shape[0], 1)
    output.create_dataset(name='enrichments', data=dynamic[indices[0], indices[1], :])
    output.create_dataset(name='chromosomes', data=chroms)
    output.create_dataset(name='chr_indices', data=bound_indices)
    output['minobservations'] = minobservations
    output['cis_scaling'] = cis_scaling
    print >> sys.stderr, ("Done\n"),
    output.close()
    return None


def write_coordinates(coordinates, starts, bound_indices, binsize, chroms, out_fname):
    print >> sys.stderr, ("Writing coordinates..."),
    output = open(out_fname, 'w')
    print >> output, "Chrom\tStart\tStop\tX\tY\tZ"
    for i in range(bound_indices.shape[0] - 1):
        for j in range(bound_indices[i], bound_indices[i + 1]):
            print >> output, "%s\t%i\t%i\t%f\t%f\t%f" % (chroms[i], starts[i, 0] + binsize * (j - bound_indices[i]),
                                                         starts[i, 0] + binsize * (j - bound_indices[i] + 1),
                                                         coordinates[j, 0], coordinates[j, 1], coordinates[j, 2])
    output.close()
    print >> sys.stderr, ("Done\n"),
    return None


def learn_model(binned, indices, cis_scaling):
    print >> sys.stderr, ("Learning model..."),
    dim = binned.shape[0]
    data_mat = numpy.zeros((dim, dim), dtype=numpy.float32)
    triu_indices = numpy.triu_indices(dim, 1)
    data_mat[triu_indices] = numpy.log(binned[triu_indices[0], triu_indices[1], 0] /
                                      binned[triu_indices[0], triu_indices[1], 1])
    data_mat[triu_indices[1], triu_indices[0]] = data_mat[triu_indices]
    for i in range(indices.shape[0] - 1):
        triu_indices = numpy.triu_indices(indices[i + 1] - indices[i], 1)
        diag = numpy.arange(indices[i], indices[i + 1])
        data_mat[diag, diag] = numpy.amax(data_mat[triu_indices[0] + indices[i], triu_indices[1] + indices[i]])
        data_mat[indices[i]:indices[i + 1], indices[i]:indices[i + 1]] += numpy.log(cis_scaling)
    pca_fast = mlpy.PCAFast(k=3)
    pca_fast.learn(data_mat)
    coordinates = pca_fast.transform(data_mat)
    print >> sys.stderr, ("Done\n"),
    return coordinates


def find_dynamic_signal(hic, chroms, starts, pair, binsize, minobservations):
    start, stop = starts[pair[0], :2]
    chrom = chroms[pair[0]]
    if len(pair) == 1:
        binned, mapping = hifive.hic.binning.bin_cis_signal(hic, chrom, datatype='fend', arraytype='upper',
                                                   binsize=binsize, start=start, stop=stop, returnmapping=True)
        mids = numpy.sum(mapping[:, 2:], axis=1).astype(numpy.int32) / 2
        dynamic = numpy.copy(binned)
        hifive.hic.binning.dynamically_bin_cis_array(binned, mids, dynamic, minobservations=minobservations,
                                                     binbounds=mapping[:, 2:])
        full = numpy.zeros((starts[pair[0], 2], starts[pair[0], 2], 2), dtype=numpy.float32)
        indices = numpy.triu_indices(starts[pair[0], 2], 1)
        full[indices[0], indices[1], :] = dynamic
        full[indices[1], indices[0], :] = dynamic
        return full
    else:
        chrom2 = chroms[pair[1]]
        start2, stop2 = starts[pair[1], :2]
        binned, map1, map2 = hifive.hic.binning.bin_trans_signal(hic, chrom, chrom2, datatype='fend',
                                                     start1=start, stop1=stop, binsize=binsize,
                                                     start2=start2, stop2=stop2, returnmapping=True)
        bounds1 = map1[:, 2:].astype(numpy.int32)
        mids1 = (numpy.sum(map1[:, 2:], axis=1) / 2).astype(numpy.int32)
        bounds2 = map2[:, 2:].astype(numpy.int32)
        mids2 = (numpy.sum(map2[:, 2:], axis=1) / 2).astype(numpy.int32)
        dynamic = numpy.copy(binned)
        hifive.hic.binning.dynamically_bin_trans_array(binned, mids1, mids2, dynamic, 
                                                       binbounds1=bounds1,
                                                       binbounds2=bounds2,
                                                       minobservations=minobservations)
        return dynamic


if __name__ == "__main__":
    main()
