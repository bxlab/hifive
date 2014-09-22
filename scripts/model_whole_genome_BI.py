#!/usr/bin/env python

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
    if len(sys.argv) < 8 and rank == 0:
        print "Usage python model_whole_genome_BI.py HIC_FILE BI_FILE OUT_PREFIX CUTOFF MIN_OBS CIS_SCALING CHROMS"
        print "HIC_FILE     File name of HiC h5dict to pull data from."
        print "BI_FILE      File name of BI h5dict to find boundaries for partitioning from."
        print "OUT_PREFIX   File prefix for all output files of script."
        print "CUTOFF       Criteria for calling BI peaks."
        print "MIN_OBS      Minimum number of observations for valid dynamic bins."
        print "CIS_SCALING  Scaling factor to adjust cis interactions by prior to modeling."
        print "CHROMS       Comma-separated list of names of chromosomes to model."
        print "This script is MPI compatible."
        return None
    elif len(sys.argv) < 8:
        return None
    hic_fname, bi_fname, out_prefix, cutoff, minobservations, cis_scaling, chroms = sys.argv[1:8]
    cutoff, minobservations, cis_scaling = float(cutoff), int(minobservations), float(cis_scaling)
    BI = hifive.BI()
    BI.load(bi_fname)
    hic = hifive.HiC(hic_fname, 'r')
    chroms = chroms.split(',')
    if rank == 0:
        bounds, bound_indices = find_bounds(BI, hic, chroms, cutoff)
        for i in range(1, num_procs):
            comm.send(bounds, dest=i, tag=11)
            comm.send(bound_indices, dest=i, tag=11)
        needed_pairs = []
        for i in range(len(chroms)):
            needed_pairs.append((i,))
            for j in range(i + 1, len(chroms)):
                needed_pairs.append((i, j))
        shuffle(needed_pairs)
        worker_size = int(ceil(len(needed_pairs) / float(num_procs)))
        for i in range(1,min(num_procs, len(needed_pairs))):
            comm.send(needed_pairs[(worker_size * (i - 1)):(worker_size * i)], dest=i, tag=11)
        for i in range(min(num_procs, len(needed_pairs)), num_procs):
            comm.send([], dest=i, tag=11)
        node_pairs = needed_pairs[(min(num_procs - 1, len(needed_pairs) - 1) * worker_size):]
    else:
        bounds = comm.recv(source=0, tag=11)
        bound_indices = comm.recv(source=0, tag=11)
        node_pairs = comm.recv(source=0, tag=11)
    dynamics = {}
    for pair in node_pairs:
        dynamics[pair] = find_dynamic_signal(hic, chroms, pair, bounds, bound_indices, minobservations)
    if rank == 0:
        for i in range(1, num_procs):
            dynamics.update(comm.recv(source=i, tag=11))
        dynamic = numpy.zeros((bound_indices[-1], bound_indices[-1], 2), dtype=numpy.float32)
        for i in range(bound_indices.shape[0] - 1):
            indices = numpy.triu_indices(bound_indices[i + 1] - bound_indices[i], 1)
            dynamic[indices[0] + bound_indices[i], indices[1] + bound_indices[i], :] = dynamics[(i,)]
            for j in range(i + 1, bound_indices.shape[0] - 1):
                dynamic[bound_indices[i]:bound_indices[i + 1], bound_indices[j]:bound_indices[j + 1], :] = (
                    dynamics[(i, j)])
        del dynamics
        coordinates = learn_model(dynamic, bound_indices, cis_scaling)
        write_coordinates(coordinates, bounds, chroms, "%s_coordinates.txt" % out_prefix)
        write_model(dynamic, coordinates, bounds, bound_indices, chroms, cutoff, minobservations,
                    cis_scaling, "%s_model.hdf5" % out_prefix)
        find_model_fits(dynamic, coordinates, bound_indices, out_prefix)
    else:
        comm.send(dynamics, dest=0, tag=11)
        del dynamics
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
    model_corr = numpy.corrcoef(log_dynamic, log_distances)[0, 1]
    output = open("%s_correlations.txt" % out_prefix, 'w')
    print >> output, "model_pearson"
    print >> output, "%f" % (model_corr)
    output.close()
    return None


def write_model(dynamic, coordinates, bounds, chr_indices, chroms, cutoff, minobservations, cis_scaling,
                out_fname):
    print >> sys.stderr, ("Writing model..."),
    output = h5py.File(out_fname, 'w')
    output.create_dataset(name='bin_coordinates', data=bounds)
    output.create_dataset(name='model_coordinates', data=coordinates)
    indices = numpy.triu_indices(dynamic.shape[0], 1)
    output.create_dataset(name='enrichments', data=dynamic[indices[0], indices[1], :])
    output.create_dataset(name='chromosomes', data=chroms)
    output.create_dataset(name='chr_indices', data=chr_indices)
    output['cutoff'] = cutoff
    output['minobservations'] = minobservations
    output['cis_scaling'] = cis_scaling
    print >> sys.stderr, ("Done\n"),
    output.close()
    return None


def write_coordinates(coordinates, bounds, chroms, out_fname):
    print >> sys.stderr, ("Writing coordinates..."),
    output = open(out_fname, 'w')
    print >> output, "Chrom\tStart\tStop\tX\tY\tZ"
    for i in range(bounds.shape[0]):
        print >> output, "%s\t%i\t%i\t%f\t%f\t%f" % (chroms[bounds[i, 0]], bounds[i, 1], bounds[i, 2],
                                                     coordinates[i, 0], coordinates[i, 1], coordinates[i, 2])
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


def find_dynamic_signal(hic, chroms, pair, bounds, bound_indices, minobservations):
    start = bound_indices[pair[0]]
    stop = bound_indices[pair[0] + 1]
    chrom = chroms[pair[0]]
    mids = numpy.mean(bounds[start:stop, 1:], axis=1).astype(numpy.int32)
    if len(pair) == 1:
        binned = hifive.hic_binning.bin_cis_signal(hic, chrom, datatype='fend', arraytype='upper',
                                                   binbounds=bounds[start:stop, 1:])
        dynamic = numpy.copy(binned)
        hifive.hic_binning.dynamically_bin_cis_array(binned, mids, dynamic, minobservations=minobservations,
                                                     binbounds=bounds[start:stop, 1:])
    else:
        start2 = bound_indices[pair[1]]
        stop2 = bound_indices[pair[1] + 1]
        chrom2 = chroms[pair[1]]
        mids2 = numpy.mean(bounds[start2:stop2, 1:], axis=1).astype(numpy.int32)
        binned = hifive.hic_binning.bin_trans_signal(hic, chrom, chrom2, datatype='fend',
                                                     binbounds1=bounds[start:stop, 1:],
                                                     binbounds2=bounds[start2:stop2, 1:])
        dynamic = numpy.copy(binned)
        hifive.hic_binning.dynamically_bin_trans_array(binned, mids, mids2, dynamic, 
                                                       minobservations=minobservations,
                                                       binbounds1=bounds[start:stop, 1:],
                                                       binbounds2=bounds[start2:stop2, 1:])
    return dynamic


def find_bounds(bi, hic, chroms, cutoff):
    bounds = bi.find_bi_bounds(cutoff, chroms)
    all_bounds = numpy.zeros((bounds.shape[0] + len(chroms), 3), dtype=numpy.int32)
    bound_indices = numpy.zeros(len(chroms) + 1, dtype=numpy.int32)
    for i, chrom in enumerate(chroms):
        chrint = hic.chr2int[chrom]
        start_fend = hic.fends['chr_indices'][chrint]
        while hic.filter[start_fend] == 0:
            start_fend += 1
        stop_fend = hic.fends['chr_indices'][chrint + 1] - 1
        while hic.filter[stop_fend] == 0:
            stop_fend -= 1
        where = numpy.where(bounds['chr'][:] == chrint)[0]
        bound_indices[i + 1] = where.shape[0] + 1 + bound_indices[i]
        all_bounds[bound_indices[i], 1] = hic.fends['fends']['mid'][start_fend] - 1
        all_bounds[bound_indices[i + 1] - 1, 2] = hic.fends['fends']['mid'][stop_fend] + 1
        all_bounds[bound_indices[i]:bound_indices[i + 1], 0] = i
        all_bounds[(bound_indices[i] + 1):bound_indices[i + 1], 1] = bounds['coord'][where]
        all_bounds[bound_indices[i]:(bound_indices[i + 1] - 1), 2] = bounds['coord'][where]
    print >> sys.stderr, ('Total bins: %i\n') % (all_bounds.shape[0]),
    return [all_bounds, bound_indices]


if __name__ == "__main__":
    main()
