#!/usr/bin/env python

import sys
import datetime
import argparse

import numpy
import hifive
import h5py

def main():
    parser = generate_parser()
    args = parser.parse_args()
    hic = hifive.HiC(args.input)
    hic.filter.fill(1)
    chromosomes = hic.fends['chromosomes'][...]
    chroms = []
    if args.chroms == '':
        for i in range(1, 23):
            if numpy.where(chromosomes == str(i))[0].shape[0] > 0:
                chroms.append(str(i))
        if numpy.where(chromosomes == 'X')[0].shape[0] > 0:
            chroms.append('X')
    else:
        for chrom in args.chroms.split(','):
            if numpy.where((chromosomes == chrom) | (chromosomes == chrom.strip('chr')))[0].shape[0] > 0:
                chroms.append(chrom)
    args.chroms = numpy.array(chroms)
    length_array = numpy.zeros(len(args.chroms), dtype=numpy.int32)
    temp = hic.fends['chrom_sizes'][...]
    for i, chrom in enumerate(args.chroms):
        if chrom in hic.chr2int:
            length_array[i] = temp[hic.chr2int[chrom]]
        else:
            length_array[i] = temp[hic.chr2int[chrom.strip('chr')]]
    order = numpy.argsort(length_array)[::-1]
    length_array = length_array[order]
    args.chroms = args.chroms[order]
    chr2int = {}
    for i, chrom in enumerate(args.chroms):
        chr2int[chrom] = i
    # find counts
    print("\r%s\rLoading initial counts" % (' ' * 80), end='', file=sys.stderr)
    chr_indices = hic.fends['chr_indices'][...]
    for i, chrom in enumerate(hic.fends['chromosomes'][...]):
        if chrom not in args.chroms and "chr%s" % chrom not in args.chroms:
            hic.filter[chr_indices[i]:chr_indices[i + 1]] = 0
    cis = hic.data['cis_data'][...].astype(numpy.int64)
    cis[numpy.where((hic.filter[cis[:, 0]] == 0) | (hic.filter[cis[:, 1]] == 0))[0], 2] = 0
    trans = hic.data['trans_data'][...].astype(numpy.int64)
    trans[numpy.where((hic.filter[trans[:, 0]] == 0) | (hic.filter[trans[:, 1]] == 0))[0], 2] = 0
    if hic.corrections is not None:
        corrections = hic.corrections
    else:
        corrections = numpy.ones(hic.filter.shape[0], dtype=numpy.float32)
    corrections[numpy.where(hic.filter == 0)[0]] = 0
    fends = hic.fends['fends'][...]
    assembly = hic.data.attrs['fendfilename'].split('/')[-1].split('_')[0]

    res = args.resolution
    outfile = h5py.File(args.output, 'w')
    print("\r%s\rChilling res: %i - binning" % (' ' * 80, res), end='', file=sys.stderr)

    # convert indices into bin numbers
    sizes = numpy.zeros(len(args.chroms) + 1, numpy.int64)
    all_indices = numpy.zeros(len(args.chroms) + 1, dtype=numpy.int64)
    all_mapping = numpy.zeros(fends.shape[0], dtype=numpy.int64)
    for j, chrom in enumerate(args.chroms):
        if chrom in hic.chr2int:
            chrint = hic.chr2int[chrom]
        else:
            chrint = hic.chr2int[chrom.strip('chr')]
        all_indices[j + 1] = all_indices[j] + chr_indices[chrint + 1] - chr_indices[chrint]
        sizes[j + 1] = (length_array[j] - 1) / res + 1 + sizes[j]
        all_mapping[chr_indices[chrint]:chr_indices[chrint + 1]] = (
            fends['mid'][chr_indices[chrint]:chr_indices[chrint + 1]] / res + sizes[j])
    correction_sums = numpy.bincount(all_mapping, weights=corrections, minlength=sizes[-1])
    where = numpy.where(correction_sums != 0)[0]
    correction_sums[where] /= numpy.mean(correction_sums[where])
    cis[:, 0] = all_mapping[cis[:, 0]]
    cis[:, 1] = all_mapping[cis[:, 1]]
    trans[:, 0] = all_mapping[trans[:, 0]]
    trans[:, 1] = all_mapping[trans[:, 1]]
    where = numpy.where(cis[:, 0] > cis[:, 1])[0]
    temp = cis[where, 0]
    cis[where, 0] = cis[where, 1]
    cis[where, 1] = temp
    where = numpy.where(trans[:, 0] > trans[:, 1])[0]
    temp = trans[where, 0]
    trans[where, 0] = trans[where, 1]
    trans[where, 1] = temp

    # find unique bin pairs
    cwhere = numpy.where(cis[:, 2])[0]
    twhere = numpy.where(trans[:, 2])[0]
    cis_indices = cis[cwhere, 0] * sizes[-1] + cis[cwhere, 1]
    trans_indices = trans[twhere, 0] * sizes[-1] + trans[twhere, 1]
    new_indices = numpy.unique(numpy.r_[cis_indices, trans_indices])
    cis_pos = numpy.searchsorted(new_indices, cis_indices)
    trans_pos = numpy.searchsorted(new_indices, trans_indices)
    counts = numpy.zeros(new_indices.shape[0], dtype=numpy.int32)
    bin1 = numpy.zeros(new_indices.shape[0], dtype=numpy.int64)
    bin2 = numpy.zeros(new_indices.shape[0], dtype=numpy.int64)
    counts[:] = numpy.bincount(cis_pos, weights=cis[cwhere, 2], minlength=counts.shape[0])
    counts += numpy.bincount(trans_pos, weights=trans[twhere, 2], minlength=counts.shape[0]).astype(numpy.int32)

    bin1[:] = new_indices / sizes[-1]
    bin2[:] = new_indices % sizes[-1]
    where = numpy.where(bin1 == bin2)[0]
    counts[where] *= 2

    # find weights
    weights = numpy.copy(correction_sums)
    where = numpy.where(weights > 0)[0]
    weights[where] = 1.0 / weights[where]
    where = numpy.where(weights == 0)[0]
    weights[where] = numpy.nan        

    # find bin1 indices
    bin_indices = numpy.r_[0, numpy.cumsum(numpy.bincount(bin1, minlength=sizes[-1]))].astype(numpy.int64)

    # create bin arrays
    chrom_array = numpy.repeat(numpy.arange(len(args.chroms)), sizes[1:] - sizes[:-1])
    start_array = numpy.zeros(sizes[-1], dtype=numpy.int32)
    stop_array = numpy.zeros(sizes[-1], dtype=numpy.int32)
    for j, chrom in enumerate(args.chroms):
        start_array[sizes[j]:sizes[j + 1]] = numpy.arange(sizes[j + 1] - sizes[j]) * res
        stop_array[sizes[j]:sizes[j + 1]] = numpy.arange(1, sizes[j + 1] - sizes[j] + 1) * res
        stop_array[sizes[j + 1] - 1] = min(stop_array[sizes[j + 1] - 1], length_array[j])

    print("\r%s\rChilling res: %i - writing" % (' ' * 80, res), end='', file=sys.stderr)
    # write resolution attributes
    outfile.attrs['bin-size'] = res
    outfile.attrs['bin-type'] = 'fixed'
    outfile.attrs['creation-date'] = str(datetime.datetime.now())
    outfile.attrs['format'] = "HDF5::Cooler"
    outfile.attrs['format-url'] = "https://github.com/mirnylab/cooler"
    outfile.attrs['format-version'] = 2
    outfile.attrs['generated-by'] = 'hifive'
    outfile.attrs['genome-assembly'] = assembly
    outfile.attrs['metadata'] = '{}'
    outfile.attrs['nbins'] = sizes[-1].astype(numpy.int64)
    outfile.attrs['nchroms'] = len(args.chroms)
    outfile.attrs['nnz'] = counts.shape[0]
    outfile.attrs['id'] = 'null'

    # write bin group data
    bin_group = outfile.create_group('bins')
    dt = h5py.special_dtype(enum=('i', chr2int))
    temp = bin_group.create_dataset('chrom', (sizes[-1],), dtype=dt,
                                    compression='gzip', compression_opts=6)
    temp[:] = chrom_array
    bin_group.create_dataset(name='start', data=start_array,
                             compression='gzip', compression_opts=6)
    bin_group.create_dataset(name='end', data=stop_array,
                             compression='gzip', compression_opts=6)
    weight_dataset = bin_group.create_dataset(name='weight', data=weights,
                                              compression='gzip', compression_opts=6)

    # write chrom group data
    chroms_group = outfile.create_group('chroms')
    chroms_group.create_dataset(name='length', data=length_array,
                                compression='gzip', compression_opts=6)
    chroms_group.create_dataset('name', (len(args.chroms),), dtype='S32',
                                compression='gzip', compression_opts=6)
    for j, chrom in enumerate(args.chroms):
        chroms_group['name'][j] = chrom

    # write indexes group data
    indexes = outfile.create_group('indexes')
    indexes.create_dataset(name='bin1_offset', data=bin_indices,
                           compression='gzip', compression_opts=6)
    indexes.create_dataset(name='chrom_offset', data=sizes,
                           compression='gzip', compression_opts=6)

    # write pixels group data
    pixels = outfile.create_group('pixels')
    pixels.create_dataset(name='bin1_id', data=bin1,
                          compression='gzip', compression_opts=6)
    pixels.create_dataset(name='bin2_id', data=bin2,
                          compression='gzip', compression_opts=6)
    pixels.create_dataset(name='count', data=counts,
                          compression='gzip', compression_opts=6)

    outfile.close()
    print("\r%s\r" % (' ' * 80), end='', file=sys.stderr)

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Output cooler formated heatmap file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="input", type=str, action='store', help="HiFive project file name")
    parser.add_argument(dest="output", type=str, action='store', help="Output prefix")
    parser.add_argument('-r', dest="resolution", type=int, action='store', default=10000,
        help="Resolution for binning")
    parser.add_argument('-c', dest="chroms", type=str, action='store', default='',
        help="A comma-separated list of chromosomes to include. Leave blank to include all chromosomes")
    return parser


if __name__ == "__main__":
    main()
