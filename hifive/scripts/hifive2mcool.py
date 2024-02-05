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
    chromosomes = hic.fends['chromosomes'][...]
    chroms = []
    if args.chroms == '':
        for i in range(1, 23):
            if numpy.where(chromosomes == str(i))[0].shape[0] > 0:
                chroms.append(str(i))
        if numpy.where(chromosomes == 'X')[0].shape[0] > 0:
            chroms.append('X')
    else:
        for chrom in args.chroms:
            if numpy.where(chromosomes == chrom.strip('chr'))[0].shape[0] > 0:
                chroms.append(chrom.strip('chr'))
    args.chroms = chroms
    chr2int = {}
    for i, chrom in enumerate(args.chroms):
        chr2int['chr%s' % chrom] = i
    # find counts
    print("\r%s\rLoading initial counts" % (' ' * 80), end='', file=sys.stderr)
    chr_indices = hic.fends['chr_indices'][...]
    for i, chrom in enumerate(hic.fends['chromosomes'][...]):
        if chrom not in args.chroms:
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
    length_array = numpy.zeros(len(args.chroms), dtype=numpy.int32)
    temp = hic.fends['chrom_sizes'][...]
    for i, chrom in enumerate(args.chroms):
        length_array[i] = temp[hic.chr2int[chrom]]

    # determine resolution zoom levels
    nres = 1
    total_length = numpy.sum(length_array)
    while 256 * args.maxres * 2 ** (nres - 1) < total_length:
        nres += 1
    minres = args.maxres * 2 ** (nres - 1)
    resolutions = (args.maxres * 2 ** numpy.arange(nres))
    outfile = h5py.File(args.output, 'w')
    outfile.attrs['max-zoom'] = nres - 1
    min_vals = numpy.zeros(resolutions.shape[0], dtype=numpy.float64)
    max_vals = numpy.zeros(resolutions.shape[0], dtype=numpy.float64)
    mean_vals = numpy.zeros(resolutions.shape[0], dtype=numpy.float64)
    sd_vals = numpy.zeros(resolutions.shape[0], dtype=numpy.float64)

    # cycle through resolutions
    for i, res in enumerate(resolutions):
        print("\r%s\rChilling res: %i - binning" % (' ' * 80, res), end='', file=sys.stderr)
        outfile.attrs[str(nres - i - 1)] = res

        # convert indices into bin numbers
        sizes = numpy.zeros(len(args.chroms) + 1, numpy.int32)
        if i == 0:
            all_indices = numpy.zeros(len(args.chroms) + 1, dtype=numpy.int32)
            all_mapping = numpy.zeros(fends.shape[0], dtype=numpy.int64)
            for j, chrom in enumerate(args.chroms):
                chrint = hic.chr2int[chrom]
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

        else:            
            all_mapping = numpy.zeros(all_indices[-1], dtype=numpy.int64)
            for j, chrom in enumerate(args.chroms):
                chrint = hic.chr2int[chrom]
                sizes[j + 1] = (length_array[j] - 1) / res + 1 + sizes[j]
                all_mapping[all_indices[j]:
                            all_indices[j + 1]] = numpy.arange(all_indices[j + 1] - all_indices[j]) / 2 + sizes[j]
            correction_sums = numpy.bincount(all_mapping, weights=correction_sums, minlength=sizes[-1])
            correction_sums /= 4.0

            # find unique bin pairs
            indices = all_mapping[prev_bin1] * sizes[-1] + all_mapping[prev_bin2]
            new_indices = numpy.unique(indices)
            pos = numpy.searchsorted(new_indices, indices)
            counts = numpy.zeros(new_indices.shape[0], dtype=numpy.int32)
            bin1 = numpy.zeros(new_indices.shape[0], dtype=numpy.int64)
            bin2 = numpy.zeros(new_indices.shape[0], dtype=numpy.int64)
            counts[:] = numpy.bincount(pos, weights=prev_counts, minlength=counts.shape[0])
        bin1[:] = new_indices / sizes[-1]
        bin2[:] = new_indices % sizes[-1]
        where = numpy.where(bin1 == bin2)[0]
        counts[where] *= 2
        if i > 0:
            where = numpy.where(prev_bin1 == prev_bin2)
            pos = numpy.searchsorted(new_indices, all_mapping[prev_bin1[where]] * (sizes[-1] + 1))
            counts -= numpy.bincount(pos, weights=prev_counts[where], minlength=counts.shape[0]).astype(numpy.int32)

        # find weights
        weights = numpy.copy(correction_sums)
        where = numpy.where(weights > 0)[0]
        weights[where] = 1.0 / weights[where]
        where = numpy.where(weights == 0)[0]
        weights[where] = numpy.nan        

        # find bin1 indices
        bin_indices = numpy.r_[0, numpy.cumsum(numpy.bincount(bin1, minlength=sizes[-1]))].astype(numpy.int32)

        # create bin arrays
        chrom_array = numpy.repeat(numpy.arange(len(args.chroms)), sizes[1:] - sizes[:-1])
        start_array = numpy.zeros(sizes[-1], dtype=numpy.int32)
        stop_array = numpy.zeros(sizes[-1], dtype=numpy.int32)
        for j, chrom in enumerate(args.chroms):
            start_array[sizes[j]:sizes[j + 1]] = numpy.arange(sizes[j + 1] - sizes[j]) * res
            stop_array[sizes[j]:sizes[j + 1]] = numpy.arange(1, sizes[j + 1] - sizes[j] + 1) * res
            stop_array[sizes[j + 1] - 1] = min(stop_array[sizes[j + 1] - 1], length_array[j])

        # find value ranges for full and cis data
        mapping = numpy.zeros(sizes[-1], dtype=numpy.int32)
        for j in range(len(args.chroms)):
            mapping[sizes[j]:sizes[j + 1]] = j
        temp1 = mapping[bin1]
        temp2 = mapping[bin2]
        where = numpy.where(temp1 == temp2)[0]
        temp3 = counts[where] * weights[bin1[where]] * weights[bin2[where]]
        mean_vals[i] = numpy.mean(temp3)
        min_vals[i] = numpy.amin(temp3)
        max_vals[i] = numpy.amax(temp3)
        sd_vals[i] = numpy.std(temp3)

        print("\r%s\rChilling res: %i - writing" % (' ' * 80, res), end='', file=sys.stderr)
        # write resolution attributes
        dataset = outfile.create_group(str(nres - i - 1))
        dataset.attrs['bin-size'] = res
        dataset.attrs['bin-type'] = 'fixed'
        dataset.attrs['creation-date'] = str(datetime.datetime.now())
        dataset.attrs['format'] = "HDF5::Cooler"
        dataset.attrs['format-url'] = "https://github.com/mirnylab/cooler"
        dataset.attrs['format-version'] = 2
        dataset.attrs['generated-by'] = 'hifive'
        dataset.attrs['genome-assembly'] = assembly
        dataset.attrs['metadata'] = '{}'
        dataset.attrs['nbins'] = sizes[-1].astype(numpy.int64)
        dataset.attrs['nchroms'] = len(args.chroms)
        dataset.attrs['nnz'] = counts.shape[0]
        dataset.attrs['sum'] = numpy.sum(counts)
        dataset.attrs['id'] = 'null'
        dataset.attrs['min_value'] = min_vals[i]
        dataset.attrs['max_value'] = max_vals[i]
        dataset.attrs['mean_value'] = mean_vals[i]
        dataset.attrs['sd_value'] = sd_vals[i]

        # write bin group data
        bin_group = dataset.create_group('bins')
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
        chroms_group = dataset.create_group('chroms')
        chroms_group.create_dataset(name='length', data=length_array,
                                    compression='gzip', compression_opts=6)
        chroms_group.create_dataset('name', (len(args.chroms),), dtype='S32',
                                    compression='gzip', compression_opts=6)
        for j, chrom in enumerate(args.chroms):
            chroms_group['name'][j] = 'chr%s' % chrom

        # write indexes group data
        indexes = dataset.create_group('indexes')
        indexes.create_dataset(name='bin1_offset', data=bin_indices,
                               compression='gzip', compression_opts=6)
        indexes.create_dataset(name='chrom_offset', data=sizes,
                               compression='gzip', compression_opts=6)

        # write pixels group data
        pixels = dataset.create_group('pixels')
        pixels.create_dataset(name='bin1_id', data=bin1,
                              compression='gzip', compression_opts=6)
        pixels.create_dataset(name='bin2_id', data=bin2,
                              compression='gzip', compression_opts=6)
        pixels.create_dataset(name='count', data=counts,
                              compression='gzip', compression_opts=6)

        # set data for next round
        all_indices = sizes
        prev_bin1 = bin1
        prev_bin2 = bin2
        prev_counts = counts

    outfile.close()
    print("\r%s\r" % (' ' * 80), end='', file=sys.stderr)

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Output multi-resolution cooler formated heatmap file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="input", type=str, action='store', help="HiFive project file name")
    parser.add_argument(dest="output", type=str, action='store', help="Output prefix")
    parser.add_argument('-r', dest="maxres", type=int, action='store', default=10000,
        help="Maximum resolution for binning")
    parser.add_argument('-c', dest="chroms", type=str, action='store', default='',
        help="A comma-separated list of chromosomes to include. Leave blank to include all chromosomes")
    return parser


if __name__ == "__main__":
    main()
