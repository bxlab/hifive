#!/usr/bin/env python

import sys
import struct
import os
import argparse

import hifive
import numpy

def main():
    parser = generate_parser()
    args = parser.parse_args()
    hic = hifive.HiC(args.input)
    if args.chroms == '':
        args.chroms = hic.fends['chromosomes'][...]
    else:
        args.chroms = numpy.array(args.chroms.split(','))
    assembly = hic.data.attrs['fendfilename'].split('/')[-1].split('_')[0]
    header_size = (4 + 16 + 4 + 4 + len(assembly) + 1 + 4 + 4 + 4 + 4 + 4 + 4)
    intrachrom_start = header_size
    #version = "hifive1.0".ljust(16, '\0')
    version = "1.3".ljust(16, '\0')
    MCV = 0.0
    chrom_sizes = numpy.zeros(args.chroms.shape[0], dtype=numpy.int32)
    data = {}
    row_locations = {}
    for i, chrom in enumerate(args.chroms):
        chrint = hic.chr2int[chrom]
        chrom_sizes[i] = hic.fends['chrom_sizes'][chrint]
    order = numpy.argsort(chrom_sizes)[::-1]
    chrom_sizes = chrom_sizes[order]
    args.chroms = args.chroms[order]
    file_list = []
    for i, chrom in enumerate(args.chroms):
        N = (chrom_sizes[i] - 1) / args.res + 1
        temp = hic.cis_heatmap(chrom, start=0, stop=(N * args.res), binsize=args.res, datatype=args.datatype,
            arraytype='upper', includediagonal=True)
        if temp is None:
            continue
        header_size += len(chrom) + 4 + 4 + 8
        where = numpy.where(temp[:, 0])[0]
        temp[where, 0] /= temp[where, 1]
        temp = temp[where, 0]
        indices = numpy.triu_indices(N, 0)
        if args.export:
            output = open("%s_chr%s.raw" % (args.output, chrom), 'w')
            file_list.append(["chr%s" % chrom, "%s/%s_chr%s.raw" % (os.getcwd(), args.output, chrom)])
            for j in range(temp.shape[0]):
                print >> output, "%i\t%i\t%f" % (indices[0][where[j]] * args.res, indices[1][where[j]] * args.res, temp[j])
            output.close()
        data[chrom] = [indices[1][where], temp]
        row_locations[chrom] = numpy.r_[0, numpy.cumsum(numpy.bincount(indices[0][where], minlength=N) * 8)].astype(numpy.int64)
    interchrom_start = 0
    interchrom_keys = []
    interchrom_data = {}
    interchrom_rows = {}
    if args.trans:
        interchrom_start = header_size
        for i, chrom in enumerate(args.chroms):
            N = (chrom_sizes[i] - 1) / args.res + 1
            for j, chrom2 in enumerate(args.chroms[:i]):
                M = (chrom_sizes[j] - 1) / args.res + 1
                temp = hic.trans_heatmap(chrom, chrom2, start1=0, stop1=(N * args.res), start2=0, stop2=(M * args.res),
                    binsize=args.res, datatype=args.datatype)
                if temp is None:
                    continue
                key = "chr%s\tchr%s\0" % (chrom, chrom2)
                interchrom_keys.append(key)
                where = numpy.where(temp[:, :, 1])
                temp = temp[where[0], where[1], 0] / temp[where[0], where[1], 1]
                interchrom_data[key] = [where[1], temp]
                interchrom_rows[key] = numpy.r_[0, numpy.cumsum(numpy.bincount(where[0], minlength=M) * 8)].astype(numpy.int64)
                file_list.append(["chr%s\tchr%s" % (chrom2, chrom),
                                  "%s/%s_chr%s_by_chr%s.raw" % (os.getcwd(), args.output, chrom2, chrom)])
                if args.export:
                    output = open("%s_chr%s_by_chr%s.raw" % (args.output, chrom2, chrom), 'w')
                    for k in range(where[0].shape[0]):
                        print >> output, "%i\t%i\t%f" % (where[1][k] * args.res, where[0][k] * args.res, temp[k])
                    output.close()
                header_size += len(key) + 8
    if args.export:
        output = open("%s.sizes" % args.output, 'w')
        for i, chrom in enumerate(args.chroms):
            print >> output, "chr%s\t%i" % (chrom, chrom_sizes[i])
        output.close()
        output = open("%s.files" % args.output, 'w')
        for i in range(len(file_list)):
            print >> output, "%s\t%s" % (file_list[i][0], file_list[i][1])
        output.close()
    output = open("%s.butlr" % args.output, 'wb')
    # write header size
    output.write(struct.pack('I', header_size))
    # write code version number
    for i in range(len(version)):
        output.write(struct.pack('c', version[i]))
    # write location of chromosome size/intrachromosomal
    output.write(struct.pack('I', intrachrom_start))
    # write location of interchromosomal
    output.write(struct.pack('I', interchrom_start))
    # write assembly
    for i in range(len(assembly)):
        output.write(struct.pack('c', assembly[i]))
    output.write(struct.pack('c', '\0'))
    # write resolution
    output.write(struct.pack('I', args.res))
    # write most common value
    output.write(struct.pack('<f', MCV))
    # write empty field holders
    output.write(struct.pack('IIII', 0, 0, 0, 0))
    pos = header_size
    # write intra-chromosome information
    for i, chrom in enumerate(args.chroms):
        if chrom not in data:
            continue
        # write chrom name
        name = 'chr%s\0' % chrom
        for j in range(len(name)):
            output.write(struct.pack('c', name[j]))
        # write chrom size
        output.write(struct.pack('I', chrom_sizes[i]))
        # write row body location
        row_locations[chrom] += pos
        row_locations[chrom][numpy.where(row_locations[chrom][1:] == row_locations[chrom][:-1])[0]] = 0
        pos += data[chrom][0].shape[0] * 8
        output.write(struct.pack('Q', pos))
        # adjust for row data
        pos += row_locations[chrom].shape[0] * 8
    # write inter-chromosome information
    for key in interchrom_keys:
        # write chrom name
        for j in range(len(key)):
            output.write(struct.pack('c', key[j]))
        # write row body location
        interchrom_rows[key] += pos
        interchrom_rows[key][numpy.where(interchrom_rows[key][1:] == interchrom_rows[key][:-1])[0]] = 0
        pos += interchrom_data[key][0].shape[0] * 8
        output.write(struct.pack('Q', pos))
        # adjust for row data
        pos += interchrom_rows[key].shape[0] * 8
    # write intra-chromosome data
    for i, chrom in enumerate(args.chroms):
        if chrom not in data:
            continue
        # write chrom columns and values
        for j in range(data[chrom][0].shape[0]):
            output.write(struct.pack('I', data[chrom][0][j]))
            output.write(struct.pack('<f', float("%f" % data[chrom][1][j])))
        # write row locations
        for j in range(row_locations[chrom].shape[0]):
            output.write(struct.pack('Q', row_locations[chrom][j]))
    for key in interchrom_keys:
        # write chrom columns and values
        for j in range(interchrom_data[key][0].shape[0]):
            output.write(struct.pack('I', interchrom_data[key][0][j]))
            output.write(struct.pack('<f', float("%f" % interchrom_data[key][1][j])))
        # write row locations
        for j in range(interchrom_rows[key].shape[0]):
            output.write(struct.pack('Q', interchrom_rows[key][j]))
    output.close()


def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Output BUTLR formated heatmap file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(dest="input", type=str, action='store', help="HiFive project file name")
    parser.add_argument(dest="output", type=str, action='store', help="Output prefix")
    parser.add_argument('-r', dest="res", type=int, action='store', default=1000000,
        help="Resolution for binning")
    parser.add_argument('-c', dest="chroms", type=str, action='store', default='',
        help="A comma-separated list of chromosomes to include. Leave blank to include all chromosomes")
    parser.add_argument('-t', dest="trans", action='store_true', help="Include trans interaction heatmaps")
    parser.add_argument('-d', dest="datatype", action='store', default='raw',
        choices=['raw', 'fend', 'enrichment'], help="Type of data to return")
    parser.add_argument('-e', dest='export', action='store_true', help='export text files for native BUTLR conversion')
    return parser


if __name__ == "__main__":
    main()
