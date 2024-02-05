#!/usr/bin/env python
"""Description: Multi-resolution heatmap plotter and data retrieval.

Copyright (c) 2015 Michael Sauria <mike.sauria@jhu.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import sys
import struct
import math
import argparse as ap

try:
    from PIL import Image
except:
    pass
import numpy


def main():
    parser = generate_parser()
    args = parser.parse_args()
    MRH = MrhSquareDataProvider( source=args.heatmap )
    temp = MRH.get_data( chromosomes=True )
    if args.chrom2 is None:
        args.chrom2 = args.chrom
    if args.chrom.strip('chr') not in temp['chromosomes'] or args.chrom2.strip('chr') not in temp['chromosomes']:
        print("Chromosome(s) don't appear to be in this MRH file.", file=sys.stderr)
        return 1
    if args.chrom != args.chrom2 and not temp['includes_trans']:
        print("There is no trans data in this MRH file.", file=sys.stderr)
        return 1
    if not args.text and 'PIL' not in sys.modules.keys():
        print("The PIL package is needed to make plots.", file=sys.stderr)
        return 1
    header = MRH.get_data( header=True, chrom1=args.chrom, chrom2=args.chrom2 )
    if args.minres is None:
        args.minres = header['lres']
    if args.maxres is None:
        args.maxres = header['hres']
    if args.chrom == args.chrom2:
        if args.start is None:
            args.start = header['start']
        if args.start2 is None:
            args.start2 = header['start']
        if args.end is None:
            args.end = header['stop']
        if args.end2 is None:
            args.end2 = header['stop']
    elif MRH.transpose:
        if args.start is None:
            args.start = header['start2']
        if args.start2 is None:
            args.start2 = header['start1']
        if args.end is None:
            args.end = header['stop2']
        if args.end2 is None:
            args.end2 = header['stop1']
    else:
        if args.start is None:
            args.start = header['start1']
        if args.start2 is None:
            args.start2 = header['start2']
        if args.end is None:
            args.end = header['stop1']
        if args.end2 is None:
            args.end2 = header['stop2']
    MRH.get_data(
        chrom1=args.chrom,
        chrom2=args.chrom2,
        start1=args.start,
        start2=args.start2,
        stop1=args.end,
        stop2=args.end2,
        min_resolution=args.minres,
        max_resolution=args.maxres,
        out_fname=args.output,
        as_text=args.text
    )

# -----------------------------------------------------------------------------
class MrhSquareDataProvider( object ):
    """
    """
    settings = {
        'chrom1'         : 'str',
        'start1'         : 'long',
        'stop1'          : 'long',
        'chrom2'         : 'str',
        'start2'         : 'long',
        'stop2'          : 'long',
        'min_resolution' : 'long',
        'max_resolution' : 'long',
        'header'         : 'bool',
        'chromosomes'    : 'bool',
    }

    def __init__( self, source, **kwargs ):
        """
        :param chrom1: chromosome for sequence1
        :type chrom1: str.

        :param start1: starting bp in sequence 1
        :type start1: long

        :param stop1: stopping bp in sequence 1
        :type stop1: long

        :param chrom2: chromosome for sequence2
        :type chrom2: str.
        
        :param start2: starting bp in sequence 2
        :type start2: long

        :param stop2: stopping bp in sequence 2
        :type stop2: long

        :param min_resolution: largest bin size to return
        :type min_resolution: long

        :param max_resolution: smallest bin size to return
        :type max_resolution: long

        :param header: T/F send chromosome-specific header data instead
        :type resolution: bool

        :param chromosomes: T/F send chromosome list instead
        :type resolution: bool
        """
        self.source = open( source, 'rb' )

    def get_data( self, chrom1='', start1=0, stop1=0, chrom2='', start2=0, stop2=0, min_resolution=0,
                  max_resolution=0, header=False, chromosomes=False, out_fname='', as_text=False ):
        self.chrom1 = chrom1.strip( 'chr' )
        self.chrom2 = chrom2.strip( 'chr' )
        self.start1 = start1
        self.start2 = start2
        self.stop1 = stop1
        self.stop2 = stop2
        self.minres = min_resolution
        self.maxres = max_resolution
        self.header_only = header
        self.chroms_only = chromosomes
        if chrom1 == chrom2 and ( stop1 > start2 and stop2 > start1):
            self.overlap = True
        else:
            self.overlap = False
        self.window1 = self.stop1 - self.start1
        self.window2 = self.stop2 - self.start2
        self.header = None
        if self.chroms_only:
            return self.load_chrom_data( self.source )
        elif self.header_only:
            return self.load_header( self.source )
        elif self.chrom1 != '' and self.chrom2 != '':
            if self.header is None:
                self.load_header( self.source )
            data = []
            if self.transpose:
                self.start1, self.start2 = self.start2, self.start1
                self.stop1, self.stop2 = self.stop2, self.stop1
                self.window1, self.window2 = self.window2, self.window1
            if self.chrom1 != self.chrom2:
                for square in self.paint_trans_canvas( self.source ):
                    data.append( self.interpolate_square( square ) )
            else:
                for square in self.paint_cis_canvas( self.source ):
                    data.append( self.interpolate_square( square ) )
            if self.transpose:
                self.start1, self.start2 = self.start2, self.start1
                self.stop1, self.stop2 = self.stop2, self.stop1
                self.window1, self.window2 = self.window2, self.window1
            if as_text:
                self.write_data( data, out_fname )
            else:
                self.plot_data( data, out_fname )

    def load_chrom_data( self, infile ):
        """
        """
        self.header = {}
        # get and validate magic number for mrh type
        mrh_magic_number = '42054205'
        mrh_magic_number_size = 4
        infile.seek( 0 )
        magic_number = ''.join( c.encode( 'hex' ) for c in infile.read( mrh_magic_number_size ).replace( '\\x', '' ) )
        if magic_number != mrh_magic_number:
            raise TypeError( 'File does not appear to be a multi-resolution heatmap file' )

        # get the number of chromosomes and whether the file includes inter-chromosome maps
        int_float32_size = 4
        trans, num_chroms = struct.unpack( '>ii', infile.read( int_float32_size * 2 ) )
        self.trans = trans
        self.num_chroms = num_chroms
        name_sizes = struct.unpack( '>' + 'i' * num_chroms, infile.read( int_float32_size * num_chroms ) )

        # create dictionary of chromosome names and indices, retrieve chrom indices for requested data
        self.chr2int = {}
        for i in range( num_chroms ):
            self.chr2int[ ''.join( struct.unpack( '>' + 'c' * name_sizes[i],
                infile.read( name_sizes[i] ) ) ).strip( ) ] = i
        return { 'chromosomes':self.chr2int.keys( ), 'includes_trans':bool( trans ) }

    def load_header( self, infile ):
        int_float32_size = 4
        short_size = 2
        self.load_chrom_data( infile )
        trans = self.trans
        num_chroms = self.num_chroms
        if self.chrom1 != self.chrom2 and not trans:
            raise TypeError( 'File does not appear to contain inter-chromosome data' )
        if trans:
            num_pairings = ( num_chroms * ( num_chroms + 1 ) ) / 2
        else:
            num_pairings = num_chroms
        if self.chrom1 in self.chr2int and self.chrom2 in self.chr2int:
            chr_index1 = self.chr2int[ self.chrom1 ]
            chr_index2 = self.chr2int[ self.chrom2 ]
        else:
            raise KeyError( 'File does not appear to contain data for the requested chromosome(s)' )
        self.transpose = False
        if trans:
            if chr_index1 > chr_index2:
                self.transpose = True
                chr_index1, chr_index2 = chr_index2, chr_index1
            pairing_index = chr_index1 * ( num_chroms - 1 ) - ( chr_index1 * ( chr_index1 - 1 ) ) / 2 + chr_index2
        else:
            pairing_index = chr_index1

        # get index of chrom/chrom-pair data starting byte and ending byte
        infile.seek( pairing_index * int_float32_size, 1 )
        start_index, stop_index = struct.unpack( '>' + 'i' * 2, infile.read( int_float32_size * 2 ) )
        self.header[ 'offset' ] = start_index
        self.header[ 'total_bytes' ] = stop_index - start_index
        skip = num_pairings - pairing_index - 1

        # get number of top-layer partitions
        if chr_index1 != chr_index2:
            # skip intra-chromosomal partitions and find inter-chromosome partitions for requested chroms
            infile.seek( ( skip + chr_index1 + num_chroms ) * int_float32_size, 1 )
            self.header[ 'n' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            infile.seek( ( chr_index2 - chr_index1 - 1 ) * int_float32_size, 1 )
            self.header[ 'm' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            skip = num_chroms - chr_index2 - 1
        else:
            # find intra-chromosome partitions for requested chrom
            infile.seek( ( skip + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'n' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            if trans:
                # skip inter-chromosomal paritions
                skip = num_chroms * 2 - chr_index1 - 1
            else:
                skip = num_chroms - chr_index1 - 1

        # get total number of data bins for chrom/chrom-pair
        infile.seek( ( skip + pairing_index ) * int_float32_size, 1 )
        self.header[ 'data_bins' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
        skip = num_pairings - pairing_index - 1

        # get total number of index bins for chrom/chrom-pair
        infile.seek( ( skip + pairing_index ) * int_float32_size, 1 )
        self.header[ 'index_bins' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
        self.header[ 'shape_bins' ] = ( self.header[ 'total_bytes' ] - ( self.header[ 'data_bins' ] +
                                        self.header[ 'index_bins' ] ) * int_float32_size ) / short_size
        self.header[ 'total_bins' ] = ( self.header[ 'index_bins' ] + self.header[ 'data_bins' ] +
                                        self.header[ 'shape_bins' ] )
        skip = num_pairings - pairing_index - 1

        # get starting coordinate(s)
        if chr_index1 != chr_index2:
            # skip intra-chromosomal start coordinates and find inter-chromosome start coordinates for requested chroms
            infile.seek( ( skip + num_chroms + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'start1' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            infile.seek( ( chr_index2 - chr_index1 - 1 ) * int_float32_size, 1 )
            self.header[ 'start2' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            skip = num_chroms - chr_index2 - 1
        else:
            # find intra-chromosome start coordinate for requested chrom
            infile.seek( ( skip + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'start' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            if trans:
                # skip inter-chromosomal start coordinates
                skip = num_chroms * 2 - chr_index1 - 1
            else:
                skip = num_chroms - chr_index1 - 1

        # get stopping coordinate(s)
        if chr_index1 != chr_index2:
            # skip intra-chromosomal stop coordinates and find inter-chromosome stop coordinates for requested chroms
            infile.seek( ( skip + num_chroms + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'stop1' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            infile.seek( ( chr_index2 - chr_index1 - 1 ) * int_float32_size, 1 )
            self.header[ 'stop2' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            skip = num_chroms - chr_index2 - 1
        else:
            # find intra-chromosome stop coordinate for requested chrom
            infile.seek( ( skip + chr_index1 ) * int_float32_size, 1 )
            self.header[ 'stop' ] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]
            if trans:
                # skip inter-chromosomal stop coordinates
                skip = num_chroms * 2 - chr_index1 - 1
            else:
                skip = num_chroms - chr_index1 - 1

        # get minimum score
        infile.seek( ( skip + pairing_index ) * int_float32_size, 1 )
        self.header[ 'minscore' ] = struct.unpack( '>f', infile.read( int_float32_size ) )[0]
        skip = num_pairings - pairing_index - 1

        # get maximum score
        infile.seek( ( skip + pairing_index ) * int_float32_size, 1 )
        self.header[ 'maxscore' ] = struct.unpack( '>f', infile.read( int_float32_size ) )[0]
        skip = num_pairings - pairing_index - 1

        # get maximum bin size (lowest resolution)
        if chr_index1 != chr_index2:
            # skip intra-chromosomal largest bin size
            skip += 1
        infile.seek( skip * int_float32_size, 1 )
        self.header[ 'lres'] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]

        # get minimum bin size (highest resolution)
        if trans:
            # skip intra- or interchromosomal value (whichever is unneeded)
            infile.seek( int_float32_size, 1 )
        self.header[ 'hres'] = struct.unpack( '>i', infile.read( int_float32_size ) )[0]

        return self.header

    def paint_cis_canvas( self, infile ):
        """
        """
        if ( min( self.stop1, self.stop2 ) <= self.header[ 'start' ] or
             max( self.start1, self.start2 ) >= self.header[ 'stop' ] ):
            return []
        # ensure positions are within bounds of data present in file
        start_pos1 = max( 0, ( self.start1 - self.header[ 'start' ] ) / self.header[ 'lres' ] )
        end_pos1 = min( self.header[ 'n' ], ( self.stop1 - self.header[ 'start' ] ) / self.header[ 'lres' ] + 1 )
        start_pos2 = max( 0, ( self.start2 - self.header[ 'start' ] ) / self.header[ 'lres' ] )
        end_pos2 = min( self.header[ 'n' ], ( self.stop2 - self.header[ 'start' ] ) / self.header[ 'lres' ] + 1 )
        
        if self.overlap:
            # if data overlap, break query into two parts, since data only covers upper-triangle
            outdata = []
            mid_coord = min(self.stop1, self.stop2)
            mid_start = min( self.header[ 'n' ], ( mid_coord - self.header[ 'start' ] ) / self.header[ 'lres' ] )
            mid_stop = min( self.header[ 'n' ], ( mid_coord - self.header[ 'start' ] ) / self.header[ 'lres' ] + 1 )
            if self.start1 <= self.start2:
                self.reverse = False
                self.eff_start1, self.eff_start2 = self.start1, self.start2
                self.eff_stop1, self.eff_stop2 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( infile, start_pos1, mid_stop, start_pos2, mid_stop )
            else:
                self.reverse = True
                self.eff_start2, self.eff_start1 = self.start1, self.start2
                self.eff_stop2, self.eff_stop1 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( infile, start_pos2, mid_stop, start_pos1, mid_stop )
            if self.stop1 > mid_coord:
                self.reverse = True
                self.eff_start2, self.eff_start1 = self.start1, self.start2
                self.eff_stop2, self.eff_stop1 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( infile, start_pos2, end_pos2, mid_start, end_pos1 )
            elif self.stop2 > mid_coord:
                self.reverse = False
                self.eff_start1, self.eff_start2 = self.start1, self.start2
                self.eff_stop1, self.eff_stop2 = self.stop1, self.stop2
                outdata += self.paint_cis_upper_level( infile, start_pos1, end_pos1, mid_start, end_pos2 )

        # if data don't overlap, figure out which range is upstream
        elif self.start2 < self.start1:
            self.reverse = True
            self.eff_start2, self.eff_start1 = self.start1, self.start2
            self.eff_stop2, self.eff_stop1 = self.stop1, self.stop2
            outdata = self.paint_cis_upper_level( infile, start_pos2, end_pos2, start_pos1, end_pos1 )
        else:
            self.reverse = False
            self.eff_start1, self.eff_start2 = self.start1, self.start2
            self.eff_stop1, self.eff_stop2 = self.stop1, self.stop2
            outdata = self.paint_cis_upper_level( infile, start_pos1, end_pos1, start_pos2, end_pos2 )
        return outdata

    def paint_cis_upper_level( self, infile, start_pos1, end_pos1, start_pos2, end_pos2 ):
        """
        """
        outdata = []
        resolution = self.header[ 'lres' ]
        for i in range( start_pos1, end_pos1 ):
            pos_index = i * ( self.header[ 'n' ] - 1 ) - ( i * ( i - 1 ) ) / 2 + max( i, start_pos2 )
            # Find position in file for data with 'i' as upstream interaction
            infile.seek( self.header[ 'offset' ] + pos_index * 4 )
            data = struct.unpack( '>' + 'f' * ( end_pos2 - max( start_pos2, i ) ),
                   infile.read( ( end_pos2 - max( start_pos2, i ) ) * 4 ) )
            # Find position in file for indices with 'i' as upstream interaction
            infile.seek( self.header[ 'offset' ] + ( pos_index + self.header[ 'data_bins' ] ) * 4 )
            indices = struct.unpack( '>' + 'i' * ( end_pos2 - max( start_pos2, i ) ),
                      infile.read( ( end_pos2 - max( start_pos2, i ) ) * 4 ) )
            # Find position in file for shapes with 'i' as upstream interaction
            infile.seek( self.header[ 'offset' ] + ( self.header[ 'data_bins' ] +
                         self.header[ 'index_bins' ] ) * 4 + pos_index * 2 )
            shapes = struct.unpack( '>' + 'h' * ( end_pos2 - max( start_pos2, i ) ),
                      infile.read( ( end_pos2 - max( start_pos2, i ) ) * 2 ) )
            for j in range( max( start_pos2, i ), end_pos2 ):
                k = j - max( start_pos2, i )
                # if data is valid, seek down to lower levels for smaller bin sizes
                if not math.isnan( data[k] ):
                    start1 = i * resolution + self.header[ 'start' ]
                    start2 = j * resolution + self.header[ 'start' ]
                    # if valid lower-level data exists, try retrieving it
                    if indices[k] != -1:
                        valid, new_outdata = self.paint_cis_lower_level( infile, indices[k], shapes[k],
                                                               resolution / 2, start1, start2 )
                    else:
                        valid = 0
                    if start1 == start2:
                        # if data is on the diagonal and not completely covered by lower-level data,
                        # only add one square
                        if resolution <= self.minres and valid < 3:
                            outdata.append( [start1, start2, resolution, data[k]] )
                    else:
                        # if data is off-diagonal and not completely covered by lower-level data, add it to dataset
                        if resolution <= self.minres and valid < 4:
                            outdata.append( [start1, start2, resolution, data[k]] )
                            # if data is off-diagonal, check if the mirror image is also needed
                            if self.overlap and start2 + resolution > self.eff_start1 and start1 < self.eff_stop2:
                                outdata.append( [start2, start1, resolution, data[k]] )
                    if valid > 0:
                        outdata += new_outdata
        # if sequence 2 was upstream of sequence 1, flip the x and y coordinates
        if self.reverse:
            for i in range( len( outdata ) ):
                outdata[i] = [ outdata[i][1], outdata[i][0] ] + outdata[i][2:]
        return outdata

    def paint_cis_lower_level( self, infile, index, shape, resolution, start1, start2 ):
        """
        """
        # don't return data the is higher resolution than requested
        if resolution < self.maxres:
            return 0, []
        outdata = []
        valid = 0
        shape = bin(shape)[2:].rjust( 4, '0' )
        num_data = shape.count('1')
        infile.seek( self.header[ 'offset' ] + index * 4 )
        data = struct.unpack( '>' + 'f' * num_data, infile.read( 4 * num_data ) )
        if index < self.header[ 'index_bins' ]:
            infile.seek( self.header[ 'offset' ] + ( index + self.header[ 'data_bins' ] ) * 4 )
            indices = struct.unpack( '>' + 'i' * num_data, infile.read( num_data * 4 ) )
            infile.seek( self.header[ 'offset' ] + ( self.header[ 'data_bins' ] + self.header[ 'index_bins' ] ) * 4 +
                         index * 2 )
            shapes = struct.unpack( '>' + 'h' * num_data, infile.read( num_data * 2 ) )
        else:
            indices = None
            shapes = None
        # for each bin, find if it has valid data, and if lower levels need to be queried
        pos = 0
        for i in range( 2 ):
            for j in range( 2 ):
                if shape[ -1 - j - i * 2 ] == '0':
                    continue
                if not math.isnan( data[pos] ):
                    start1b = start1 + i * resolution
                    start2b = start2 + j * resolution
                    # if the bin is valid but out of the query range, do return data
                    if ( start1b > self.eff_stop1 or start1b + resolution < self.eff_start1 or
                         start2b > self.eff_stop2 or start2b + resolution < self.eff_start2 ):
                        valid += 1
                    else:
                        # if the bin is valid and has valid data at a lower level, go down a level
                        if indices is not None and indices[pos] != -1:
                            new_valid, new_outdata = self.paint_cis_lower_level( infile, indices[pos], shapes[pos],
                                                     resolution / 2, start1b, start2b )
                        else:
                            new_valid = 0
                        # if data is high enough resolution and isn't completely covered by
                        # lower-level bins, return bin data
                        if resolution <= self.minres:
                            if start1b == start2b:
                                if new_valid < 3:
                                    outdata.append( [start1b, start2b, resolution, data[pos]] )
                            else:
                                if new_valid < 4:
                                    outdata.append( [start1b, start2b, resolution, data[pos]] )
                                if ( self.overlap and start2b + resolution > self.eff_start1 and
                                     start1b < self.eff_stop2 ):
                                    outdata.append( [start2b, start1b, resolution, data[pos]] )
                        if new_valid > 0:
                            outdata += new_outdata
                        valid += 1
                pos += 1
        return valid, outdata

    def paint_trans_canvas( self, infile ):
        """
        """
        outdata = []
        start_pos1 = max( 0, ( self.start1 - self.header[ 'start1' ] ) / self.header[ 'lres' ] )
        end_pos1 = min( self.header[ 'n' ], ( self.stop1 - self.header[ 'start1' ] ) / self.header[ 'lres' ] + 1 )
        start_pos2 = max( 0, ( self.start2 - self.header[ 'start2' ] ) / self.header[ 'lres' ] )
        end_pos2 = min( self.header[ 'm' ], ( self.stop2 - self.header[ 'start2' ] ) / self.header[ 'lres' ] + 1 )
        resolution = self.header[ 'lres' ]

        for i in range( start_pos1, end_pos1 ):
            pos_index = i * self.header[ 'm' ] + start_pos2
            # Find position in file for data with 'i' as upstream interaction
            infile.seek( self.header[ 'offset' ] + pos_index * 4 )
            data = struct.unpack( '>' + 'f' * ( end_pos2 - start_pos2 ), infile.read( ( end_pos2 - start_pos2 ) * 4 ) )
            # Find position in file for indices with 'i' as upstream interaction
            infile.seek( self.header[ 'offset' ] + ( pos_index + self.header[ 'data_bins' ] ) * 4 )
            indices = struct.unpack( '>' + 'i' * ( end_pos2 - start_pos2 ),
                infile.read( ( end_pos2 - start_pos2 ) * 4 ) )
            # Find position in file for shapes with 'i' as upstream interaction
            infile.seek( self.header[ 'offset' ] + ( self.header[ 'data_bins' ] + self.header[ 'index_bins' ] ) * 4 +
                         pos_index * 2 )
            shapes = struct.unpack( '>' + 'h' * ( end_pos2 - start_pos2 ),
                infile.read( ( end_pos2 - start_pos2 ) * 2 ) )
            for j in range( start_pos2, end_pos2 ):
                k = j - start_pos2
                # if data is valid, seek down to lower levels for smaller bin sizes
                if not math.isnan( data[k] ):
                    start1 = i * resolution + self.header[ 'start1' ]
                    start2 = j * resolution + self.header[ 'start2' ]
                    # if valid lower-level data exists, try retrieving it
                    if indices[k] != -1:
                        new_valid, new_outdata = self.paint_trans_lower_level( infile, indices[k], shapes[k],
                                                 resolution / 2, start1, start2 )
                    else:
                        new_valid = 0
                    # if top-level square not completely covered by lower-level squares, add it to dataset
                    if resolution <= self.minres and new_valid < 4:
                        outdata.append( [start1, start2, resolution, data[k]] )
                    if new_valid > 0:
                        outdata += new_outdata
        return outdata

    def paint_trans_lower_level( self, infile, index, shape, resolution, start1, start2 ):
        """
        """
        # don't return data the is higher resolution than requested
        if resolution < self.maxres:
            return 0, []
        outdata = []
        valid = 0
        shape = bin(shape)[2:].rjust( 4, '0' )
        num_data = shape.count('1')
        infile.seek( self.header[ 'offset' ] + index * 4 )
        data = struct.unpack( '>' + 'f' * num_data, infile.read( 4 * num_data ) )
        if index < self.header[ 'index_bins' ]:
            infile.seek( self.header[ 'offset' ] + ( index + self.header[ 'data_bins' ] ) * 4 )
            indices = struct.unpack( '>' + 'i' * num_data, infile.read( num_data * 4 ) )
            infile.seek( self.header[ 'offset' ] + ( self.header[ 'data_bins' ] + self.header[ 'index_bins' ] ) * 4 +
                         index * 2 )
            shapes = struct.unpack( '>' + 'h' * num_data, infile.read( num_data * 2 ) )
        else:
            indices = None
            shapes = None
        # for each bin, find if it has valid data, and if lower levels need to be queried
        pos = 0
        for i in range( 2 ):
            for j in range( 2 ):
                if shape[ -1 - ( i * 2 + j ) ] == '0':
                    continue
                if not math.isnan( data[pos] ):
                    start1b = start1 + i * resolution
                    start2b = start2 + j * resolution
                    # if the bin is valid but out of the query range, do return data
                    if ( start1b > self.stop1 or start1b + resolution < self.start1 or
                         start2b > self.stop2 or start2b + resolution < self.start2 ):
                        valid += 1
                    else:
                        if indices is not None and indices[pos] != -1:
                            new_valid, new_outdata = self.paint_trans_lower_level( infile, indices[pos], shapes[pos],
                                                     resolution / 2, start1b, start2b )
                        else:
                            new_valid = 0
                        # if data is high enough resolution and isn't completely covered by
                        # lower-level bins, return bin data
                        if resolution <= self.minres and new_valid < 4:
                            outdata.append( [start1b, start2b, resolution, data[pos]] )
                        if new_valid > 0:
                            outdata += new_outdata
                        valid += 1
                pos += 1
        return valid, outdata

    def interpolate_square( self, square ):
        """
        """
        if self.start1 == self.stop1 or self.start2 == self.stop2:
            return
        square_dict = {
            'x1' : square[0],
            'y1' : square[1],
            'x2' : square[0] + square[2],
            'y2' : square[1] + square[2],
            'color' : ( ( square[3] - self.header[ 'minscore' ] ) /
                        ( self.header[ 'maxscore' ] - self.header[ 'minscore' ] ) * 2.0 - 1.0 ),
            'value' : square[3]
        }
        if self.transpose:
            square_dict.update({
                'x1' : square_dict['y1'], 'y1' : square_dict['x1'],
                'x2' : square_dict['y2'], 'y2' : square_dict['x2'], 'value' : square[3]
            })
        return square_dict

    def plot_data( self, data, out_fname ):
        width = int( math.ceil( float(self.stop1 - self.start1) / self.maxres ) )
        height = int( math.ceil( float(self.stop2 - self.start2) / self.maxres ) )
        canvas = numpy.empty( ( height, width ), dtype=numpy.uint32 )
        canvas.shape = ( height, width )
        canvas.fill( int( 'ff000000', 16 ) )
        for square in data:
            value = int( round( square[ 'color' ] * 255 ) )
            if value >= 0.0:
                color = int( 'ff%02x%02xff' % ( 255 - value, 255 - value ), 16 )
            else:
                color = int( 'ffff%02x%02x' % ( 255 + value, 255 + value ), 16 )
            x1 = max( 0, min( width, int( math.floor( float( square[ 'x1' ] - self.start1 ) / self.maxres ) ) ) )
            x2 = max( 0, min( width, int( math.ceil( float( square[ 'x2' ] - self.start1 ) / self.maxres ) ) ) )
            y1 = max( 0, min( height, int( math.floor( float( square[ 'y1' ] - self.start2 ) / self.maxres ) ) ) )
            y2 = max( 0, min( height, int( math.ceil( float( square[ 'y2' ] - self.start2 ) / self.maxres ) ) ) )
            canvas[y1:y2, x1:x2] = color
        pilImage = Image.frombuffer( 'RGBA', ( width, height ), canvas, 'raw', 'RGBA', 0, 1 )
        pilImage.save( out_fname )

    def write_data( self, data, out_fname ):
        output = open( out_fname, 'w' )
        print("chr1\tstart1\tend1\tchr2\tstart2\tend2\tscore", file=output)
        for square in data:
            print("%s\t%i\t%i\t%s\t%i\t%i\t%f" % (self.chrom1, square['x1'], square['x2'],
                    self.chrom2, square['y1'], square['y2'], square['value']), file=output)
        output.close()


def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Produce a multi-resolution heatmap image or interval file from a HiC multi-resolution heatmap file."
    parser = ap.ArgumentParser(description=description)
    parser.add_argument("-c", "--chrom", dest="chrom", required=True, type=str, default=None,
        action='store', help="The first region chromosome.")
    parser.add_argument("-C", "--chrom2", dest="chrom2", required=False, type=str, default=None,
        action='store', help="The first region chromosome. If no value is passed, this will be set to the same value as 'chrom'.")
    parser.add_argument("-s", "--start", dest="start", required=False, type=int, default=None,
        action='store', help="The first region start coordinate to plot. If no value is passed, this will be set to the first bin position in the heatmap. [default: %(default)s]")
    parser.add_argument("-S", "--start2", dest="start2", required=False, type=int, default=None,
        action='store', help="The second region start coordinate to plot. If no value is passed, this will be set to the same value as start. [default: %(default)s]")
    parser.add_argument("-e", "--end", dest="end", required=False, type=int, default=None,
        action='store', help="The first region stop coordinate to plot. If no value is passed, this will be set to the last bin position in the heatmap [default: %(default)s]")
    parser.add_argument("-E", "--end2", dest="end2", required=False, type=int, default=None,
        action='store', help="The second region stop coordinate to plot. If no value is passed, this will be set to the last bin position in the heatmap [default: %(default)s]")
    parser.add_argument("-R", "--max-resolution", dest="maxres", required=False, type=int, default=None,
        action='store', help="The maximum resolution bound for returned data. If no value is passed, this will be set to the highest resolution available in the heatmap for the chromosome(s). [default: %(default)s]")
    parser.add_argument("-r", "--min-resolution", dest="minres", required=False, type=int, default=None,
        action='store', help="The minimum resolution bound for returned data. If no value is passed, this will be set to the lowest resolution available in the heatmap for the chromosome(s). [default: %(default)s]")
    parser.add_argument("-t", "--text", dest="text", required=False, default=False,
        action='store_true', help="Write a genomic interval text file instead of an image.")
    parser.add_argument(dest="heatmap", type=str,
        help="The name of a HiFive multi-resolution heatmap file to construct the image from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write the multi-resolution HiFive heatmap image to.")
    return parser

if __name__ == "__main__":
    main()
