#!/usr/bin/env python

"""
This is a module contains scripts for generating plots from compact and full matrices of interaction data.

Concepts
----------

These functions take either compact, upper-triangle, or full 3d data matrices.

Data can either be arranged in compact, upper-triangle, or complete (rectangular) arrays. With HiC data, compact arrays are N x M x 2, where N is the number of fends or bins, and M is the maximum distance between fends or bins. This is useful for working with sets of short interactions. When using 5C data, the compact format is an N x M x 2 array where N is the number of forward primers and M is the number of reverse primers. Data can be raw, fend-corrected, distance-dependence removed, or enrichment values. Arrays are 3-dimensional with observed values in the first layer of d3, expected values in the second layer of d3.

API Documentation
------------------
"""

import sys
from math import exp

import h5py
import numpy
try:
    from PIL import Image
except:
    pass
try:
    from pyx import canvas, path, style, bitmap, text, trafo
except:
    pass


def plot_compact_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                       min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                       returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a HiC compact array.

    :param data: A three-dimensional compact array of HiC interaction data.
    :type data: numpy array
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :param diagonal_included: If true, adjust output shape as necessary.
    :type diagonal_included: bool.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting compact array..."),
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    dim = data.shape[0]
    scaled = numpy.copy(data).astype(numpy.float64)
    where = numpy.where(data[:, :, 1] > 0)
    scaled[where[0], where[1], 0] /= scaled[where[0], where[1], 1]
    if logged:
        where = numpy.where(scaled[:, :, 0] <= 0)
        scaled[where[0], where[1], 1] = 0
        where = numpy.where(scaled[:, :, 1] > 0)
        scaled[where[0], where[1], 0] = numpy.log2(scaled[where[0], where[1], 0])
    scaled[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(scaled[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(scaled[where[0], where[1], 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        scaled[where[0], where[1], 0] /= maxscore
    else:
        scaled[where[0], where[1], 0] -= minscore
        scaled[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        scaled[where[0], where[1], 0] -= 1.0
    scaled = numpy.minimum(1.0, numpy.maximum(-1.0, scaled))
    where1 = numpy.where((scaled[:, :, 1] == 1) * (scaled[:, :, 0] >= 0))
    where2 = numpy.where((scaled[:, :, 1] == 1) * (scaled[:, :, 0] < 0))
    temp0 = scaled[where1[0], where1[1], 0]
    temp1 = 1.0 - temp0
    scaled[where1[0], where1[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    temp0 = -scaled[where2[0], where2[1], 0]
    temp1 = 1.0 - temp0
    scaled[where2[0], where2[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    scaled = numpy.round(scaled).astype(numpy.uint32)
    img = numpy.empty((dim, dim), dtype=numpy.uint32)
    img[:, :] = int('ff999999', 16)
    for i in range(dim - 1 + int(diagonal_included)):
        where = numpy.where(scaled[i, :min(scaled.shape[1], dim - i - 1 + int(diagonal_included)), 1] > 0)[0]
        img[i, where + i + 1 - int(diagonal_included)] = scaled[i, where, 0]
        img[where + i + 1 - int(diagonal_included), i] = img[i, where + i + 1 - int(diagonal_included)]
    pilImage = Image.frombuffer('RGBA', (dim, dim), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_full_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                    min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                    returnscale=False, **kwargs):
    """
    Fill in and rescale bitmap from a 5C or HiC full array.

    :param data: A three-dimensional compact array of interaction data.
    :type data: numpy array
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting full array..."),
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    xdim, ydim = data.shape[:2]
    scaled = numpy.copy(data).astype(numpy.float64)
    where = numpy.where(data[:, :, 1] > 0)
    scaled[where[0], where[1], 0] /= scaled[where[0], where[1], 1]
    if logged:
        where = numpy.where(scaled[:, :, 0] <= 0)
        scaled[where[0], where[1], 1] = 0
        where = numpy.where(scaled[:, :, 1] > 0)
        scaled[where[0], where[1], 0] = numpy.log2(scaled[where[0], where[1], 0])
    scaled[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(scaled[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(scaled[where[0], where[1], 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        scaled[where[0], where[1], 0] /= maxscore
    else:
        scaled[where[0], where[1], 0] -= minscore
        scaled[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        scaled[where[0], where[1], 0] -= 1.0
    scaled = numpy.minimum(1.0, numpy.maximum(-1.0, scaled))
    img = numpy.empty((xdim, ydim), dtype=numpy.uint32)
    img.shape = (ydim, xdim)
    img[:, :] = int('ff999999', 16)
    where1 = numpy.where((scaled[:, :, 1] == 1) * (scaled[:, :, 0] >= 0))
    temp0 = scaled[where1[0], where1[1], 0]
    temp1 = 1.0 - temp0
    img[where1[1], where1[0]] = (255 * 256 ** 3 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])).astype(numpy.int64) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])).astype(numpy.int64) * 256 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])).astype(numpy.int64) * 256 ** 2)
    del where1
    where2 = numpy.where((scaled[:, :, 1] == 1) * (scaled[:, :, 0] < 0))
    temp0 = -scaled[where2[0], where2[1], 0]
    temp1 = 1.0 - temp0
    img[where2[1], where2[0]] = (255 * 256 ** 3 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])).astype(numpy.int64) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])).astype(numpy.int64) * 256 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])).astype(numpy.int64) * 256 ** 2)
    pilImage = Image.frombuffer('RGBA', (xdim, ydim), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_upper_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                     min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                     returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a 5C or HiC upper array.

    :param data: A two-dimensional compact array of interaction data.
    :type data: numpy array
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :param diagonal_included: If true, adjust output shape as necessary.
    :type diagonal_included: bool.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting upper array..."),
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    dim = int(0.5 + (0.25 + 2 * data.shape[0]) ** 0.5) - int(diagonal_included)
    scaled = numpy.copy(data).astype(numpy.float64)
    where = numpy.where(data[:, 1] > 0)[0]
    scaled[where, 0] /= scaled[where, 1]
    if logged:
        where = numpy.where(scaled[:, 0] <= 0)[0]
        scaled[where, 1] = 0
        where = numpy.where(scaled[:, 1] > 0)[0]
        scaled[where, 0] = numpy.log2(scaled[where, 0])
    scaled[where, 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(scaled[where, 0])
    if minscore is None:
        minscore = numpy.amin(scaled[where, 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        scaled[where, 0] /= maxscore
    else:
        scaled[where, 0] -= minscore
        scaled[where, 0] /= (maxscore - minscore) * 0.5
        scaled[where, 0] -= 1.0
    scaled = numpy.minimum(1.0, numpy.maximum(-1.0, scaled))
    where1 = numpy.where((scaled[:, 1] == 1) * (scaled[:, 0] >= 0))
    where2 = numpy.where((scaled[:, 1] == 1) * (scaled[:, 0] < 0))
    temp0 = scaled[where1[0], 0]
    temp1 = 1.0 - temp0
    scaled[where1[0], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    temp0 = -scaled[where2[0], 0]
    temp1 = 1.0 - temp0
    scaled[where2[0], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    scaled = numpy.round(scaled).astype(numpy.uint32)
    img = numpy.empty((dim, dim), dtype=numpy.uint32)
    img[:, :] = int('ff999999', 16)
    for i in range(dim - 1 + int(diagonal_included)):
        index = i * dim - i * (i + 1 - 2 * int(diagonal_included)) / 2
        where = numpy.where(scaled[index:(index + dim - i - 1 + int(diagonal_included)), 1] > 0)[0]
        img[i, where + i + 1 - int(diagonal_included)] = scaled[where + index, 0]
        img[where + i + 1 - int(diagonal_included), i] = img[i, where + i + 1 - int(diagonal_included)]
    pilImage = Image.frombuffer('RGBA', (dim, dim), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_hic_heatmap(filename, maxscore=None, minscore=None, symmetricscaling=True, logged=True, chroms=[],
                     min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                     returnscale=False, format='hdf5', **kwargs):
    """
    Fill in and rescale bitmap from a HiC heatmap h5dict file.

    :param filename: File name of heatmap h5dict containing binned data arrays.
    :type filename: str.
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param chroms: A list of chromosome names to include in the plot. If left empty, all chromosomes present in the heatmap file will be plotted.
    :type chroms: list
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :param format: Format of the heatmap.
    :type format: str.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.\n")
        return None
    if format not in ['hdf5', 'npz']:
        if not silent:
            print >> sys.stderr, ("The format is not recognized or supported.\n")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting heatmap dict..."),
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    if format == 'hdf5':
        input = h5py.File(filename, 'r')
    else:
        input = numpy.load(filename)
    if len(chroms) == 0:
        chroms = list(input['chromosomes'][:])
    starts = [0]
    sizes = []
    for chrom in chroms:
        sizes.append(input['%s.positions' % chrom].shape[0])
        starts.append(starts[-1] + sizes[-1] + 1)
    data = numpy.zeros((starts[-1] - 1, starts[-1] - 1, 2), dtype=numpy.float64)
    for i in range(len(chroms)):
        if '%s.counts' % chroms[i] in input and '%s.expected' % chroms[i] in input:
            n = input['%s.positions' % chroms[i]].shape[0]
            if n * (n - 1) / 2 == input['%s.counts' % chroms[i]].shape[0]:
                indices = numpy.triu_indices(sizes[i], 1)
            else:
                indices = numpy.triu_indices(sizes[i], 0)
            data[indices[0] + starts[i], indices[1] + starts[i], 0] = input['%s.counts' % chroms[i]][:]
            data[indices[0] + starts[i], indices[1] + starts[i], 1] = input['%s.expected' % chroms[i]][:]
    for i in range(len(chroms) - 1):
        for j in range(i + 1, len(chroms)):
            if '%s_by_%s.counts' % (chroms[i], chroms[j]) in input.keys():
                data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), 0] = (
                    input['%s_by_%s.counts' % (chroms[i], chroms[j])][:, :])
                data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), 1] = (
                    input['%s_by_%s.expected' % (chroms[i], chroms[j])][:, :])
                data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), :] = numpy.transpose(
                    data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), :], axes=[1, 0, 2])
            elif '%s_by_%s.counts' % (chroms[j], chroms[i]) in input.keys():
                data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), 0] = (
                    input['%s_by_%s.counts' % (chroms[j], chroms[i])][:, :])
                data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), 1] = (
                    input['%s_by_%s.expected' % (chroms[j], chroms[i])][:, :])
                data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), :] = numpy.transpose(
                    data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), :], axes=[1, 0, 2])
    where = numpy.where(data[:, :, 1] > 0)
    data[where[0], where[1], 0] /= data[where[0], where[1], 1]
    if logged:
        where = numpy.where(data[:, :, 0] <= 0)
        data[where[0], where[1], 1] = 0
        where = numpy.where(data[:, :, 1] > 0)
        data[where[0], where[1], 0] = numpy.log2(data[where[0], where[1], 0])
    data[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(data[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(data[where[0], where[1], 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        data[where[0], where[1], 0] /= maxscore
    else:
        data[where[0], where[1], 0] -= minscore
        data[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        data[where[0], where[1], 0] -= 1.0
    data = numpy.minimum(1.0, numpy.maximum(-1.0, data))
    where1 = numpy.where((data[:, :, 1] == 1) * (data[:, :, 0] >= 0))
    where2 = numpy.where((data[:, :, 1] == 1) * (data[:, :, 0] < 0))
    temp0 = data[where1[0], where1[1], 0]
    temp1 = 1.0 - temp0
    data[where1[0], where1[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    temp0 = -data[where2[0], where2[1], 0]
    temp1 = 1.0 - temp0
    data[where2[0], where2[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    data = numpy.round(data).astype(numpy.uint32)
    img = numpy.empty((data.shape[0], data.shape[0]), dtype=numpy.uint32)
    img[:, :] = int('ff999999', 16)
    where = numpy.where(data[:, :, 1] > 0)
    img[where] = data[where[0], where[1], 0]
    img[where[1], where[0]] = img[where]
    pilImage = Image.frombuffer('RGBA', (data.shape[0], data.shape[0]), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_fivec_heatmap(filename, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                       regions=[], min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                       returnscale=False, **kwargs):
    """
    Fill in and rescale bitmap in a full format from a 5C heatmap h5dict file.

    :param filename: Location of a heatmap h5dict containing 5C data arrays.
    :type data: str.
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param regions: If specified, only the indicated regions are plotted. Otherwise all regions present in the h5dict are plotted.
    :type regions: list
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting heatmap dict..."),
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    input = h5py.File(filename, 'r')
    if len(regions) == 0:
        regions = []
        for key in input.keys():
            if key.count('positions') > 0 and key.count('reverse') == 0:
                regions.append(int(key.split('.')[0]))
        regions.sort()
    else:
        for i in range(len(regions))[::-1]:
            if ("%i.positions" % regions[i] not in input.keys() and
                "%i.forward_positions" % regions[i] not in input.keys()):
                if not silent:
                    print >> sys.stderr, ("Region %i not in heatmap, skipping\n"),
                del regions[i]
    if len(regions) == 0:
        if not silent:
            print >> sys.stderr, ("No valid data found.\n"),
        return None
    if "%i.positions" % regions[0] in input.keys():
        starts = [0]
    else:
        starts = [[0, 0]]
    sizes = []
    for region in regions:
        if "%i.positions" % region in input:
            sizes.append(input['%i.positions' % region].shape[0])
            starts.append(starts[-1] + sizes[-1] + 1)
        else:
            sizes.append([input['%i.forward_positions' % region].shape[0],
                          input['%i.reverse_positions' % region].shape[0]])
            starts.append([starts[-1][0] + sizes[-1][0] + 1, starts[-1][1] + sizes[-1][1] + 1])
    sizes = numpy.array(sizes, dtype=numpy.int32)
    starts = numpy.array(starts, dtype=numpy.int32)
    if len(sizes.shape) == 1:
        data = numpy.zeros((starts[-1] - 1, starts[-1] - 1, 2), dtype=numpy.float64)
    else:
        data = numpy.zeros((starts[-1, 0] - 1, starts[-1, 1] - 1, 2), dtype=numpy.float64)
    for i in range(len(regions)):
        if '%i.counts' % regions[i] in input and '%i.expected' % regions[i] in input:
            if len(sizes.shape) == 1:
                data[starts[i]:(starts[i + 1] - 1),
                     starts[i]:(starts[i + 1]  - 1), 0] = input['%i.counts' % regions[i]][...]
                data[starts[i]:(starts[i + 1] - 1),
                     starts[i]:(starts[i + 1] - 1), 1] = input['%i.expected' % regions[i]][...]
            else:
                data[starts[i, 0]:(starts[i + 1, 0] - 1),
                     starts[i, 1]:(starts[i + 1, 1]  - 1), 0] = input['%i.counts' % regions[i]][...]
                data[starts[i, 0]:(starts[i + 1, 0] - 1),
                     starts[i, 1]:(starts[i + 1, 1] - 1), 1] = input['%i.expected' % regions[i]][...]
    for i in range(len(regions) - 1):
        for j in range(i + 1, len(regions)):
            if '%i_by_%i.counts' % (regions[i], regions[j]) in input.keys():
                if len(sizes.shape) == 1:
                    data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), 0] = (
                        input['%i_by_%i.counts' % (regions[i], regions[j])][:, :])
                    data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), 1] = (
                        input['%i_by_%i.expected' % (regions[i], regions[j])][:, :])
                    data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), :] = numpy.transpose(
                        data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), :], axes=[1, 0, 2])
                else:
                    data[starts[i, 0]:(starts[i + 1, 0] - 1), starts[j, 1]:(starts[j + 1, 1] - 1), 0] = (
                        input['%i_by_%i.counts' % (regions[i], regions[j])][:, :])
                    data[starts[i, 0]:(starts[i + 1, 0] - 1), starts[j, 1]:(starts[j + 1, 1] - 1), 1] = (
                        input['%i_by_%i.expected' % (regions[i], regions[j])][:, :])
            if '%i_by_%i.counts' % (regions[j], regions[i]) in input.keys():
                if len(sizes.shape) == 1:
                    data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), 0] = (
                        input['%i_by_%i.counts' % (regions[j], regions[i])][:, :])
                    data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), 1] = (
                        input['%i_by_%i.expected' % (regions[j], regions[i])][:, :])
                    data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), :] = numpy.transpose(
                        data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), :], axes=[1, 0, 2])
                else:
                    data[starts[j, 0]:(starts[j + 1, 0] - 1), starts[i, 1]:(starts[i + 1, 1] - 1), 0] = (
                        input['%i_by_%i.counts' % (regions[j], regions[i])][:, :])
                    data[starts[j, 0]:(starts[j + 1, 0] - 1), starts[i, 1]:(starts[i + 1, 1] - 1), 1] = (
                        input['%i_by_%i.expected' % (regions[j], regions[i])][:, :])
    where = numpy.where(data[:, :, 1] > 0)
    data[where[0], where[1], 0] /= data[where[0], where[1], 1]
    if logged:
        where = numpy.where(data[:, :, 0] <= 0)
        data[where[0], where[1], 1] = 0
        where = numpy.where(data[:, :, 1] > 0)
        data[where[0], where[1], 0] = numpy.log2(data[where[0], where[1], 0])
    data[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(data[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(data[where[0], where[1], 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        data[where[0], where[1], 0] /= maxscore
    else:
        data[where[0], where[1], 0] -= minscore
        data[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        data[where[0], where[1], 0] -= 1.0
    data = numpy.minimum(1.0, numpy.maximum(-1.0, data))
    where1 = numpy.where((data[:, :, 1] == 1) * (data[:, :, 0] >= 0))
    where2 = numpy.where((data[:, :, 1] == 1) * (data[:, :, 0] < 0))
    temp0 = data[where1[0], where1[1], 0]
    temp1 = 1.0 - temp0
    data[where1[0], where1[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    temp0 = -data[where2[0], where2[1], 0]
    temp1 = 1.0 - temp0
    data[where2[0], where2[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    data = numpy.round(data).astype(numpy.uint32)
    img = numpy.empty((data.shape[0], data.shape[1]), dtype=numpy.uint32)
    img.shape = (img.shape[1], img.shape[0])
    img[:, :] = int('ff999999', 16)
    where = numpy.where(data[:, :, 1] > 0)
    img[where[1], where[0]] = data[where[0], where[1], 0]
    black = int('ff000000', 16)
    if len(sizes.shape) == 1:
        for i in range(1, starts.shape[0] - 1):
            img[starts[i] - 1, :] = black
            img[:, starts[i] - 1] = black
    else:
        for i in range(1, starts.shape[0] - 1):
            img[starts[i, 1] - 1, :] = black
            img[:, starts[i, 0] - 1] = black
    pilImage = Image.frombuffer('RGBA', (data.shape[0], data.shape[1]), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_fivec_compact_heatmap_dict(filename, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                    regions=[], min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                                    returnscale=False, **kwargs):
    """
    Fill in and rescale bitmap in a compact from a 5C heatmap h5dict file.

    This plots the data in a 5C compact format such that the rows correspond to positive-strand primers and columns correspond to negative-strand primers.

    :param filename: Location of a heatmap h5dict containing 5C data arrays.
    :type data: str.
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param regions: If specified, only the indicated regions are plotted. Otherwise all regions present in the h5dict are plotted.
    :type regions: list
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting heatmap dict..."),
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    input = h5py.File(filename, 'r')
    if len(regions) == 0:
        regions = []
        for key in input.keys():
            if key.count('forward_fragments') > 0:
                regions.append(int(key.split('.')[0]))
        regions.sort()
    xstarts = [0]
    ystarts = [0]
    xsizes = []
    ysizes = []
    for region in regions:
        xsizes.append(input['%i.forward_fragments' % region].shape[0])
        xstarts.append(xstarts[-1] + xsizes[-1] + 1)
        ysizes.append(input['%i.reverse_fragments' % region].shape[0])
        ystarts.append(ystarts[-1] + ysizes[-1] + 1)
    data = numpy.zeros((xstarts[-1] - 1, ystarts[-1] - 1, 2), dtype=numpy.float64)
    for i in range(len(regions)):
        if '%i.counts' % regions[i] in input and '%i.expected' % regions[i] in input:
            data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[i]:(ystarts[i + 1]  - 1), 0] = (
                    input['%i.counts' % regions[i]][...])
            data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[i]:(ystarts[i + 1]  - 1), 1] = (
                    input['%i.expected' % regions[i]][...])
    for i in range(len(regions) - 1):
        for j in range(i + 1, len(regions)):
            if ('%i_by_%i.counts' % (regions[i], regions[j]) in input and
                '%i_by_%i.expected' % (regions[i], regions[j]) in input):
                data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[j]:(ystarts[j + 1] - 1), 0] = (
                    input['%i_by_%i.counts' % (regions[i], regions[j])][:, :])
                data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[j]:(ystarts[j + 1] - 1), 1] = (
                    input['%i_by_%i.expected' % (regions[i], regions[j])][:, :])
            if ('%i_by_%i.counts' % (regions[j], regions[i]) in input and
                '%i_by_%i.expected' % (regions[j], regions[i]) in input):
                data[xstarts[j]:(xstarts[j + 1] - 1), ystarts[i]:(ystarts[i + 1] - 1), 0] = (
                    input['%i_by_%i.counts' % (regions[j], regions[i])][:, :])
                data[xstarts[j]:(xstarts[j + 1] - 1), ystarts[i]:(ystarts[i + 1] - 1), 1] = (
                    input['%i_by_%i.expected' % (regions[j], regions[i])][:, :])
    where = numpy.where(data[:, :, 1] > 0)
    data[where[0], where[1], 0] /= data[where[0], where[1], 1]
    if logged:
        where = numpy.where(data[:, :, 0] <= 0)
        data[where[0], where[1], 1] = 0
        where = numpy.where(data[:, :, 1] > 0)
        data[where[0], where[1], 0] = numpy.log2(data[where[0], where[1], 0])
    data[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(data[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(data[where[0], where[1], 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        data[where[0], where[1], 0] /= maxscore
    else:
        data[where[0], where[1], 0] -= minscore
        data[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        data[where[0], where[1], 0] -= 1.0
    data = numpy.minimum(1.0, numpy.maximum(-1.0, data))
    where1 = numpy.where((data[:, :, 1] == 1) * (data[:, :, 0] >= 0))
    where2 = numpy.where((data[:, :, 1] == 1) * (data[:, :, 0] < 0))
    temp0 = data[where1[0], where1[1], 0]
    temp1 = 1.0 - temp0
    data[where1[0], where1[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    temp0 = -data[where2[0], where2[1], 0]
    temp1 = 1.0 - temp0
    data[where2[0], where2[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    data = numpy.round(data).astype(numpy.uint32)
    img = numpy.empty((data.shape[0], data.shape[1]), dtype=numpy.uint32)
    img.shape = (data.shape[1], data.shape[0])
    img[:, :] = int('ff999999', 16)
    where = numpy.where(data[:, :, 1] > 0)
    img[where[1], where[0]] = data[where[0], where[1], 0]
    black = int('ff000000', 16)
    for i in range(1, len(xstarts) - 1):
        img[:, xstarts[i] - 1] = black
    for i in range(1, len(ystarts) - 1):
        img[ystarts[i] - 1, :] = black
    pilImage = Image.frombuffer('RGBA', (data.shape[0], data.shape[1]), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_diagonal_from_upper_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                   min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                                   returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a HiC upper array, plotting only the upper triangle rotated 45 degrees counter-clockwise.

    :param data: A three-dimensional upper array of HiC interaction data.
    :type data: numpy array
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :param diagonal_included: If true, adjust output shape as necessary.
    :type diagonal_included: bool.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting rotated compact array..."),
    diag = int(diagonal_included)
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    dim = int(0.5 + (0.25 + data.shape[0] * 2) ** 0.5) - diag
    rotated = numpy.zeros((dim * 2 - 1 + diag * 2, dim + diag, 2), dtype=numpy.float64)
    xdim, ydim = rotated.shape[:2]
    for i in range(dim - 1 + diag):
        start = i * dim - (i * (i + 1 - 2 * diag) / 2)
        stop = start + dim - i - 1 + diag
        span = stop - start
        rotated[numpy.arange(i * 2 + 1, i * 2 + span + 1), ydim - numpy.arange(span) - 1, :] += data[start:stop, :]
        rotated[numpy.arange(i * 2, i * 2 + span), ydim - numpy.arange(span) - 1, :] += data[start:stop, :] * 0.25
        rotated[numpy.arange(i * 2 + 2, i * 2 + span + 2), ydim - numpy.arange(span) - 1, :] += data[start:stop, :] * 0.25
        rotated[numpy.arange(i * 2 + 1, i * 2 + span + 1), ydim - numpy.arange(span) - 2, :] += data[start:stop, :] * 0.25
        rotated[numpy.arange(i * 2 + 2, i * 2 + span + 1), ydim - numpy.arange(span - 1) - 1, :] += data[(start + 1):stop, :] * 0.25
    where = numpy.where(rotated[:, :, 1] > 0)
    rotated[where[0], where[1], 0] /= rotated[where[0], where[1], 1]
    if logged:
        where = numpy.where(rotated[:, :, 0] <= 0)
        rotated[where[0], where[1], 1] = 0
        where = numpy.where(rotated[:, :, 1] > 0)
        rotated[where[0], where[1], 0] = numpy.log2(rotated[where[0], where[1], 0])
    rotated[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(rotated[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(rotated[where[0], where[1], 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        rotated[where[0], where[1], 0] /= maxscore
    else:
        rotated[where[0], where[1], 0] -= minscore
        rotated[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        rotated[where[0], where[1], 0] -= 1.0
    rotated = numpy.minimum(1.0, numpy.maximum(-1.0, rotated))
    where1 = numpy.where((rotated[:, :, 1] == 1) * (rotated[:, :, 0] >= 0))
    where2 = numpy.where((rotated[:, :, 1] == 1) * (rotated[:, :, 0] < 0))
    temp0 = rotated[where1[0], where1[1], 0]
    temp1 = 1.0 - temp0
    rotated[where1[0], where1[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    temp0 = -rotated[where2[0], where2[1], 0]
    temp1 = 1.0 - temp0
    rotated[where2[0], where2[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    rotated = numpy.round(rotated).astype(numpy.uint32)
    img = numpy.empty((xdim, ydim), dtype=numpy.uint32)
    img.shape = (ydim, xdim)
    img[:, :] = int('ff999999', 16)
    where = numpy.where(rotated[:, :, 1] > 0)
    img[where[1], where[0]] = rotated[where[0], where[1], 0]
    pilImage = Image.frombuffer('RGBA', (xdim, ydim), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_diagonal_from_compact_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                     min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                                     returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a HiC compact array, plotting only the upper triangle rotated 45 degrees counter-clockwise.

    :param data: A three-dimensional compact array of HiC interaction data.
    :type data: numpy array
    :param maxscore: A ceiling value to cutoff scores at for plot color.
    :type maxscore: float
    :param minscore: A floor value to cutoff scores at for plot color.
    :type minscore: float
    :param symmetricscaling: Indicates whether to recenter data for scaling or maintain scores about zero.
    :type symmetricscaling: bool.
    :param logged: Indicates whether to use log2 values of scores for color values.
    :type logged: bool.
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param returnscale: Indicates whether to return a list containing the bitmap, minimum score, and maximum score, or just the bitmap.
    :type returnscale: bool.
    :param diagonal_included: If true, adjust output shape as necessary.
    :type diagonal_included: bool.
    :returns: :mod:`PIL` bitmap object and if requested, minimum and maximum scores.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    if not silent:
        print >> sys.stderr, ("Plotting rotated compact array..."),
    diag = int(diagonal_included)
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    rotated = numpy.zeros((data.shape[0] * 2 - 1 + diag * 2, data.shape[1] + 1, 2), dtype=numpy.float64)
    xdim, ydim = rotated.shape[:2]
    for i in range(data.shape[1]):
        rotated[(i + 1):(2 * data.shape[0] - i - 1 + diag * 2):2, ydim - i - 1, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :])
        rotated[i:(2 * data.shape[0] - i - 2 + diag * 2):2, ydim - i - 1, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
        rotated[(i + 2):(2 * data.shape[0] - i + diag * 2):2, ydim - i - 1, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
        rotated[(i + 1):(2 * data.shape[0] - i - 1 + diag * 2):2, ydim - i - 2, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
        if i > 0:
            rotated[(i + 1):(2 * data.shape[0] - i - 1 + diag * 2):2, ydim - i, :] += (
                data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
    where = numpy.where(rotated[:, :, 1] > 0)
    rotated[where[0], where[1], 0] /= rotated[where[0], where[1], 1]
    if logged:
        where = numpy.where(rotated[:, :, 0] <= 0)
        rotated[where[0], where[1], 1] = 0
        where = numpy.where(rotated[:, :, 1] > 0)
        rotated[where[0], where[1], 0] = numpy.log2(rotated[where[0], where[1], 0])
    rotated[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(rotated[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(rotated[where[0], where[1], 0])
    if symmetricscaling:
        maxscore = max(abs(maxscore), abs(minscore))
        minscore = -maxscore
        rotated[where[0], where[1], 0] /= maxscore
    else:
        rotated[where[0], where[1], 0] -= minscore
        rotated[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        rotated[where[0], where[1], 0] -= 1.0
    rotated = numpy.minimum(1.0, numpy.maximum(-1.0, rotated))
    where1 = numpy.where((rotated[:, :, 1] == 1) * (rotated[:, :, 0] >= 0))
    where2 = numpy.where((rotated[:, :, 1] == 1) * (rotated[:, :, 0] < 0))
    temp0 = rotated[where1[0], where1[1], 0]
    temp1 = 1.0 - temp0
    rotated[where1[0], where1[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    temp0 = -rotated[where2[0], where2[1], 0]
    temp1 = 1.0 - temp0
    rotated[where2[0], where2[1], 0] = (255.0 * 256.0 ** 3.0 +
        numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
        numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256.0 +
        numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256.0 ** 2.0)
    rotated = numpy.round(rotated).astype(numpy.uint32)
    img = numpy.empty((xdim, ydim), dtype=numpy.uint32)
    img.shape = (ydim, xdim)
    img[:, :] = int('ff999999', 16)
    where = numpy.where(rotated[:, :, 1] > 0)
    img[where[1], where[0]] = rotated[where[0], where[1], 0]
    pilImage = Image.frombuffer('RGBA', (xdim, ydim), img, 'raw', 'RGBA', 0, 1)
    if not silent:
        print >> sys.stderr, ("Done\n"),
    if returnscale:
        return [pilImage, minscore, maxscore]
    else:
        return pilImage


def plot_key(min_score, max_score, height, width, labelformat='%0.2f', orientation='left', num_ticks=5,
             min_color="0000ff", mid_color="ffffff", max_color="ff0000", labelattr=None, log_display=True, **kwargs):
    """
    Create a key including color gradient and labels indicating associated values, returning a :mod:`pyx` canvas.

    :param min_score: The minimum value of the key scale.
    :type min_score: float
    :param max_score: The maximum value of the key scale.
    :type max_score: float
    :param height: The height of the gradient bar in whatever units :mod:`pyx` is using.
    :type height: float
    :param width: The width of the gradient bar in whatever units :mod:`pyx` is using.
    :type width: float
    :param labelformat: A string denoting the format for displaying number labels using the string formatting style from Python <= 2.6.
    :type labelformat: str.
    :param orientation: Indicates where labels are placed relative to gradient bar. This parameter will accept 'left', 'right', 'top', and 'bottom'.
    :type type: str.
    :param num_ticks: Indicates how many evenly-spaced tick marks and associated labels to insert. This can be zero for no labels or greater than one. Labels are inserted at the minimum and maximum values first with remaining ticks occuring evenly distributed between the extremes.
    :type num_ticks: int.
    :param min_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This variable is used to create a color gradient for plotting along with max_color and mid_color.
    :type min_color: str.
    :param mid_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the minimum plot value. This can be set to None to create a gradient ranging from min_color to max_color or to a hex color to create a divergent gradient.
    :type mid_color: str.
    :param max_color: This is a hex color code ("rrggbb" where each pair ranges from 00-ff) specifying the color associated with the maximum plot value. This variable is used to create a color gradient for plotting along with min_color and mid_color.
    :type max_color: str.
    :param labelattr: A list of pyx attributes to be passed to the text function.
    :type labelattr: str.
    :param log_display: If True, min_score and max_score are taken to be logged values and so labels are evenly spaced in log space but converted to normal space for display.
    :type log_display: bool.
    :returns: :mod:`Pxy` canvas object.
    """
    if 'silent' in kwargs and kwargs['silent']:
        silent = True
    else:
        silent = False
    if 'PIL' not in sys.modules.keys() or 'pyx' not in sys.modules.keys():
        if not silent:
            print >> sys.stderr, ("The PIL and pyx modules must be installed to use this function.")
        return None
    height = float(height)
    width = float(width)
    min_score = float(min_score)
    max_score = float(max_score)
    c = canvas.canvas()
    if orientation in ['right', 'left']:
        img = numpy.zeros( (int(round(511.0 * width / height)), 511), dtype=numpy.uint32 )
    else:
        img = numpy.zeros( (int(round(511.0 * height / width)), 511), dtype=numpy.uint32 )
    img.shape = (img.shape[1],img.shape[0])
    min_color = (int(min_color[:2], 16) / 255.0, int(min_color[2:4], 16) / 255.0, int(min_color[4:6], 16) / 255.0)
    max_color = (int(max_color[:2], 16) / 255.0, int(max_color[2:4], 16) / 255.0, int(max_color[4:6], 16) / 255.0)
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[0] + max_color[0]) / 2.0,
                     (min_color[0] + max_color[0]) / 2.0)
    else:
        mid_color = (int(mid_color[:2], 16) / 255.0, int(mid_color[2:4], 16) / 255.0, int(mid_color[4:6], 16) / 255.0)
    temp0 = numpy.linspace(1.0, 0.0, 256)[:255].astype(numpy.float64).reshape(-1, 1)
    temp1 = 1.0 - temp0
    img[:255, :] = (255 * 256 ** 3.0 + numpy.round(255 * (temp0 * max_color[0] + temp1 * mid_color[0])) +
                    numpy.round(255 * (temp0 * max_color[1] + temp1 * mid_color[1])) * 256 +
                    numpy.round(255 * (temp0 * max_color[2] + temp1 * mid_color[2])) * 256 ** 2).astype(numpy.uint32)
    temp0 = numpy.linspace(0.0, 1.0, 256).astype(numpy.float64).reshape(-1, 1)
    temp1 = 1.0 - temp0
    img[255:, :] = (255 * 256 ** 3.0 + numpy.round(255 * (temp0 * min_color[0] + temp1 * mid_color[0])) +
                    numpy.round(255 * (temp0 * min_color[1] + temp1 * mid_color[1])) * 256 +
                    numpy.round(255 * (temp0 * min_color[2] + temp1 * mid_color[2])) * 256 ** 2).astype(numpy.uint32)
    pilImage = Image.frombuffer( 'RGBA',(img.shape[1],img.shape[0]),img,'raw','RGBA',0,1)
    if orientation in ['left', 'right']:
        c.insert(bitmap.bitmap( 0, 0, pilImage, height=height))
    else:
        c.insert(bitmap.bitmap( 0, 0, pilImage, height=width), [trafo.rotate(-90), trafo.translate(0, height)])
    c.stroke( path.rect( 0, 0, width, height ), [style.linewidth.THin] )
    labels = []
    if num_ticks > 1:
        ticks = numpy.zeros((num_ticks, 4), dtype=numpy.float32)
        if log_display:
            for i in range(num_ticks):
                labels.append(labelformat % exp(min_score + i / (num_ticks - 1.0) * (max_score - min_score)))
        else:
            for i in range(num_ticks):
                labels.append(labelformat % (min_score + i / (num_ticks - 1.0) * (max_score - min_score)))
        if orientation == 'left':
            ypos = numpy.linspace(0.0, height, num_ticks)
            xpos = numpy.ones(num_ticks, dtype=numpy.float32) * -width * 0.5
            ticks[:, 2] = -width * 0.4
            ticks[:, 1] = numpy.linspace(0, height, num_ticks)
            ticks[:, 3] = numpy.linspace(0, height, num_ticks)
            if labelattr is None:
                labelattr = [text.halign.right, text.valign.middle, text.size(-3)]
        elif orientation == 'right':
            ypos = numpy.linspace(0.0, height, num_ticks)
            xpos = numpy.ones(num_ticks, dtype=numpy.float32) * width * 1.5
            ticks[:, 0] = width
            ticks[:, 2] = width * 1.4
            ticks[:, 1] = numpy.linspace(0, height, num_ticks)
            ticks[:, 3] = numpy.linspace(0, height, num_ticks)
            if labelattr is None:
                labelattr = [text.halign.left, text.valign.middle, text.size(-3)]
        elif orientation == 'bottom':
            xpos = numpy.linspace(0.0, width, num_ticks)
            ypos = numpy.ones(num_ticks, dtype=numpy.float32) * -height * 0.5
            ticks[:, 3] = -height * 0.4
            ticks[:, 0] = numpy.linspace(0, width, num_ticks)
            ticks[:, 2] = numpy.linspace(0, width, num_ticks)
            if labelattr is None:
                labelattr = [text.halign.center, text.valign.top, text.size(-3)]
        else:
            xpos = numpy.linspace(0.0, width, num_ticks)
            ypos = numpy.ones(num_ticks, dtype=numpy.float32) * height * 1.5
            ticks[:, 1] = height
            ticks[:, 3] = height * 1.4
            ticks[:, 0] = numpy.linspace(0, width, num_ticks)
            ticks[:, 2] = numpy.linspace(0, width, num_ticks)
            if labelattr is None:
                labelattr = [text.halign.center, text.valign.bottom, text.size(-3)]
        for i in range(num_ticks):
            c.stroke( path.line(ticks[i, 0], ticks[i, 1], ticks[i, 2], ticks[i, 3]), [style.linewidth.THin] )
            c.text(xpos[i], ypos[i], r"%s" % labels[i], labelattr)
    return c
