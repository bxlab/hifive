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
import matplotlib.pyplot as plt
import matplotlib.colors
try:
    from PIL import Image
except:
    pass
try:
    from pyx import canvas, path, style, bitmap, text, trafo
except:
    pass


def plot_compact_array(fname, data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                       min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                       returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a HiC compact array.

    :param fname: A file name to save plot to.
    :type fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting compact array...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
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
    if not diagonal_included:
        scaled[numpy.arange(scores.shape[0]), numpy.arange(scores.shape[0]), 0] = numpy.nan
    fig, ax = plt.subplots(1, 1, figsize=(scaled.shape[0]/100, scaled.shape[0]/100))
    ax.imshow(scaled[:, :, 0], norm=plt.Normalize(-1, 1), cmap=cmap)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_full_array(fname, data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                    min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                    returnscale=False, **kwargs):
    """
    Fill in and rescale bitmap from a 5C or HiC full array.

    :param fname: A file name to save plot to.
    :type fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting full array...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
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
    fig, ax = plt.subplots(1, 1, figsize=(scaled.shape[0]/100, scaled.shape[1]/100))
    ax.imshow(scaled[:, :, 0].T, norm=plt.Normalize(-1, 1), cmap=cmap)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_upper_array(fname, data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                     min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                     returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a 5C or HiC upper array.

    :param fname: A file name to save plot to.
    :type fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting upper array...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
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
    img = numpy.full((dim, dim), numpy.nan, dtype=numpy.float32)
    for i in range(dim - 1 + int(diagonal_included)):
        index = i * dim - i * (i + 1 - 2 * int(diagonal_included)) / 2
        where = numpy.where(scaled[index:(index + dim - i - 1 + int(diagonal_included)), 1] > 0)[0]
        img[i, where + i + 1 - int(diagonal_included)] = scaled[where + index, 0]
        img[where + i + 1 - int(diagonal_included), i] = img[i, where + i + 1 - int(diagonal_included)]
    fig, ax = plt.subplots(1, 1, figsize=(img.shape[0]/100, img.shape[1]/100))
    ax.imshow(img.T, norm=plt.Normalize(-1, 1), cmap=cmap)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_hic_heatmap(fname, in_fname, maxscore=None, minscore=None, symmetricscaling=True, logged=True, chroms=[],
                     min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                     returnscale=False, format='hdf5', **kwargs):
    """
    Fill in and rescale bitmap from a HiC heatmap h5dict file.

    :param fname: A file name to save plot to.
    :type fname: str.
    :param in_fname: File name of heatmap h5dict or numpy npz containing binned data arrays.
    :type in_fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if format not in ['hdf5', 'npz']:
        if not silent:
            print("The format is not recognized or supported.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting heatmap dict...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
    if format == 'hdf5':
        input = h5py.File(in_fname, 'r')
        if len(chroms) == 0:
            chroms = [x.decode('utf8') for x in input['chromosomes'][:]]
    else:
        input = numpy.load(in_fname)
        if len(chroms) == 0:
            chroms = [x.decode('utf8') for x in input['chromosomes']]
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
    img = numpy.full((data.shape[0], data.shape[0]), numpy.nan, dtype=numpy.float32)
    where = numpy.where(data[:, :, 1] > 0)
    img[where] = data[where[0], where[1], 0]
    img[where[1], where[0]] = img[where]
    fig, ax = plt.subplots(1, 1, figsize=(img.shape[0]/100, img.shape[1]/100))
    ax.imshow(img.T, norm=plt.Normalize(-1, 1), cmap=cmap)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_fivec_heatmap(fname, hdf_fname, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                       regions=[], min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                       returnscale=False, **kwargs):
    """
    Fill in and rescale bitmap in a full format from a 5C heatmap h5dict file.

    :param fname: A file name to save plot to.
    :type fname: str.
    :param hdf_fname: Location of a heatmap h5dict containing 5C data arrays.
    :type hdf_fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", end='', file=sys.stderr)
        return None
    if not silent:
        print("Plotting heatmap dict...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
    input = h5py.File(hdf_fname, 'r')
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
                    print("Region %i not in heatmap, skipping" % regions[i], file=sys.stderr)
                del regions[i]
    if len(regions) == 0:
        if not silent:
            print("No valid data found.", file=sys.stderr)
        return None
    if "%i.positions" % regions[0] in input.keys():
        starts = [0]
    else:
        starts = [[0, 0]]
    sizes = []
    for region in regions:
        if "%i.positions" % region in input:
            sizes.append(input['%i.positions' % region].shape[0])
            starts.append(starts[-1] + sizes[-1])
        else:
            sizes.append([input['%i.forward_positions' % region].shape[0],
                          input['%i.reverse_positions' % region].shape[0]])
            starts.append([starts[-1][0] + sizes[-1][0], starts[-1][1] + sizes[-1][1]])
    sizes = numpy.array(sizes, dtype=numpy.int32)
    starts = numpy.array(starts, dtype=numpy.int32)
    if len(sizes.shape) == 1:
        data = numpy.zeros((starts[-1], starts[-1], 2), dtype=numpy.float64)
    else:
        data = numpy.zeros((starts[-1, 0], starts[-1, 1], 2), dtype=numpy.float64)
    for i in range(len(regions)):
        if '%i.counts' % regions[i] in input and '%i.expected' % regions[i] in input:
            if len(sizes.shape) == 1:
                data[starts[i]:(starts[i + 1]),
                     starts[i]:(starts[i + 1]), 0] = input['%i.counts' % regions[i]][...]
                data[starts[i]:(starts[i + 1]),
                     starts[i]:(starts[i + 1]), 1] = input['%i.expected' % regions[i]][...]
            else:
                data[starts[i, 0]:(starts[i + 1, 0]),
                     starts[i, 1]:(starts[i + 1, 1]), 0] = input['%i.counts' % regions[i]][...]
                data[starts[i, 0]:(starts[i + 1, 0]),
                     starts[i, 1]:(starts[i + 1, 1]), 1] = input['%i.expected' % regions[i]][...]
    for i in range(len(regions) - 1):
        for j in range(i + 1, len(regions)):
            if '%i_by_%i.counts' % (regions[i], regions[j]) in input.keys():
                if len(sizes.shape) == 1:
                    data[starts[i]:(starts[i + 1]), starts[j]:(starts[j + 1]), 0] = (
                        input['%i_by_%i.counts' % (regions[i], regions[j])][:, :])
                    data[starts[i]:(starts[i + 1]), starts[j]:(starts[j + 1]), 1] = (
                        input['%i_by_%i.expected' % (regions[i], regions[j])][:, :])
                    data[starts[j]:(starts[j + 1]), starts[i]:(starts[i + 1]), :] = numpy.transpose(
                        data[starts[i]:(starts[i + 1]), starts[j]:(starts[j + 1]), :], axes=[1, 0, 2])
                else:
                    data[starts[i, 0]:(starts[i + 1, 0]), starts[j, 1]:(starts[j + 1, 1]), 0] = (
                        input['%i_by_%i.counts' % (regions[i], regions[j])][:, :])
                    data[starts[i, 0]:(starts[i + 1, 0]), starts[j, 1]:(starts[j + 1, 1]), 1] = (
                        input['%i_by_%i.expected' % (regions[i], regions[j])][:, :])
            if '%i_by_%i.counts' % (regions[j], regions[i]) in input.keys():
                if len(sizes.shape) == 1:
                    data[starts[j]:(starts[j + 1]), starts[i]:(starts[i + 1]), 0] = (
                        input['%i_by_%i.counts' % (regions[j], regions[i])][:, :])
                    data[starts[j]:(starts[j + 1]), starts[i]:(starts[i + 1]), 1] = (
                        input['%i_by_%i.expected' % (regions[j], regions[i])][:, :])
                    data[starts[i]:(starts[i + 1]), starts[j]:(starts[j + 1]), :] = numpy.transpose(
                        data[starts[j]:(starts[j + 1]), starts[i]:(starts[i + 1]), :], axes=[1, 0, 2])
                else:
                    data[starts[j, 0]:(starts[j + 1, 0]), starts[i, 1]:(starts[i + 1, 1]), 0] = (
                        input['%i_by_%i.counts' % (regions[j], regions[i])][:, :])
                    data[starts[j, 0]:(starts[j + 1, 0]), starts[i, 1]:(starts[i + 1, 1]), 1] = (
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
    img = numpy.full((data.shape[0], data.shape[1]), numpy.nan, dtype=numpy.float32)
    where = numpy.where(data[:, :, 1] > 0)
    img[where] = data[where[0], where[1], 0]
    fig, ax = plt.subplots(1, 1, figsize=(img.shape[0]/100, img.shape[1]/100))
    ax.imshow(img.T, norm=plt.Normalize(-1, 1), cmap=cmap)
    if len(starts.shape) == 1:
        for i in range(1, starts.shape[0] - 1):
            ax.plot([-0.5, starts[-1] - 0.5], [starts[i] - 0.5, starts[i] - 0.5],
                    color='black', linewidth=1)
            ax.plot([starts[i] - 0.5, starts[i] - 0.5], [-0.5, starts[-1] - 0.5],
                    color='black', linewidth=1)
    else:
        for i in range(1, starts.shape[0] - 1):
            ax.plot([-0.5, starts[-1, 1] - 0.5], [starts[i, 0] - 0.5, starts[i, 0] - 0.5],
                    color='black', linewidth=1)
            ax.plot([starts[i, 1] - 0.5, starts[i, 1] - 0.5], [-0.5, starts[-1, 0] - 0.5],
                    color='black', linewidth=1)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_fivec_dict(fname, array_dict, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                       regions=[], min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                       returnscale=False, **kwargs):
    """
    Fill in and bitmap from a 5C heatmap dictionary.

    :param fname: A file name to save plot to.
    :type fname: str.
    :param array_dict: Dictionary containing 5C data arrays and either 1 or two 2D mapping arrays (full or compact, respectively) with tuples are region name or pairs as keys.
    :type array_dict: dict.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting heatmap dict...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
    if len(regions) == 0:
        regions = []
        for key in array_dict.keys():
            regions.append(key[0])
            if len(key) > 1:
                regions.append(key[1])
        regions = list(set(regions))
        regions.sort()
    else:
        for i in range(len(regions))[::-1]:
            if (i,) not in array_dict:
                if not silent:
                    print("Region %i not in heatmap, skipping" % regions[i], file=sys.stderr)
                del regions[i]
    if len(regions) == 0:
        if not silent:
            print("No valid data found.", file=sys.stderr)
        return None
    if len(array_dict[(regions[0],)]) == 2:
        starts = [0]
        arraytype = 'full'
    else:
        starts = [[0, 0]]
        arraytype = 'compact'
    sizes = []
    for region in regions:
        if arraytype == 'full':
            sizes.append(array_dict[(region,)][1].shape[0])
            starts.append(starts[-1] + sizes[-1] + 1)
        else:
            sizes.append([array_dict[(region,)][1].shape[0],
                          array_dict[(region,)][2].shape[0]])
            starts.append([starts[-1][0] + sizes[-1][0] + 1, starts[-1][1] + sizes[-1][1] + 1])
    sizes = numpy.array(sizes, dtype=numpy.int32)
    starts = numpy.array(starts, dtype=numpy.int32)
    if len(sizes.shape) == 1:
        data = numpy.zeros((starts[-1] - 1, starts[-1] - 1, 2), dtype=numpy.float64)
    else:
        data = numpy.zeros((starts[-1, 0] - 1, starts[-1, 1] - 1, 2), dtype=numpy.float64)
    for i in range(len(regions)):
        if arraytype == 'full':
            data[starts[i]:(starts[i + 1] - 1),
                 starts[i]:(starts[i + 1]  - 1), :] = array_dict[(regions[i],)][0]
        else:
            data[starts[i, 0]:(starts[i + 1, 0] - 1),
                 starts[i, 1]:(starts[i + 1, 1]  - 1), :] = array_dict[(regions[i],)][0]
        for j in range(i + 1, len(regions)):
            key = (regions[i], regions[j])
            if key in array_dict:
                if arraytype == 'full':
                    data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), :] = array_dict[key]
                    data[starts[j]:(starts[j + 1] - 1),
                         starts[i]:(starts[i + 1] - 1), :] = array_dict[key].transpose(1, 0, 2)
                else:
                    data[starts[i, 0]:(starts[i + 1, 0] - 1), starts[j, 1]:(starts[j + 1, 1] - 1), :] = array_dict[key]
                    data[starts[j, 0]:(starts[j + 1, 0] - 1),
                         starts[i, 1]:(starts[i + 1, 1] - 1), :] = array_dict[(regions[j], regions[i])]
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
    img = numpy.full((data.shape[0], data.shape[1]), numpy.nan, dtype=numpy.float32)
    where = numpy.where(data[:, :, 1] > 0)
    img[where] = data[where[0], where[1], 0]
    fig, ax = plt.subplots(1, 1, figsize=(img.shape[0]/100, img.shape[1]/100))
    ax.imshow(img.T, norm=plt.Normalize(-1, 1), cmap=cmap)
    if len(starts.shape) == 1:
        for i in range(1, starts.shape[0] - 1):
            ax.plot([-0.5, starts[-1] - 0.5], [starts[i] - 0.5, starts[i] - 0.5],
                    color='black', linewidth=1)
            ax.plot([starts[i] - 0.5, starts[i] - 0.5], [-0.5, starts[-1] - 0.5],
                    color='black', linewidth=1)
    else:
        for i in range(1, starts.shape[0] - 1):
            ax.plot([-0.5, starts[-1, 1] - 0.5], [starts[i, 0] - 0.5, starts[i, 0] - 0.5],
                    color='black', linewidth=1)
            ax.plot([starts[i, 1] - 0.5, starts[i, 1] - 0.5], [-0.5, starts[-1, 0] - 0.5],
                    color='black', linewidth=1)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_fivec_compact_heatmap_dict(fname, hdf_fname, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                    regions=[], min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                                    returnscale=False, **kwargs):
    """
    Fill in and rescale bitmap in a compact format from a 5C heatmap h5dict file.

    This plots the data in a 5C compact format such that the rows correspond to positive-strand primers and columns correspond to negative-strand primers.

    :param fname: A file name to save plot to.
    :type fname: str.
    :param hdf_fname: Location of a heatmap h5dict containing 5C data arrays.
    :type hdf_fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting heatmap dict...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
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
    img = numpy.full((data.shape[0], data.shape[1]), numpy.nan, dtype=numpy.float32)
    where = numpy.where(data[:, :, 1] > 0)
    img[where] = data[where[0], where[1], 0]
    fig, ax = plt.subplots(1, 1, figsize=(img.shape[0]/100, img.shape[1]/100))
    ax.imshow(img.T, norm=plt.Normalize(-1, 1), cmap=cmap)
    if len(starts.shape) == 1:
        for i in range(1, starts.shape[0] - 1):
            ax.plot([-0.5, starts[-1] - 0.5], [starts[i] - 0.5, starts[i] - 0.5],
                    color='black', linewidth=1)
            ax.plot([starts[i] - 0.5, starts[i] - 0.5], [-0.5, starts[-1] - 0.5],
                    color='black', linewidth=1)
    else:
        for i in range(1, starts.shape[0] - 1):
            ax.plot([-0.5, starts[-1, 1] - 0.5], [starts[i, 0] - 0.5, starts[i, 0] - 0.5],
                    color='black', linewidth=1)
            ax.plot([starts[i, 1] - 0.5, starts[i, 1] - 0.5], [-0.5, starts[-1, 0] - 0.5],
                    color='black', linewidth=1)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_diagonal_from_upper_array(fname, data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                   min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                                   returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a HiC upper array, plotting only the upper triangle rotated 45 degrees counter-clockwise.

    :param fname: A file name to save plot to.
    :type fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting rotated compact array...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
    dim = int(0.5 + (0.25 + data.shape[0] * 2) ** 0.5) - diag
    scaled = numpy.zeros((dim, dim, 2), dtype=numpy.float64)
    triu = numpy.triu_indices(dim, 1 - diag)
    scaled[triu[0], triu[1], :] = data
    where = numpy.where(scaled[:, :, 1] > 0)
    scaled[where[0], where[1], 0] /= scaled[where[0], where[1], 1]
    if logged:
        where = numpy.where(scaled[:, :, 0] <= 0)
        scalde[where[0], where[1], 1] = 0
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
    img = numpy.full((scaled.shape[0], scaled.shape[1]), numpy.nan, dtype=numpy.float32)
    where = numpy.where(scaled[:, :, 1] > 0)
    img[where] = scaled[where[0], where[1], 0]
    fig, ax = plt.subplots(1, 1, figsize=(img.shape[0]/100 * 2**0.5, img.shape[0]/200 * 2**0.5))
    ax.imshow(img.T, norm=plt.Normalize(-1, 1), cmap=cmap, origin='upper',
              transform=transforms.Affine2D().rotate_deg(45))
    ax.set_ylim(-img.shape[0] * 2**0.5, 0)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_diagonal_from_compact_array(fname, data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                     min_color="0000ff", mid_color="ffffff", max_color="ff0000",
                                     returnscale=False, diagonal_included=False, **kwargs):
    """
    Fill in and rescale bitmap from a HiC compact array, plotting only the upper triangle rotated 45 degrees counter-clockwise.

    :param fname: A file name to save plot to.
    :type fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    if not silent:
        print("Plotting rotated compact array...", end='', file=sys.stderr)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
    scaled = numpy.zeros((data.shape[0] * 2 - 1 + diag * 2, data.shape[1] + 1, 2), dtype=numpy.float64)
    xdim, ydim = scaled.shape[:2]
    for i in range(data.shape[1]):
        scaled[(i + 1):(2 * data.shape[0] - i - 1 + diag * 2):2, ydim - i - 1, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :])
        scaled[i:(2 * data.shape[0] - i - 2 + diag * 2):2, ydim - i - 1, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
        scaled[(i + 2):(2 * data.shape[0] - i + diag * 2):2, ydim - i - 1, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
        scaled[(i + 1):(2 * data.shape[0] - i - 1 + diag * 2):2, ydim - i - 2, :] += (
            data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
        if i > 0:
            scaled[(i + 1):(2 * data.shape[0] - i - 1 + diag * 2):2, ydim - i, :] += (
                data[:(data.shape[0] - i - 1 + diag), i, :]) * 0.25
    where = numpy.where(scaled[:, :, 1] > 0)
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
    img = numpy.full((scaled.shape[0], scaled.shape[1]), numpy.nan, dtype=numpy.float32)
    where = numpy.where(scaled[:, :, 1] > 0)
    img[where] = scaled[where[0], where[1], 0]
    fig, ax = plt.subplots(1, 1, figsize=(img.shape[0]/100 * 2**0.5, img.shape[0]/200 * 2**0.5))
    ax.imshow(img.T, norm=plt.Normalize(-1, 1), cmap=cmap, origin='upper',
              transform=transforms.Affine2D().rotate_deg(45))
    ax.set_ylim(-img.shape[0] * 2**0.5, 0)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return


def plot_key(fname, min_score, max_score, height, width, labelformat='0.2f', orientation='left', num_ticks=5,
             min_color="0000ff", mid_color="ffffff", max_color="ff0000", labelattr=None, log_display=True, **kwargs):
    """
    Create a key including color gradient and labels indicating associated values, returning a :mod:`pyx` canvas.

    :param fname: A file name to save plot to.
    :type fname: str.
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
    if 'matplotlib' not in sys.modules.keys():
        if not silent:
            print("Matplotlib must be installed to use this function.", file=sys.stderr)
        return None
    height = float(height)
    width = float(width)
    min_score = float(min_score)
    max_score = float(max_score)
    if mid_color is None:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list([f"#{min_color}", f"#{max_color}"])
    else:
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
            [f"#{min_color}", f"#{mid_color}", f"#{max_color}"])
    img = numpy.linspace(minscore, maxscore, 512).reshape(-1, 1)
    fig, ax = plt.subplots(1, 1, figsize=(width, height))
    if orientation in ['right', 'left']:
        img = img.T
    ticks = numpy.linspace(minscore, maxscore, num_ticks)
    labels = [f"{num:{format}}".format(num=x, format=labelformat) for x in ticks]
    ax.imshow(img, aspect='auto', norm=plt.Normalize(0, 511), cmap=cmap)
    if orientation == 'left':
        ax.get_xaxis().set_visible(False)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
    elif orientation == 'bottom':
        ax.get_yaxis().set_visible(False)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    elif orientation == 'right':
        ax.get_xaxis().set_visible(False)
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
    else:
        ax.get_yaxis().set_visible(False)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    fig.tight_layout()
    plt.savefig(fname)
    plt.close()
    if not silent:
        print("Done", file=sys.stderr)
    if returnscale:
        return [minscore, maxscore]
    else:
        return

