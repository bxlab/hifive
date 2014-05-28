#!/usr/bin/env python
#(c) 2014 Michael Sauria (mike.sauria@gmail.com)

"""
This is a module contains scripts for generating plots from compact and full
matrices of interaction data.

Input data
----------

These functions take either compact or full 3d data matrices.

Concepts
--------

Data can either be arranged in compact or complete arrays. Compact arrays
are N x M, where N is the number of fends or bins, and M is the maximum
distance between fends or bins. This is useful for working with sets of short
interactions. Data can be raw, fend-corrected, distance-dependence removed,
or enrichment values. Arrays are 3-dimensional with observed values in the
first layer of d3, expected values in the second layer of d3.

-----------------------------------------------------------------------------

API documentation
-----------------



"""

import os
import sys
from math import log, exp

import h5py
import numpy
try:
    from PIL import Image
except:
    pass
try:
    import pyx
except:
    pass


def plot_compact_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                       min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0), max_color=(1.0, 0.0, 0.0)):
    """
    plot_compact_array method

    Rescale and fill in bitmap from a HiC compact array.

    Parameters
    ----------
    data : 3d numpy array
        A compact array of interaction data.
    maxscore : float, optional
        A ceiling value for plotting.
    minscore : float, optional
        A floot value for plotting.
    symmetricscaling : bool, optional
        Indicates whether to recenter data for scaling or maintain scores about zero.
    logged : bool, optional
        Indicates whether to use log values of scores.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle plot value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    """
    if 'PIL' not in sys.modules.keys():
        print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    print >> sys.stderr, ("Plotting compact array..."),
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * max_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * max_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * max_color[2], 256))
                 ).astype(numpy.uint32)
    dim = data.shape[0]
    scaled = numpy.copy(data)
    where = numpy.where(data[:, :, 1] > 0)
    scaled[where[0], where[1], 0] /= scaled[where[0], where[1], 1]
    if logged:
        where = numpy.where(scaled[:, :, 0] <= 0)
        scaled[where[0], where[1], 1] = 0
        where = numpy.where(scaled[:, :, 1] > 0)
        scaled[where[0], where[1], 0] = numpy.log(scaled[where[0], where[1], 0])
    scaled[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(scaled[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(scaled[where[0], where[1], 0])
    if symmetricscaling:
        scaled[where[0], where[1], 0] /= max(abs(maxscore), abs(minscore))
    else:
        scaled[where[0], where[1], 0] -= minscore
        scaled[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        scaled[where[0], where[1], 0] -= 1.0
    scaled = numpy.minimum(1.0, numpy.maximum(-1.0, scaled))
    scaled[where[0], where[1], 0] *= 255
    scaled = numpy.round(scaled).astype(numpy.int32)
    img = numpy.empty((dim, dim), dtype=numpy.uint32)
    img[:, :] = int('ff999999', 16)
    for i in range(dim - 1):
        where = numpy.where((scaled[i, :min(scaled.shape[1], dim - i - 1), 1] > 0) *
                            (scaled[i, :min(scaled.shape[1], dim - i - 1), 0] >= 0))[0]
        img[i, where + i + 1] = gradient2[scaled[i, where, 0]]
        img[where + i + 1, i] = img[i, where + i + 1]
        where = numpy.where((scaled[i, :min(scaled.shape[1], dim - i - 1), 1] > 0) *
                            (scaled[i, :min(scaled.shape[1], dim - i - 1), 0] < 0))[0]
        img[i, where + i + 1] = gradient1[-scaled[i, where, 0]]
        img[where + i + 1, i] = img[i, where + i + 1]
    pilImage = Image.frombuffer('RGBA', (dim, dim), img, 'raw', 'RGBA', 0, 1)
    print >> sys.stderr, ("Done\n"),
    return pilImage


def plot_full_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                    min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0), max_color=(1.0, 0.0, 0.0)):
    """
    plot_full_array method

    Rescale and fill in bitmap from a full array.

    Parameters
    ----------
    data : 3d numpy array
        A full array of interaction data.
    maxscore : float, optional
        A ceiling value for plotting.
    minscore : float, optional
        A floot value for plotting.
    symmetricscaling : bool, optional
        Indicates whether to recenter data for scaling or maintain scores about zero.
    logged : bool, optional
        Indicates whether to use log values of scores.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle plot value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    """
    if 'PIL' not in sys.modules.keys():
        print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    print >> sys.stderr, ("Plotting full array..."),
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * max_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * max_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * max_color[2], 256))
                 ).astype(numpy.uint32)
    xdim, ydim = data.shape[:2]
    scaled = numpy.copy(data)
    where = numpy.where(data[:, :, 1] > 0)
    scaled[where[0], where[1], 0] /= scaled[where[0], where[1], 1]
    if logged:
        where = numpy.where(scaled[:, :, 0] <= 0)
        scaled[where[0], where[1], 1] = 0
        where = numpy.where(scaled[:, :, 1] > 0)
        scaled[where[0], where[1], 0] = numpy.log(scaled[where[0], where[1], 0])
    scaled[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(scaled[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(scaled[where[0], where[1], 0])
    if symmetricscaling:
        scaled[where[0], where[1], 0] /= max(abs(maxscore), abs(minscore))
    else:
        scaled[where[0], where[1], 0] -= minscore
        scaled[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        scaled[where[0], where[1], 0] -= 1.0
    scaled = numpy.minimum(1.0, numpy.maximum(-1.0, scaled))
    scaled[where[0], where[1], 0] *= 255
    scaled = numpy.round(scaled).astype(numpy.int32)
    img = numpy.empty((xdim, ydim), dtype=numpy.uint32)
    img.shape = (ydim, xdim)
    img[:, :] = int('ff999999', 16)
    where = numpy.where((scaled[:, :, 1] > 0) * (scaled[:, :, 0] >= 0))
    img[where[1], where[0]] = gradient2[scaled[where[0], where[1], 0]]
    where = numpy.where((scaled[:, :, 1] > 0) * (scaled[:, :, 0] < 0))
    img[where[1], where[0]] = gradient1[-scaled[where[0], where[1], 0]]
    pilImage = Image.frombuffer('RGBA', (xdim, ydim), img, 'raw', 'RGBA', 0, 1)
    print >> sys.stderr, ("Done\n"),
    return pilImage


def plot_upper_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                     min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0), max_color=(1.0, 0.0, 0.0)):
    """
    plot_upper_array method

    Rescale and fill in bitmap from a full array.

    Parameters
    ----------
    data : 2d numpy array
        A flattened upper-triangle array of interaction data.
    maxscore : float, optional
        A ceiling value for plotting.
    minscore : float, optional
        A floot value for plotting.
    symmetricscaling : bool, optional
        Indicates whether to recenter data for scaling or maintain scores about zero.
    logged : bool, optional
        Indicates whether to use log values of scores.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle plot value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    """
    if 'PIL' not in sys.modules.keys():
        print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    print >> sys.stderr, ("Plotting upper array..."),
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * max_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * max_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * max_color[2], 256))
                 ).astype(numpy.uint32)
    dim = int(0.5 + (0.25 + 2 * data.shape[0]) ** 0.5)
    scaled = numpy.copy(data)
    where = numpy.where(data[:, 1] > 0)[0]
    scaled[where, 0] /= scaled[where, 1]
    if logged:
        where = numpy.where(scaled[:, 0] <= 0)[0]
        scaled[where, 1] = 0
        where = numpy.where(scaled[:, 1] > 0)[0]
        scaled[where, 0] = numpy.log(scaled[where, 0])
    scaled[where, 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(scaled[where, 0])
    if minscore is None:
        minscore = numpy.amin(scaled[where, 0])
    if symmetricscaling:
        scaled[where, 0] /= max(abs(maxscore), abs(minscore))
    else:
        scaled[where, 0] -= minscore
        scaled[where, 0] /= (maxscore - minscore) * 0.5
        scaled[where, 0] -= 1.0
    scaled = numpy.minimum(1.0, numpy.maximum(-1.0, scaled))
    scaled[where, 0] *= 255
    scaled = numpy.round(scaled).astype(numpy.int32)
    img = numpy.empty((dim, dim), dtype=numpy.uint32)
    img[:, :] = int('ff999999', 16)
    for i in range(dim - 1):
        index = i * dim - i * (i + 1) / 2
        where = numpy.where((scaled[index:(index + dim - i - 1), 1] > 0) *
                            (scaled[index:(index + dim - i - 1), 0] >= 0))[0]
        img[i, where + i + 1] = gradient2[scaled[where + index, 0]]
        img[where + i + 1, i] = img[i, where + i + 1]
        where = numpy.where((scaled[index:(index + dim - i - 1), 1] > 0) *
                            (scaled[index:(index + dim - i - 1), 0] < 0))[0]
        img[i, where + i + 1] = gradient1[-scaled[where + index, 0]]
        img[where + i + 1, i] = img[i, where + i + 1]
    pilImage = Image.frombuffer('RGBA', (dim, dim), img, 'raw', 'RGBA', 0, 1)
    print >> sys.stderr, ("Done\n"),
    return pilImage


def plot_hic_heatmap_dict(filename, maxscore=None, minscore=None, symmetricscaling=True, logged=True, chroms=[],
                          min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0), max_color=(1.0, 0.0, 0.0)):
    """
    plot_hic_heatmap_dict method

    Rescale and fill in bitmap from a binned array h5dict.

    Parameters
    ----------
    filename : str
        Location of h5dict containing binned data arrays.
    maxscore : float, optional
        A ceiling value for plotting.
    minscore : float, optional
        A floot value for plotting.
    symmetricscaling : bool, optional
        Indicates whether to recenter data for scaling or maintain scores about zero.
    logged : bool, optional
        Indicates whether to use log values of scores.
    chroms : list, optional
        If specified, only the indicated chromosomes are plotted. Otherwise all chromosomes present in the h5dict are
        plotted.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle plot value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    """
    if 'PIL' not in sys.modules.keys():
        print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    print >> sys.stderr, ("Plotting heatmap dict..."),
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * max_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * max_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * max_color[2], 256))
                 ).astype(numpy.uint32)
    input = h5py.File(filename, 'r')
    if len(chroms) == 0:
        chroms = list(input['chromosomes'][...])
    starts = [0]
    sizes = []
    for chrom in chroms:
        sizes.append(input['%s.positions' % chrom].shape[0])
        starts.append(starts[-1] + sizes[-1] + 1)
    data = numpy.zeros((starts[-1] - 1, starts[-1] - 1, 2), dtype=numpy.float32)
    for i in range(len(chroms)):
        indices = numpy.triu_indices(sizes[i], 1)
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
        data[where[0], where[1], 0] = numpy.log(data[where[0], where[1], 0])
    data[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(data[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(data[where[0], where[1], 0])
    if symmetricscaling:
        data[where[0], where[1], 0] /= max(abs(maxscore), abs(minscore))
    else:
        data[where[0], where[1], 0] -= minscore
        data[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        data[where[0], where[1], 0] -= 1.0
    data = numpy.minimum(1.0, numpy.maximum(-1.0, data))
    data[where[0], where[1], 0] *= 255
    data = numpy.round(data).astype(numpy.int32)
    img = numpy.empty((data.shape[0], data.shape[0]), dtype=numpy.uint32)
    img[:, :] = int('ff999999', 16)
    where = numpy.where((data[:, :, 1] > 0) * (data[:, :, 0] >= 0))
    img[where] = gradient2[data[where[0], where[1], 0]]
    img[where[1], where[0]] = img[where]
    where = numpy.where((data[:, :, 1] > 0) * (data[:, :, 0] < 0))
    img[where] = gradient1[-data[where[0], where[1], 0]]
    img[where[1], where[0]] = img[where]
    pilImage = Image.frombuffer('RGBA', (data.shape[0], data.shape[0]), img, 'raw', 'RGBA', 0, 1)
    print >> sys.stderr, ("Done\n"),
    return pilImage


def plot_fivec_full_heatmap_dict(filename, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                 regions=[], min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0),
                                 max_color=(1.0, 0.0, 0.0)):
    """
    plot_fivec_full_heatmap_dict method

    Rescale and fill in bitmap from a binned array h5dict.

    Parameters
    ----------
    filename : str
        Location of h5dict containing binned data arrays.
    maxscore : float, optional
        A ceiling value for plotting.
    minscore : float, optional
        A floot value for plotting.
    symmetricscaling : bool, optional
        Indicates whether to recenter data for scaling or maintain scores about zero.
    logged : bool, optional
        Indicates whether to use log values of scores.
    regions : list, optional
        If specified, only the indicated regions are plotted. Otherwise all regions present in the h5dict are
        plotted.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle plot value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    """
    if 'PIL' not in sys.modules.keys():
        print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    print >> sys.stderr, ("Plotting heatmap dict..."),
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * max_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * max_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * max_color[2], 256))
                 ).astype(numpy.uint32)
    input = h5py.File(filename, 'r')
    if len(regions) == 0:
        regions = []
        for key in input.keys():
            if key.count('positions') > 0 or key.count('fragments') > 0:
                regions.append(int(key.split('.')[0]))
        regions.sort()
    starts = [0]
    sizes = []
    for region in regions:
        if "%i.positions" % region in input:
            sizes.append(input['%i.positions' % region].shape[0])
        else:
            sizes.append(input['%i.fragments' % region].shape[0])
        starts.append(starts[-1] + sizes[-1] + 1)
    data = numpy.zeros((starts[-1] - 1, starts[-1] - 1, 2), dtype=numpy.float32)
    for i in range(len(regions)):
        data[starts[i]:(starts[i + 1] - 1), starts[i]:(starts[i + 1]  - 1), 0] = input['%i.counts' % regions[i]][...]
        data[starts[i]:(starts[i + 1] - 1), starts[i]:(starts[i + 1]  - 1), 1] = input['%i.expected' % regions[i]][...]
    for i in range(len(regions) - 1):
        for j in range(i + 1, len(regions)):
            if '%i_by_%i.counts' % (regions[i], regions[j]) in input.keys():
                data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), 0] = (
                    input['%i_by_%i.counts' % (regions[i], regions[j])][:, :])
                data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), 1] = (
                    input['%i_by_%i.expected' % (regions[i], regions[j])][:, :])
                data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), :] = numpy.transpose(
                    data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), :], axes=[1, 0, 2])
            elif '%i_by_%i.counts' % (regions[j], regions[i]) in input.keys():
                data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), 0] = (
                    input['%i_by_%i.counts' % (regions[j], regions[i])][:, :])
                data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), 1] = (
                    input['%i_by_%i.expected' % (regions[j], regions[i])][:, :])
                data[starts[i]:(starts[i + 1] - 1), starts[j]:(starts[j + 1] - 1), :] = numpy.transpose(
                    data[starts[j]:(starts[j + 1] - 1), starts[i]:(starts[i + 1] - 1), :], axes=[1, 0, 2])
    where = numpy.where(data[:, :, 1] > 0)
    data[where[0], where[1], 0] /= data[where[0], where[1], 1]
    if logged:
        where = numpy.where(data[:, :, 0] <= 0)
        data[where[0], where[1], 1] = 0
        where = numpy.where(data[:, :, 1] > 0)
        data[where[0], where[1], 0] = numpy.log(data[where[0], where[1], 0])
    data[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(data[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(data[where[0], where[1], 0])
    if symmetricscaling:
        data[where[0], where[1], 0] /= max(abs(maxscore), abs(minscore))
    else:
        data[where[0], where[1], 0] -= minscore
        data[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        data[where[0], where[1], 0] -= 1.0
    data = numpy.minimum(1.0, numpy.maximum(-1.0, data))
    data[where[0], where[1], 0] *= 255
    data = numpy.round(data).astype(numpy.int32)
    img = numpy.empty((data.shape[0], data.shape[0]), dtype=numpy.uint32)
    img[:, :] = int('ff999999', 16)
    where = numpy.where((data[:, :, 1] > 0) * (data[:, :, 0] >= 0))
    img[where] = gradient2[data[where[0], where[1], 0]]
    img[where[1], where[0]] = img[where]
    where = numpy.where((data[:, :, 1] > 0) * (data[:, :, 0] < 0))
    img[where] = gradient1[-data[where[0], where[1], 0]]
    img[where[1], where[0]] = img[where]
    black = int('ff000000', 16)
    for i in range(1, len(starts) - 1):
        img[starts[i] - 1, :] = black
        img[:, starts[i] - 1] = black
    pilImage = Image.frombuffer('RGBA', (data.shape[0], data.shape[0]), img, 'raw', 'RGBA', 0, 1)
    print >> sys.stderr, ("Done\n"),
    return pilImage


def plot_fivec_compact_heatmap_dict(filename, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                    regions=[], min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0),
                                    max_color=(1.0, 0.0, 0.0)):
    """
    plot_fivec_compact_heatmap_dict method

    Rescale and fill in bitmap from a binned array h5dict.

    Parameters
    ----------
    filename : str
        Location of h5dict containing binned data arrays.
    maxscore : float, optional
        A ceiling value for plotting.
    minscore : float, optional
        A floot value for plotting.
    symmetricscaling : bool, optional
        Indicates whether to recenter data for scaling or maintain scores about zero.
    logged : bool, optional
        Indicates whether to use log values of scores.
    regions : list, optional
        If specified, only the indicated regions are plotted. Otherwise all regions present in the h5dict are
        plotted.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle plot value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    """
    if 'PIL' not in sys.modules.keys():
        print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    print >> sys.stderr, ("Plotting heatmap dict..."),
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * max_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * max_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * max_color[2], 256))
                 ).astype(numpy.uint32)
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
    data = numpy.zeros((xstarts[-1] - 1, ystarts[-1] - 1, 2), dtype=numpy.float32)
    for i in range(len(regions)):
        data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[i]:(ystarts[i + 1]  - 1), 0] = (
                input['%i.counts' % regions[i]][...])
        data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[i]:(ystarts[i + 1]  - 1), 1] = (
                input['%i.expected' % regions[i]][...])
    for i in range(len(regions) - 1):
        for j in range(i + 1, len(regions)):
            data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[j]:(ystarts[j + 1] - 1), 0] = (
                input['%i_by_%i.counts' % (regions[i], regions[j])][:, :])
            data[xstarts[i]:(xstarts[i + 1] - 1), ystarts[j]:(ystarts[j + 1] - 1), 1] = (
                input['%i_by_%i.expected' % (regions[i], regions[j])][:, :])
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
        data[where[0], where[1], 0] = numpy.log(data[where[0], where[1], 0])
    data[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(data[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(data[where[0], where[1], 0])
    if symmetricscaling:
        data[where[0], where[1], 0] /= max(abs(maxscore), abs(minscore))
    else:
        data[where[0], where[1], 0] -= minscore
        data[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        data[where[0], where[1], 0] -= 1.0
    data = numpy.minimum(1.0, numpy.maximum(-1.0, data))
    data[where[0], where[1], 0] *= 255
    data = numpy.round(data).astype(numpy.int32)
    img = numpy.empty((data.shape[0], data.shape[1]), dtype=numpy.uint32)
    img.shape = (data.shape[1], data.shape[0])
    img[:, :] = int('ff999999', 16)
    where = numpy.where((data[:, :, 1] > 0) * (data[:, :, 0] >= 0))
    img[where[1], where[0]] = gradient2[data[where[0], where[1], 0]]
    where = numpy.where((data[:, :, 1] > 0) * (data[:, :, 0] < 0))
    img[where[1], where[0]] = gradient1[-data[where[0], where[1], 0]]
    black = int('ff000000', 16)
    for i in range(1, len(xstarts) - 1):
        img[:, xstarts[i] - 1] = black
    for i in range(1, len(ystarts) - 1):
        img[ystarts[i] - 1, :] = black
    pilImage = Image.frombuffer('RGBA', (data.shape[0], data.shape[1]), img, 'raw', 'RGBA', 0, 1)
    print >> sys.stderr, ("Done\n"),
    return pilImage


def plot_diagonal_from_compact_array(data, maxscore=None, minscore=None, symmetricscaling=True, logged=True,
                                     min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0), max_color=(1.0, 0.0, 0.0)):
    """
    plot_diagonal_from_compact_array method

    Rotate 45 degrees, rescale, and fill in bitmap from a HiC compact array.

    Parameters
    ----------
    data : 3d numpy array
        A compact array of interaction data.
    maxscore : float, optional
        A ceiling value for plotting.
    minscore : float, optional
        A floot value for plotting.
    symmetricscaling : bool, optional
        Indicates whether to recenter data for scaling or maintain scores about zero.
    logged : bool, optional
        Indicates whether to use log values of scores.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle plot value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum plot value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    """
    if 'PIL' not in sys.modules.keys():
        print >> sys.stderr, ("The PIL module must be installed to use this function.")
        return None
    print >> sys.stderr, ("Plotting rotated compact array..."),
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256 ** 2 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * max_color[0], 256)) +
                 256 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * max_color[1], 256)) +
                 256 ** 2 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * max_color[2], 256))
                 ).astype(numpy.uint32)
    rotated = numpy.zeros((data.shape[0] * 2 - 1, data.shape[1] + 2, 2), dtype=numpy.float32)
    xdim, ydim = rotated.shape[:2]
    for i in range(data.shape[1]):
        rotated[(i + 1):(2 * data.shape[0] - i - 1):2, ydim - i - 2, :] += data[:(data.shape[0] - i - 1), i, :]
        rotated[i:(2 * data.shape[0] - i - 2):2, ydim - i - 2, :] += data[:(data.shape[0] - i - 1), i, :] * 0.25
        rotated[(i + 2):(2 * data.shape[0] - i):2, ydim - i - 2, :] += data[:(data.shape[0] - i - 1), i, :] * 0.25
        rotated[(i + 1):(2 * data.shape[0] - i - 1):2, ydim - i - 1, :] += data[:(data.shape[0] - i - 1), i, :] * 0.25
        rotated[(i + 1):(2 * data.shape[0] - i - 1):2, ydim - i - 3, :] += data[:(data.shape[0] - i - 1), i, :] * 0.25
    where = numpy.where(rotated[:, :, 1] > 0)
    rotated[where[0], where[1], 0] /= rotated[where[0], where[1], 1]
    if logged:
        where = numpy.where(rotated[:, :, 0] <= 0)
        rotated[where[0], where[1], 1] = 0
        where = numpy.where(rotated[:, :, 1] > 0)
        rotated[where[0], where[1], 0] = numpy.log(rotated[where[0], where[1], 0])
    rotated[where[0], where[1], 1] = 1
    if maxscore is None:
        maxscore = numpy.amax(rotated[where[0], where[1], 0])
    if minscore is None:
        minscore = numpy.amin(rotated[where[0], where[1], 0])
    if symmetricscaling:
        rotated[where[0], where[1], 0] /= max(abs(maxscore), abs(minscore))
    else:
        rotated[where[0], where[1], 0] -= minscore
        rotated[where[0], where[1], 0] /= (maxscore - minscore) * 0.5
        rotated[where[0], where[1], 0] -= 1.0
    rotated = numpy.minimum(1.0, numpy.maximum(-1.0, rotated))
    rotated[where[0], where[1], 0] *= 255
    rotated = numpy.round(rotated).astype(numpy.int32)
    img = numpy.empty((xdim, ydim), dtype=numpy.uint32)
    img.shape = (ydim, xdim)
    img[:, :] = int('ff999999', 16)
    where = numpy.where((rotated[:, :, 0] >= 0) * (rotated[:, :, 1] > 0))
    img[where[1], where[0]] = gradient2[rotated[where[0], where[1], 0]]
    where = numpy.where((rotated[:, :, 0] < 0) * (rotated[:, :, 1] > 0))
    img[where[1], where[0]] = gradient1[-rotated[where[0], where[1], 0]]
    pilImage = Image.frombuffer('RGBA', (xdim, ydim), img, 'raw', 'RGBA', 0, 1)
    print >> sys.stderr, ("Done\n"),
    return pilImage


def plot_key(min_score, max_score, height, width, labelformat='%0.2f', orientation='left', num_ticks=5,
             min_color=(0.0, 0.0, 1.0), mid_color=(1.0, 1.0, 1.0), max_color=(1.0, 0.0, 0.0), labelattr=None,
             log_display=True):
    """
    plot_key method

    Create a key including color gradient and labels indicating asscociated values, returning a pyx canvas.

    Parameters
    ----------
    min_score : float
        The minimum value of the key scale.
    max_score : float
        The maximum value of the key scale.
    height : float
        The height of the gradient bar.
    width : float
        The width of the gradient bar.
    labelformat : string
        A string denoting the format for displaying number labels using the string formatting style from Python <= 2.6.
    orientation : string
        Indicates where labels are placed relative to gradient bar. This parameter will accept 'left', 'right', 'top',
        and 'bottom'.
    num_ticks : int
        Indicates how many evenly-spaced tick marks and associated labels to insert. This can be zero for no labels or
        greater than one. Labels are inserted at the minimum and maximum values first with remaining ticks occuring
        evenly distributed between the extremes.
    min_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the minimum key value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with max_color and possibly mid_color.
    mid_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the middle key value, respectively. Numbers range from 0 to 1. This can be set to None to
        create a gradient ranging from min_color to max_color or to a tuple to create a divergent gradient.
    max_color : tuple
        This is a tuple containing three numbers representing the red, green, and blue component of the color
        associated with the maximum key value, respectively. Numbers range from 0 to 1. This variable is used to
        create a color gradient for plotting along with min_color and possibly mid_color.
    labelattr : list
        A list of pyx attributes to be passed to the text function.
    log_display : bool
        If True, min_score and max_score are taken to be logged values and so labels are evenly spaced in log space
        but converted to normal space for display.
    """
    if 'PIL' not in sys.modules.keys() or 'pyx' not in sys.modules.keys():
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
    if mid_color is None:
        mid_color = ((min_color[0] + max_color[0]) / 2.0, (min_color[1] + max_color[1]) / 2.0,
                     (min_color[2] + max_color[2]) / 2.0)
    gradient1 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * mid_color[0], 255.0 * min_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * mid_color[1], 255.0 * min_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * mid_color[2], 255.0 * min_color[2], 256))
                 ).astype(numpy.uint32)
    gradient2 = (256**3 * 255 + numpy.round(numpy.linspace(255.0 * max_color[0], 255.0 * mid_color[0], 256)) +
                 256.0 * numpy.round(numpy.linspace(255.0 * max_color[1], 255.0 * mid_color[1], 256)) +
                 256.0 ** 2.0 * numpy.round(numpy.linspace(255.0 * max_color[2], 255.0 * mid_color[2], 256))
                 ).astype(numpy.uint32)
    img[:255, :] = gradient2[:255].reshape(-1, 1)
    img[255:, :] = gradient1.reshape(-1, 1)
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
