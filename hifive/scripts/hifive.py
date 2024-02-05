#!/usr/bin/env python

"""Description: HiFive main executable.

Copyright (c) 2015 Michael Sauria <mike.sauria@jhu.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

import sys
import argparse as ap

try:
    from mpi4py import MPI
except:
    pass

from hifive.version import __version__ as VERSION

def main():
    global comm, rank, num_procs
    if "mpi4py" in sys.modules.keys():
        try:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()
        except:
            comm = None
            rank = 0
            num_procs = 1
    else:
        comm = None
        rank = 0
        num_procs = 1

    """The main function/pipeline for HiFive."""
    parser = generate_parser()
    args = None
    if not comm is None:
        try:
            if rank == 0:
                args = parser.parse_args()
        finally:
            args = comm.bcast(args, root=0)
    else:
        args = parser.parse_args()

    subcommand = args.subcommand

    if subcommand == "connect" and rank == 0:
        from hifive.commands.connect_files import run
        run(args)
    if subcommand == "fragments" and rank == 0:
        from hifive.commands.create_fragmentset import run
        run(args)
    elif subcommand == "5c-data" and rank == 0:
        from hifive.commands.create_fivec_dataset import run
        run(args)
    elif subcommand == "5c-project" and rank == 0:
        from hifive.commands.create_fivec_project import run
        run(args)
    elif subcommand == "5c-normalize" and rank == 0:
        from hifive.commands.normalize_fivec_project import run
        run(args)
    elif subcommand == "5c-complete":
        from hifive.commands.complete_fivec_project import run
        run(args)
    elif subcommand == "5c-heatmap":
        from hifive.commands.create_fivec_heatmap import run
        run(args)
    elif subcommand == "5c-interval":
        from hifive.commands.get_fivec_interval import run
        run(args)
    elif subcommand == "5c-combine-replicates" and rank == 0:
        from hifive.commands.combine_fivec_replicates import run
        run(args)
    elif subcommand == "fends" and rank == 0:
        from hifive.commands.create_fendset import run
        run(args)
    elif subcommand == "hic-data" and rank == 0:
        from hifive.commands.create_hic_dataset import run
        run(args)
    elif subcommand == "hic-project":
        from hifive.commands.create_hic_project import run
        run(args)
    elif subcommand == "hic-normalize":
        from hifive.commands.normalize_hic_project import run
        run(args)
    elif subcommand == "hic-complete":
        from hifive.commands.complete_hic_project import run
        run(args)
    elif subcommand == "hic-heatmap":
        from hifive.commands.create_hic_heatmap import run
        run(args)
    elif subcommand == "hic-mrheatmap":
        from hifive.commands.create_hic_mrheatmap import run
        run(args)
    elif subcommand == "hic-interval":
        from hifive.commands.get_hic_interval import run
        run(args)
    elif subcommand == "hic-combine-replicates" and rank == 0:
        from hifive.commands.combine_hic_replicates import run
        run(args)
    elif subcommand == "quasar":
        from hifive.commands.find_quasar_scores import run
        run(args)

def generate_parser():
    """Generate an argument parser."""
    description = "%(prog)s -- Data handling, normalization, manipulation, and plotting for HiC and 5C experimental data"
    epilog = "For command line options of each command, type: %(prog)s <COMMAND> -h"
    parser = ap.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("--version", action="version", version="%(prog)s %(version_num)s" % {'prog':parser.prog, 'version_num':VERSION})
    subparsers = parser.add_subparsers(dest='subcommand', required=True)

    add_connect_subparser(subparsers)
    add_fragments_subparser(subparsers)
    add_fivecdataset_subparser(subparsers)
    add_fivecproject_subparser(subparsers)
    add_fivecnormalize_subparser(subparsers)
    add_complete_fivec_subparser(subparsers)
    add_fivec_heatmap_subparser(subparsers)
    add_fivec_interval_subparser(subparsers)
    add_fivec_combine_replicates_subparser(subparsers)
    add_fends_subparser(subparsers)
    add_hicdataset_subparser(subparsers)
    add_hicproject_subparser(subparsers)
    add_hicnormalize_subparser(subparsers)
    add_complete_hic_subparser(subparsers)
    add_hic_heatmap_subparser(subparsers)
    add_hic_mrheatmap_subparser(subparsers)
    add_hic_interval_subparser(subparsers)
    add_hic_combine_replicates_subparser(subparsers)
    add_quasar_subparser(subparsers)
    return parser

def add_connect_subparser(subparsers):
    """Add command 'connect' arguments to parser."""
    parser = subparsers.add_parser("connect", help="Link (or re-link) HiFive files, partition to data or data to project. Partition file links will automatically be updated in project files with linking data files to them.")
    parser.add_argument(dest="type", type=str, choices=['fragments', 'fends', '5c-data', 'hic-data'],
        help="The type of file being linked to.")
    parser.add_argument(dest="target", type=str,
        help="The file to be linked to (partition or data).")
    parser.add_argument(dest="source", type=str,
        help="The file getting an updated link (data or project).")
    add_silent_argument(parser)
    return

def add_fragments_subparser(subparsers):
    """Add command 'fragments' arguments to parser."""
    parser = subparsers.add_parser("fragments", help="HiFive Fragment Creation Function: Read target RE fragment data from a BED file and create a HiFive 'Fragment' file.")
    parser.add_argument("-r", "--re", dest="re", required=False, default=None, type=str,
        help="The name of the restriction enzyme used to produce the RE fragments.")
    parser.add_argument("-g", "--genome", dest="genome", required=False, default=None, type=str,
        help="The name of the genome RE fragments correspond to.")
    parser.add_argument(dest="bed", type=str,
        help="The name of the BED file containing 5C fragment positions.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiFive Fragments to.")
    add_silent_argument(parser)
    return

def add_fivecdataset_subparser(subparsers):
    """Add command '5c-data' arguments to parser."""
    parser = subparsers.add_parser("5c-data", help="HiFive 5C Dataset Creation Function: Read 5C data from either sets of paired-end BAM files or tabular counts files and create a HiFive 'FiveCData' file.")
    infile_group = parser.add_mutually_exclusive_group(required=True)
    infile_group.add_argument("-B", "--bam", dest="bam", nargs=2, action='append',
        help="A pair of bam read end files from a single sequencing run. For multiple runs, this option can be passed multiple times.")
    infile_group.add_argument("-C", "--count", dest="count", action='append',
        help="A tab-separated text file containing a pair of fragment names and the number of observed reads for that pair (fragment1 fragment2 count), one per line. For multiple files, this option can be passed multiple times.")
    parser.add_argument(dest="fragment", type=str,
        help="The file name of an appropriate HiFive Fragment file.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiFive FiveCData to.")
    add_silent_argument(parser)
    return

def add_fivecproject_subparser(subparsers):
    """Add command '5c-project' arguments to parser."""
    parser = subparsers.add_parser("5c-project", help="HiFive 5C Project Creation Function: Creata a HiFive 5C project associated with a specific 5C dataset, filter out low-coverage fragments and calculate the distance-dependence function.")
    parser.add_argument("-f", "--min-interactions", dest="minint", required=False, type=int, default=20,
        action='store', help="The minimum number of interactions needed for valid fragment. [default: %(default)s]")
    parser.add_argument("-m", "--min-distance", dest="mindist", required=False, type=int, default=0,
        action='store', help="The minimum interaction distance to include in fragment filtering. [default: %(default)s]")
    parser.add_argument("-x", "--max-distance", dest="maxdist", required=False, type=int, default=None,
        action='store', help="The maximum interaction distance to include in fragment filtering (None or zero indicate no maximum). [default: %(default)s]")
    parser.add_argument(dest="data", type=str,
        help="The file name of an appropriate HiFive FiveCData file.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiFive FiveC project to.")
    add_silent_argument(parser)
    return

def add_fivecnormalize_subparser(subparsers):
    """Add command '5c-normalize' arguments to parser."""
    parser = subparsers.add_parser("5c-normalize", help="HiFive 5C Project Normalization Function: Learn correction parameters for a HiFive 5C project.", epilog="For command line options of each normalization approach, type: %(prog)s <ALGORITHM> -h")
    subparser2 = parser.add_subparsers(dest='algorithm')
    prob_parser = subparser2.add_parser("probability", help="A probability model based approach for learning correction values associated with each fragment. Learning is accomplished using gradient descent.")
    exp_parser = subparser2.add_parser("express", help="An appoximation based approach for learning correction values associated with each fragment. Learning is accomplished using a variation of matrix balancing.")
    bin_parser = subparser2.add_parser("binning", help="A multivariate binning probability model-based approach for learning correction values associated with fragment characteristics. Learning is accomplished using the Broyden-Fletcher-Goldfarb-Shanno algorithm.")
    binprob_parser = subparser2.add_parser("binning-probability", help="A chained-correction approach first learning fragment characteristic corrections and applying them prior to learning fragment-associated correction values via a probability model.")
    binexp_parser = subparser2.add_parser("binning-express", help="A chained-correction approach first learning fragment characteristic corrections and applying them prior to learning fragment-associated correction values via a matrix-balancing approximation.")
    probbin_parser = subparser2.add_parser("probability-binning", help="A chained-correction approach first learning fragment-associated correction values via a probability model and applying them prior to learning fragment characteristic corrections.")
    expbin_parser = subparser2.add_parser("express-binning", help="A chained-correction approach first learning fragment-associated correction values via a matrix-balancing approximation and applying them prior to learning fragment characteristic corrections.")
    add_fivec_normalize_group(prob_parser)
    add_fivec_probability_group(prob_parser)
    add_fivec_normalize_group(exp_parser)
    add_fivec_express_group(exp_parser)
    add_fivec_normalize_group(bin_parser)
    add_fivec_binning_group(bin_parser)
    add_fivec_normalize_group(binprob_parser)
    add_fivec_binning_group(binprob_parser)
    add_fivec_probability_group(binprob_parser)
    add_fivec_normalize_group(binexp_parser)
    add_fivec_binning_group(binexp_parser)
    add_fivec_express_group(binexp_parser)
    add_fivec_normalize_group(probbin_parser)
    add_fivec_probability_group(probbin_parser)
    add_fivec_binning_group(probbin_parser)
    add_fivec_normalize_group(expbin_parser)
    add_fivec_express_group(expbin_parser)
    add_fivec_binning_group(expbin_parser)
    return

def add_complete_fivec_subparser(subparsers):
    """Add command '5c-complete' arguments to parser."""
    parser = subparsers.add_parser("5c-complete", help="HiFive 5C Project Complete Analysis Function: Create all necessary files (Fragment, Data, and Project) and learn correction parameters for a HiFive 5C project.", epilog="For command line options of each normalization approach, type: %(prog)s <ALGORITHM> -h")
    subparser2 = parser.add_subparsers(dest='algorithm')
    prob_parser = subparser2.add_parser("probability", help="A probability model based approach for learning correction values associated with each fragment. Learning is accomplished using gradient descent.")
    exp_parser = subparser2.add_parser("express", help="An appoximation based approach for learning correction values associated with each fragment. Learning is accomplished using a variation of matrix balancing.")
    bin_parser = subparser2.add_parser("binning", help="A multivariate binning probability model-based approach for learning correction values associated with fragment characteristics. Learning is accomplished using the Broyden-Fletcher-Goldfarb-Shanno algorithm.")
    binprob_parser = subparser2.add_parser("binning-probability", help="A chained-correction approach first learning fragment characteristic corrections and applying them prior to learning fragment-associated correction values via a probability model.")
    binexp_parser = subparser2.add_parser("binning-express", help="A chained-correction approach first learning fragment characteristic corrections and applying them prior to learning fragment-associated correction values via a matrix-balancing approximation.")
    probbin_parser = subparser2.add_parser("probability-binning", help="A chained-correction approach first learning fragment-associated correction values via a probability model and applying them prior to learning fragment characteristic corrections.")
    expbin_parser = subparser2.add_parser("express-binning", help="A chained-correction approach first learning fragment-associated correction values via a matrix-balancing approximation and applying them prior to learning fragment characteristic corrections.")
    add_fivec_complete_group(prob_parser)
    add_fivec_probability_group(prob_parser)
    add_fivec_complete_group(exp_parser)
    add_fivec_express_group(exp_parser)
    add_fivec_complete_group(bin_parser)
    add_fivec_binning_group(bin_parser)
    add_fivec_complete_group(binprob_parser)
    add_fivec_binning_group(binprob_parser)
    add_fivec_probability_group(binprob_parser)
    add_fivec_complete_group(binexp_parser)
    add_fivec_binning_group(binexp_parser)
    add_fivec_express_group(binexp_parser)
    add_fivec_complete_group(probbin_parser)
    add_fivec_probability_group(probbin_parser)
    add_fivec_binning_group(probbin_parser)
    add_fivec_complete_group(expbin_parser)
    add_fivec_express_group(expbin_parser)
    add_fivec_binning_group(expbin_parser)
    return

def add_fivec_normalize_group(subparser):
    """Add common 5C normalization options to subparser."""
    subparser.add_argument("-m", "--min-distance", dest="mindist", required=False, type=int, default=0,
        action='store', help="The minimum interaction distance to include in normalization. [default: %(default)s]")
    subparser.add_argument("-x", "--max-distance", dest="maxdist", required=False, type=int, default=None,
        action='store', help="The maximum interaction distance to include in normalization (None or zero indicate no maximum). [default: %(default)s]")
    subparser.add_argument("-r", "--regions", dest="regions", required=False, type=str, default=None,
        action='store', help="A comma-separated list of regions to learn correction values for (None indicates all regions). [default: %(default)s]")
    subparser.add_argument("-o", "--output-file", dest="output", required=False, type=str, default=None,
        action='store', help="An alternate filename to save the normalized project to. If not given, the original project file will be overwritten. [default: %(default)s]")
    subparser.add_argument(dest="project", type=str,
        help="The name of the HiFive FiveC project to normalize.")
    add_silent_argument(subparser)
    return

def add_fivec_complete_group(subparser):
    """Add common 5C complete analysis options to subparser."""
    subparser.add_argument("--re", dest="re", required=False, default=None, type=str,
        help="The name of the restriction enzyme used to produce the RE fragments.")
    subparser.add_argument("--genome", dest="genome", required=False, default=None, type=str,
        help="The name of the genome RE fragments correspond to.")
    subparser.add_argument(dest="bed", type=str,
        help="The name of the BED file containing 5C fragment positions.")
    infile_group = subparser.add_mutually_exclusive_group(required=True)
    infile_group.add_argument("-B", "--bam", dest="bam", nargs=2, action='append',
        help="A pair of bam read end files from a single sequencing run. For multiple runs, this option can be passed multiple times.")
    infile_group.add_argument("-C", "--count", dest="count", action='append',
        help="A tab-separated text file containing a pair of fragment names and the number of observed reads for that pair (fragment1 fragment2 count), one per line. For multiple files, this option can be passed multiple times.")
    subparser.add_argument("-f", "--min-interactions", dest="minint", required=False, type=int, default=20,
        action='store', help="The minimum number of interactions needed for valid fragment. [default: %(default)s]")
    subparser.add_argument("-m", "--min-distance", dest="mindist", required=False, type=int, default=0,
        action='store', help="The minimum interaction distance to include in normalization. [default: %(default)s]")
    subparser.add_argument("-x", "--max-distance", dest="maxdist", required=False, type=int, default=None,
        action='store', help="The maximum interaction distance to include in normalization (None or zero indicate no maximum). [default: %(default)s]")
    subparser.add_argument("-r", "--regions", dest="regions", required=False, type=str, default=None,
        action='store', help="A comma-separated list of regions to learn correction values for (None indicates all regions). [default: %(default)s]")
    add_silent_argument(subparser)
    outfile_group = subparser.add_mutually_exclusive_group(required=True)
    outfile_group.add_argument("-P", "--prefix", dest="prefix", type=str, default=None,
        help="A prefix for the output files (*.frags, *.fcd, and *.fcp).")
    outfile_group.add_argument("-o", "--output", dest="output", type=str, nargs=3, default=None,
        help="File names for the Fragment, FiveCData, and FiveC project files to be created.")
    return

def add_fivec_probability_group(subparser):
    """Add 5C probability normalization options to subparser."""
    subparser.add_argument("-b", "--max-iterations", dest="probiter", required=False, type=int, default=1000,
        action='store', help="The maximum number of iterations to carry on gradient descent for. [default: %(default)s]")
    subparser.add_argument("-g", "--min-change", dest="change", required=False, type=float, default=0.0005,
        action='store', help="The minimum allowable absolute gradient size to coninute learning process. [default: %(default)s]")
    subparser.add_argument("-p", "--precalculate", dest="precalc", required=False, default=False,
        action='store_true', help="Precalculate correction values from fragment means for the probability algorithm. [default: %(default)s]")
    subparser.add_argument("-l", "--learning-step", dest="step", required=False, type=float, default=0.5,
        action='store', help="The scaling factor for decreasing learning rate by if step doesn't meet armijo criterion. [default: %(default)s]")
    return

def add_fivec_express_group(subparser):
    """Add 5C express normalization options to subparser."""
    subparser.add_argument("-e", "--express-iterations", dest="expiter", required=False, type=int, default=1000,
        action='store', help="The number of iterations to run the express learning phase for. [default: %(default)s]")
    subparser.add_argument("-d", "--remove-distance", dest="nodist", required=False, default=False,
        action='store_true', help="Remove the distant-dependent portion of the signal prior to learning corrections with the express algorithm. [default: %(default)s]")
    subparser.add_argument("-w", "--express-reads", dest="expreads", required=False, type=str, default='cis',
        choices=['cis', 'trans', 'all'], help="Which set of reads to use for express normalization. [default: %(default)s]")
    subparser.add_argument("-k", "--logged", dest="logged", required=False, default=False,
        action='store_true', help="Use log counts. [default: %(default)s]")
    subparser.add_argument("-z", "--knight-ruiz", dest="kr", required=False, default=False,
        action='store_true', help="Use Knight-Ruiz algorithm for unweighted matrix balancing. [default: %(default)s]")
    return

def add_fivec_binning_group(subparser):
    """Add 5C binning normalization options to subparser."""
    subparser.add_argument("-i", "--binning-iterations", dest="biniter", required=False, type=int, default=1000,
        action='store', help="The maximum number of iterations to run binning modeling for. [default: %(default)s]")
    subparser.add_argument("-t", "--learning-threshold", dest="threshold", required=False, type=float, default=1.0,
        action='store', help="The learning threshold cutoff for binning algorithm. [default: %(default)s]")
    subparser.add_argument("-y", "--binning-reads", dest="binreads", required=False, type=str, default='cis',
        choices=['cis', 'trans', 'all'], help="Which set of reads to use for binning normalization. [default: %(default)s]")
    subparser.add_argument("-v", "--model", dest="model", required=False, type=str, default="len",
        action='store', help="A comma-separated list of parameters to include in binning model. [default: %(default)s]")
    subparser.add_argument("-n", "--model-bins", dest="modelbins", required=False, type=str, default="10",
        action='store', help="A comma-separated list of the number of bins to separate model parameters into. [default: %(default)s]")
    subparser.add_argument("-u", "--parameter-types", dest="parameters", required=False, type=str, default="even",
        action='store', help="A comma-separated list parameter types ('even' or 'fixed', depending on whether bins should contain equal numbers of fends or be equally spaced along the parameter range, and a '-const' suffix to indicate the values are not to be optimized). [default: %(default)s]")
    return

def add_fivec_combine_replicates_subparser(subparsers):
    """Add command '5c-combine-replicates' arguments to parser."""
    parser = subparsers.add_parser("5c-combine-replicates",
        help="HiFive Data Function: Combine two HiFive 5C datasets into a single dataset without needing to reload data.")
    parser.add_argument(dest="output", type=str,
        help="Name of file in which to write combined data.")
    parser.add_argument(dest="replicate", nargs="+",
        help="Two of more HiFive 5C data files to be combined. These files need to have been generated using the same Fragment file.")
    add_silent_argument(parser)
    return

def add_fivec_heatmap_subparser(subparsers):
    """Add command '5c-heatmap' arguments to parser."""
    parser = subparsers.add_parser("5c-heatmap",
        help="HiFive Binning Function: Create a heatmap HDF5 file containing data from a HiFive FiveC project.")
    parser.add_argument("-b", "--binsize", dest="binsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into. [default: %(default)s]")
    parser.add_argument("-t", "--trans", dest="trans", default=False, action="store_true",
        help="Calculate and include trans interactions in heatmaps. [default: %(default)s]")
    parser.add_argument("-r", "--regions", dest="regions", default=None, type=str,
        help="A comma-separated list of regions to include in the heatmaps (None indicates all regions). [default: %(default)s]")
    parser.add_argument("-d", "--datatype", dest="datatype", default="fragment",
        help="Which corrections (if any) to apply to counts. [default: %(default)s]",
        choices=["raw", "fragment", "distance", "enrichment", "expected"])
    parser.add_argument("-a", "--arraytype", dest="arraytype", default="full",
        help="The type of array layout to store data in ('compact' only valid for unbinned data). [default: %(default)s]",
        choices=["compact", "full"])
    parser.add_argument("-F", "--format", dest="format", default="hdf5",
        help="Format of output. [default: %(default)s]", choices=["hdf5", "txt", "npz"])
    parser.add_argument("-y", "--dynamically-bin", dest="dynamic", default=False, action="store_true",
        help="Dynamically bin heatmap.")
    parser.add_argument("-x", "--expansion-binsize", dest="expbinsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into for expanding under-populated bins. [default: %(default)s]")
    parser.add_argument("-f", "--minobservations", dest="minobs", default=20, type=int,
        help="The minimum number of observed reads in a bin for it to be considered valid. [default: %(default)s]")
    parser.add_argument("-g", "--search-distance", dest="search", default=0, type=int,
        help="The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [default: %(default)s]")
    parser.add_argument("-v", "--remove-failed", dest="remove", default=False, action="store_true",
        help="If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.")
    parser.add_argument("-i", "--image-file", dest="image", default=None, type=str,
        help="Save the data as an image to this file.")
    parser.add_argument("-p", "--pdf", dest="pdf", default=False, action="store_true",
        help="Format the image in PDF format. [default: %(default)s]")
    parser.add_argument("-l", "--legend", dest="legend", default=False, action="store_true",
        help="Add color scale to the plot (pdf format only). [default: %(default)s]")
    parser.add_argument("-n", "--names", dest="names", default=False, action="store_true",
        help="Add region labels to the plot (pdf format only). [default: %(default)s]")
    parser.add_argument("-k", "--keyword", dest="keywords", default=[], type=str, action='append',
        help="Additional keyword arguments to pass to plotting function.")
    add_silent_argument(parser)
    parser.add_argument(dest="project", type=str,
        help="The name of a HiFive FiveC project file to pull data from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write 5C heatmaps to.")
    return

def add_fivec_interval_subparser(subparsers):
    """Add command '5c-interval' arguments to parser."""
    parser = subparsers.add_parser("5c-interval",
        help="HiFive Binning Function: Create a tabular interaction file containing data from a HiFive FiveC project. Data are a genomic-interval format (chr1 start1 stop1 chr2 start2 stop2).")
    parser.add_argument("-c", "--region", dest="region", default=None, required=True, type=int,
        help="The region from which to pull interaction data from.")
    parser.add_argument("-s", "--start", dest="start", default=None, required=False, type=int,
        help="The start coordinate of the pulled region to return. (None indicates the first valid bin in the region) [default %(default)s]")
    parser.add_argument("-e", "--stop", dest="stop", default=None, required=False, type=int,
        help="The stop coordinate of the pulled region to return. (None indicates the last valid bin in the region) [default %(default)s]")
    parser.add_argument("--region2", dest="region2", default=None, required=False, type=int,
        help="The second region from which to pull interaction data from if trans interactions are desired.")
    parser.add_argument("--start2", dest="start2", default=None, required=False, type=int,
        help="The start coordinate of the second pulled region to return. (None indicates the first valid bin in the region) [default %(default)s]")
    parser.add_argument("--stop2", dest="stop2", default=None, required=False, type=int,
        help="The stop coordinate of the second pulled region to return. (None indicates the last valid bin in the region) [default %(default)s]")
    parser.add_argument("-b", "--binsize", dest="binsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into. [default: %(default)s]")
    parser.add_argument("-d", "--data-type", dest="datatype", default="fragment",
        help="Which corrections (if any) to apply to counts. [default: %(default)s]",
        choices=["raw", "fragment", "distance", "enrichment", "expected"])
    parser.add_argument("-y", "--dynamically-bin", dest="dynamic", default=False, action="store_true",
        help="Dynamically bin heatmap.")
    parser.add_argument("-x", "--expansion-binsize", dest="expbinsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into for expanding under-populated bins. [default: %(default)s]")
    parser.add_argument("-f", "--minobservations", dest="minobs", default=20, type=int,
        help="The minimum number of observed reads in a bin for it to be considered valid. [default: %(default)s]")
    parser.add_argument("-g", "--search-distance", dest="search", default=0, type=int,
        help="The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [default: %(default)s]")
    parser.add_argument("-v", "--remove-failed", dest="remove", default=False, action="store_true",
        help="If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.")
    parser.add_argument("-i", "--image-file", dest="image", default=None, type=str,
        help="Save the data as an image to this file.")
    parser.add_argument("-p", "--pdf", dest="pdf", default=False, action="store_true",
        help="Format the image in PDF format. [default: %(default)s]")
    parser.add_argument("-r", "--rotate", dest="rotate", default=False, action="store_true",
        help="Rotate the plot 45 degrees (binned only). [default: %(default)s]")
    parser.add_argument("-t", "--ticks", dest="ticks", default=False, action="store_true",
        help="Add tick marks and labels to the plot (pdf format and binned only). [default: %(default)s]")
    parser.add_argument("-l", "--legend", dest="legend", default=False, action="store_true",
        help="Add color scale to the plot (pdf format only). [default: %(default)s]")
    parser.add_argument("-k", "--keyword", dest="keywords", default=[], type=str, action='append',
        help="Additional keyword arguments to pass to plotting function.")
    add_silent_argument(parser)
    parser.add_argument(dest="project", type=str,
        help="The name of a HiFive FiveC project file to pull data from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write 5C interval to.")
    return

def add_fends_subparser(subparsers):
    """Add command 'fends' arguments to parser."""
    parser = subparsers.add_parser("fends", help="HiFive Fend Creation Function: Read RE fragment data from either a BED file or a HiCPipe-style tabular fend file and create a HiFive 'Fend' file.")
    infile_group = parser.add_mutually_exclusive_group(required=True)
    infile_group.add_argument("-F", "--fend", dest="fend", type=str, default=None,
        help="Fend file in HiCPipe-compatible tabular format, optionally containing fend characteristics (GC content and/or mappability).")
    infile_group.add_argument("-B", "--bed", dest="bed", type=str, default=None,
        help="Bed file containing either RE fragment boundary data or RE cutsites.")
    infile_group.add_argument("-L", "--length", dest="length", type=str, default=None,
        help="Text file containing chromosome names and lengths.")
    parser.add_argument("--binned", dest="binned", type=int, default=None,
        help="Interval to bin data into.")
    parser.add_argument("-r", "--re", dest="re", required=False, default=None, type=str,
        help="The name of the restriction enzyme used to produce the RE fragments.")
    parser.add_argument("-g", "--genome", dest="genome", required=False, default=None, type=str,
        help="The name of the genome RE fragments correspond to.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiFive Fends to.")
    add_silent_argument(parser)
    return

def add_hicdataset_subparser(subparsers):
    """Add command 'hic-data' arguments to parser."""
    parser = subparsers.add_parser("hic-data", help="HiFive HiC Dataset Creation Function: Read HiC data from sets of paired-end BAM files, tabular read position (RAW) files, or HiCPipe-compatible fend-pair counts (MAT) files and create a HiFive 'HiCData' file.")
    infile_group = parser.add_mutually_exclusive_group(required=True)
    infile_group.add_argument("-S", "--bam", dest="bam", nargs=2, action='append',
        help="A pair of BAM read end files from a single sequencing run. For multiple runs, this option can be passed multiple times.")
    infile_group.add_argument("-R", "--raw", dest="raw", action='append',
        help="A tab-separated text file containing pairs of read ends (chr1 pos1 strand1 chr2 pos2 strand2), one per line. For multiple files, this option can be passed multiple times.")
    infile_group.add_argument("-M", "--mat", dest="mat", action='store',
        help="A HiCPipe-style tabular MAT file containing fend pair counts.")
    infile_group.add_argument("-X", "--matrix", dest="matrix", action='store',
        help="Binned matrix containing summed fend interactions.")
    parser.add_argument("-i", "--insert", dest="insert", required=False, type=int, default=500,
        help="The maximum allowable distance sum between both fend ends and cutsites. [default: %(default)s]")
    parser.add_argument("--skip-duplicate-filtering", dest="skipdups", required=False, default=False,
        action='store_true', help="Skip filtering of PCR duplicates. [default: %(default)s]")
    parser.add_argument(dest="fend", type=str,
        help="The file name of an appropriate HiFive Fend file.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiFive HiCData to.")
    add_silent_argument(parser)
    return

def add_hicproject_subparser(subparsers):
    """Add command 'hic-project' arguments to parser. This command is MPI-compatible."""
    parser = subparsers.add_parser("hic-project", help="HiFive HiC Project Creation Function: Creata a HiFive HiC project associated with a specific HiC dataset, filter out low-coverage fends and calculate the distance-dependence function.")
    parser.add_argument("-c", "--chromosomes", dest="chroms", required=False, type=str, default=None,
        action='store', help="A comma-separated list of chromosomes to learn correction values for (None indicates all chromosomes). [default: %(default)s]")
    parser.add_argument("-f", "--min-interactions", dest="minint", required=False, type=int, default=20,
        action='store', help="The minimum number of interactions needed for valid fragment. [default: %(default)s]")
    parser.add_argument("-m", "--min-distance", dest="mindist", required=False, type=int, default=0,
        action='store', help="The minimum interaction distance to include in fragment filtering. [default: %(default)s]")
    parser.add_argument("-x", "--max-distance", dest="maxdist", required=False, type=int, default=0,
        action='store', help="The maximum interaction distance to include in fragment filtering (None or zero indicate no maximum). [default: %(default)s]")
    parser.add_argument("-j", "--min-binsize", dest="minbin", required=False, type=int, default=1000,
        action='store', help="The smallest interaction distance bin size for the distance-dependence function. [default: %(default)s]")
    parser.add_argument("-n", "--num-bins", dest="numbins", required=False, type=int, default=100,
        action='store', help="The number of bins to partion the interaction distance range into for distance-dependence function. A value of zero indicates that finding the distance dependence function should be skipped. [default: %(default)s]")
    parser.add_argument(dest="data", type=str,
        help="The file name of an appropriate HiFive HiCData file.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiFive HiC project to.")
    add_silent_argument(parser)
    return

def add_hicnormalize_subparser(subparsers):
    """Add command 'hic-normalize' arguments to parser."""
    parser = subparsers.add_parser("hic-normalize", help="HiFive HiC Project Normalization Function: Learn correction parameters for a HiFive HiC project.", epilog="For command line options of each normalization approach, type: %(prog)s <ALGORITHM> -h")
    subparser2 = parser.add_subparsers(dest='algorithm')
    prob_parser = subparser2.add_parser("probability", help="A probability model based approach for learning correction values associated with each fend. Learning is accomplished using gradient descent.")
    exp_parser = subparser2.add_parser("express", help="An appoximation based approach for learning correction values associated with each fend. Learning is accomplished using a variation of matrix balancing.")
    bin_parser = subparser2.add_parser("binning", help="A multivariate binning probability model-based approach for learning correction values associated with fend characteristics. Learning is accomplished using the Broyden-Fletcher-Goldfarb-Shanno algorithm.")
    binprob_parser = subparser2.add_parser("binning-probability", help="A chained-correction approach first learning fend characteristic corrections and applying them prior to learning fend-associated correction values via a probability model.")
    binexp_parser = subparser2.add_parser("binning-express", help="A chained-correction approach first learning fend characteristic corrections and applying them prior to learning fend-associated correction values via a matrix-balancing approximation.")
    add_hic_normalize_group(prob_parser)
    add_hic_probability_group(prob_parser)
    add_hic_normalize_group(exp_parser)
    add_hic_express_group(exp_parser)
    exp_parser.add_argument("-f", "--min-interactions", dest="minint", required=False, type=int, default=20,
        action='store', help="The minimum number of interactions for fend filtering, if refiltering is required due to distance cutoff or selected reads. [default: %(default)s]")
    add_hic_normalize_group(bin_parser)
    add_hic_binning_group(bin_parser)
    add_hic_normalize_group(binprob_parser)
    add_hic_probability_group(binprob_parser)
    add_hic_binning_group(binprob_parser)
    add_hic_normalize_group(binexp_parser)
    add_hic_express_group(binexp_parser)
    binexp_parser.add_argument("-f", "--min-interactions", dest="minint", required=False, type=int, default=20,
        action='store', help="The minimum number of interactions for fend filtering, if refiltering is required due to distance cutoff or selected reads. [default: %(default)s]")
    add_hic_binning_group(binexp_parser)
    return

def add_complete_hic_subparser(subparsers):
    """Add command 'hic-complete' arguments to parser."""
    parser = subparsers.add_parser("hic-complete", help="HiFive HiC Project Complete Analysis Function: Create all necessary files (Fend, Data, and Project) and learn correction parameters for a HiFive HiC project.", epilog="For command line options of each normalization approach, type: %(prog)s <ALGORITHM> -h")
    subparser2 = parser.add_subparsers(dest='algorithm')
    prob_parser = subparser2.add_parser("probability", help="A probability model based approach for learning correction values associated with each fend. Learning is accomplished using gradient descent.")
    exp_parser = subparser2.add_parser("express", help="An appoximation based approach for learning correction values associated with each fend. Learning is accomplished using a variation of matrix balancing.")
    bin_parser = subparser2.add_parser("binning", help="A multivariate binning probability model-based approach for learning correction values associated with fend characteristics. Learning is accomplished using the Broyden-Fletcher-Goldfarb-Shanno algorithm.")
    binprob_parser = subparser2.add_parser("binning-probability", help="A chained-correction approach first learning fend characteristic corrections and applying them prior to learning fend-associated correction values via a probability model.")
    binexp_parser = subparser2.add_parser("binning-express", help="A chained-correction approach first learning fend characteristic corrections and applying them prior to learning fend-associated correction values via a matrix-balancing approximation.")
    add_complete_hic_group(prob_parser)
    add_hic_probability_group(prob_parser)
    add_complete_hic_group(exp_parser)
    add_hic_express_group(exp_parser)
    add_complete_hic_group(bin_parser)
    add_hic_binning_group(bin_parser)
    add_complete_hic_group(binprob_parser)
    add_hic_probability_group(binprob_parser)
    add_hic_binning_group(binprob_parser)
    add_complete_hic_group(binexp_parser)
    add_hic_express_group(binexp_parser)
    add_hic_binning_group(binexp_parser)
    return

def add_complete_hic_group(subparser):
    """Add common HiC complete analysis options to subparser."""
    fendfile_group = subparser.add_mutually_exclusive_group(required=True)
    fendfile_group.add_argument("-F", "--fend", dest="fend", type=str, default=None,
        help="Fend file in HiCPipe-compatible tabular format, optionally containing fend characteristics (GC content and/or mappability).")
    fendfile_group.add_argument("-B", "--bed", dest="bed", type=str, default=None,
        help="Bed file containing either RE fragment boundary data or RE cutsites.")
    fendfile_group.add_argument("-L", "--length", dest="length", type=str, default=None,
        help="Text file containing chromosome names and lengths.")
    subparser.add_argument("--binned", dest="binned", type=int, default=None,
        help="Interval to bin data into. If set to zero, indicates non-uniform binning to be read from bed file.")
    subparser.add_argument("--re", dest="re", required=False, default=None, type=str,
        help="The name of the restriction enzyme used to produce the RE fragments.")
    subparser.add_argument("--genome", dest="genome", required=False, default=None, type=str,
        help="The name of the genome RE fragments correspond to.")
    infile_group = subparser.add_mutually_exclusive_group(required=True)
    infile_group.add_argument("-S", "--bam", dest="bam", nargs=2, action='append',
        help="A pair of BAM read end files from a single sequencing run. For multiple runs, this option can be passed multiple times.")
    infile_group.add_argument("-R", "--raw", dest="raw", action='append',
        help="A tab-separated text file containing pairs of read ends (chr1 pos1 strand1 chr2 pos2 strand2), one per line. For multiple files, this option can be passed multiple times.")
    infile_group.add_argument("-M", "--mat", dest="mat", type=str, action='store',
        help="A HiCPipe-style tabular MAT file containing fend pair counts.")
    infile_group.add_argument("-X", "--matrix", dest="matrix", action='append',
        help="Binned matrix containing summed fend interactions.")
    subparser.add_argument("-i", "--insert", dest="insert", required=False, type=int, default=500,
        help="The maximum allowable distance sum between both fend ends and cutsites. [default: %(default)s]")
    subparser.add_argument("--skip-duplicate-filtering", dest="skipdups", required=False, default=False,
        action='store_true', help="Skip filtering of PCR duplicates. [default: %(default)s]")
    subparser.add_argument("-f", "--min-interactions", dest="minint", required=False, type=int, default=20,
        action='store', help="The minimum number of interactions needed for valid fragment. [default: %(default)s]")
    subparser.add_argument("-m", "--min-distance", dest="mindist", required=False, type=int, default=0,
        action='store', help="The minimum interaction distance to include in fragment filtering. [default: %(default)s]")
    subparser.add_argument("-x", "--max-distance", dest="maxdist", required=False, type=int, default=None,
        action='store', help="The maximum interaction distance to include in fragment filtering (None or zero indicate no maximum). [default: %(default)s]")
    subparser.add_argument("-c", "--chromosomes", dest="chroms", required=False, type=str, default=None,
        action='store', help="A comma-separated list of chromosomes to learn correction values for (None indicates all chromosomes). [default: %(default)s]")
    subparser.add_argument("-j", "--min-binsize", dest="minbin", required=False, type=int, default=1000,
        action='store', help="The smallest interaction distance bin size for the distance-dependence function. [default: %(default)s]")
    subparser.add_argument("-n", "--num-bins", dest="numbins", required=False, type=int, default=100,
        action='store', help="The number of bins to partion the interaction distance range into for distance-dependence function. A value of zero indicates that finding the distance dependence function should be skipped. [default: %(default)s]")
    add_silent_argument(subparser)
    outfile_group = subparser.add_mutually_exclusive_group(required=True)
    outfile_group.add_argument("-P", "--prefix", dest="prefix", type=str, default=None,
        help="A prefix for the output files (*.fends, *.hcd, and *.hcp).")
    outfile_group.add_argument("-o", "--output", dest="output", type=str, nargs=3, default=None,
        help="File names for the Fend, HiCData, and HiC project files to be created.")
    return

def add_hic_normalize_group(subparser):
    """Add common HiC normalization options to subparser."""
    subparser.add_argument("-m", "--min-distance", dest="mindist", required=False, type=int, default=0,
        action='store', help="The minimum interaction distance to include in normalization. [default: %(default)s]")
    subparser.add_argument("-x", "--max-distance", dest="maxdist", required=False, type=int, default=None,
        action='store', help="The maximum interaction distance to include in normalization (None or zero indicate no maximum). [default: %(default)s]")
    subparser.add_argument("-c", "--chromosomes", dest="chroms", required=False, type=str, default=None,
        action='store', help="A comma-separated list of chromosomes to learn correction values for (None indicates all chromosomes). [default: %(default)s]")
    subparser.add_argument("-o", "--output-file", dest="output", required=False, type=str, default=None,
        action='store', help="An alternate filename to save the normalized project to. If not given, the original project file will be overwritten. [default: %(default)s]")
    subparser.add_argument(dest="project", type=str,
        help="The name of the HiFive HiC project to normalize.")
    add_silent_argument(subparser)
    return

def add_hic_probability_group(subparser):
    """Add HiC probability normalization options to subparser."""
    subparser.add_argument("-b", "--max-iterations", dest="probiter", required=False, type=int, default=1000,
        action='store', help="The maximum number of iterations to carry on gradient descent for. [default: %(default)s]")
    subparser.add_argument("-g", "--min-change", dest="change", required=False, type=float, default=0.0005,
        action='store', help="The minimum allowable absolute gradient size to coninute learning process. [default: %(default)s]")
    subparser.add_argument("-p", "--precalculate", dest="precalc", required=False, default=False,
        action='store_true', help="Precalculate correction values from fend means for the probability algorithm. [default: %(default)s]")
    subparser.add_argument("-l", "--learning-step", dest="step", required=False, type=float, default=0.4,
        action='store', help="The scaling factor for decreasing learning rate by if step doesn't meet armijo criterion. [default: %(default)s]")
    subparser.add_argument("-a", "--probability-model", dest="probmodel", required=False, type=str, default='binomial',
        choices=['binomial', 'poisson'], help="Which probability model to use for normalization. [default: %(default)s]")
    return

def add_hic_express_group(subparser):
    """Add HiC express normalization options to subparser."""
    subparser.add_argument("-e", "--express-iterations", dest="expiter", required=False, type=int, default=1000,
        action='store', help="The minimum number of iterations to run the express learning phase for. [default: %(default)s]")
    subparser.add_argument("-g", "--min-change", dest="change", required=False, type=float, default=0.0001,
        action='store', help="The maximum allowable change per iteration in fend correction parameter values allowable to terminate learning. [default: %(default)s]")
    subparser.add_argument("-d", "--remove-distance", dest="nodist", required=False, default=False,
        action='store_true', help="Remove the distant-dependent portion of the signal prior to learning corrections with the express algorithm. [default: %(default)s]")
    subparser.add_argument("-w", "--express-reads", dest="expreads", required=False, type=str, default='cis',
        choices=['cis', 'trans', 'all'], help="Which set of reads to use for express normalization. [default: %(default)s]")
    subparser.add_argument("-k", "--binary", dest="binary", required=False, default=False,
        action='store_true', help="Use binary indicator instead of counts. [default: %(default)s]")
    subparser.add_argument("-z", "--knight-ruiz", dest="kr", required=False, default=False,
        action='store_true', help="Use Knight-Ruiz algorithm for unweighted matrix balancing. [default: %(default)s]")
    return

def add_hic_binning_group(subparser):
    """Add HiC binning normalization options to subparser."""
    subparser.add_argument("-r", "--binning-iterations", dest="biniter", required=False, type=int, default=1000,
        action='store', help="The maximum number of iterations to run binning modeling for. [default: %(default)s]")
    subparser.add_argument("-t", "--learning-threshold", dest="threshold", required=False, type=float, default=1.0,
        action='store', help="The learning threshold cutoff for binning algorithm. [default: %(default)s]")
    subparser.add_argument("-y", "--binning-reads", dest="binreads", required=False, type=str, default='cis',
        choices=['cis', 'trans', 'all'], help="Which set of reads to use for binning normalization. [default: %(default)s]")
    subparser.add_argument("-v", "--model", dest="model", required=False, type=str, default="len,distance",
        action='store', help="A comma-separated list of parameters to include in binning model. [default: %(default)s]")
    subparser.add_argument("-s", "--model-bins", dest="modelbins", required=False, type=str, default="20,20",
        action='store', help="A comma-separated list of the number of bins to separate model parameters into. [default: %(default)s]")
    subparser.add_argument("-u", "--parameter-types", dest="parameters", required=False, type=str, default="even,fixed-const",
        action='store', help="A comma-separated list parameter types ('even' or 'fixed', depending on whether bins should contain equal numbers of fends or be equally spaced along the parameter range, and a '-const' suffix to indicate the values are not to be optimized). [default: %(default)s]")
    subparser.add_argument("--pseudocounts", dest="pseudo", required=False, type=int, default=None,
        action='store', help="The number of pseudo-counts to add to each bin prior to learning. [default: %(default)s]")
    return

def add_hic_combine_replicates_subparser(subparsers):
    """Add command 'hic-combine-replicates' arguments to parser."""
    parser = subparsers.add_parser("hic-combine-replicates",
        help="HiFive Data Function: Combine two HiFive HiC datasets into a single dataset without needing to reload data. These files need to have been generated using the same HiFive Fend file.")
    parser.add_argument(dest="replicate1", type=str,
        help="Name of the first HiFive HiCDataset replicate file to be combined.")
    parser.add_argument(dest="replicate2", type=str,
        help="Name of the second HiFive HiCDataset replicate file to be combined.")
    parser.add_argument(dest="output", type=str,
        help="Name of file in which to write combined data.")
    add_silent_argument(parser)
    return

def add_hic_heatmap_subparser(subparsers):
    """Add command 'hic-heatmap' arguments to parser."""
    parser = subparsers.add_parser("hic-heatmap",
        help="HiFive Binning Function: Create a heatmap HDF5 file containing data from a HiFive HiC project.")
    add_common_heatmap_arguments(parser)
    parser.add_argument("-b", "--binsize", dest="binsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into. [default: %(default)s]")
    parser.add_argument("-d", "--datatype", dest="datatype", default="fend",
        help="Which corrections (if any) to apply to counts. [default: %(default)s]",
        choices=["raw", "fend", "distance", "enrichment", "expected"])
    parser.add_argument("-F", "--format", dest="format", default="hdf5",
        help="Format of output. [default: %(default)s]",
        choices=["hdf5", "txt", "npz"])
    parser.add_argument("-y", "--dynamically-bin", dest="dynamic", default=False, action="store_true",
        help="Dynamically bin heatmap.")
    parser.add_argument("-x", "--expansion-binsize", dest="expbinsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into for expanding under-populated bins. [default: %(default)s]")
    parser.add_argument("-a", "--search-distance", dest="search", default=0, type=int,
        help="The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [default: %(default)s]")
    parser.add_argument("-v", "--remove-failed", dest="remove", default=False, action="store_true",
        help="If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.")
    parser.add_argument("-i", "--image-file", dest="image", default=None, type=str,
        help="Save the data as an image to this file.")
    parser.add_argument("-p", "--pdf", dest="pdf", default=False, action="store_true",
        help="Format the image in PDF format. [default: %(default)s]")
    parser.add_argument("-l", "--legend", dest="legend", default=False, action="store_true",
        help="Add color scale to the plot (pdf format only). [default: %(default)s]")
    parser.add_argument("-n", "--names", dest="names", default=False, action="store_true",
        help="Add chromosome labels to the plot (pdf format only). [default: %(default)s]")
    parser.add_argument("-k", "--keyword", dest="keywords", default=[], type=str, action='append',
        help="Additional keyword arguments to pass to plotting function.")
    add_silent_argument(parser)
    parser.add_argument(dest="project", type=str,
        help="The name of a HiFive HiC project file to pull data from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiC heatmaps to.")
    return

def add_hic_mrheatmap_subparser(subparsers):
    """Add command 'hic-mrheatmap' arguments to parser."""
    parser = subparsers.add_parser("hic-mrheatmap",
        help="HiFive Binning Function: Create a multi-resolution heatmap binary file containing data from a HiFive HiC project.")
    add_common_heatmap_arguments(parser)
    parser.add_argument("-B", "--maximum-binsize", dest="maxbin", default=1280000, type=int,
        help="The largest sized bin to use (minimum resolution), in base pairs. [default: %(default)s]")
    parser.add_argument("-b", "--minimum-binsize", dest="minbin", default=10000, type=int,
        help="The smallest sized bin to use (maximum resolution), in base pairs. [default: %(default)s]")
    parser.add_argument("-R", "--maximum-trans-binsize", dest="maxtransbin", default=None, type=int,
        help="The largest sized bin to use for trans interactions (minimum resolution), in base pairs. [default: %(default)s]")
    parser.add_argument("-r", "--minimum-trans-binsize", dest="mintransbin", default=None, type=int,
        help="The smallest sized bin to use for trans interactions (maximum resolution), in base pairs. [default: %(default)s]")
    parser.add_argument("-m", "--mid-binsize", dest="midbin", default=40000, type=int,
        help="The smallest sized bin to use for the entire chromosome(s), in base pairs. This is used to balance memory usage vs. speed. [default: %(default)s]")
    parser.add_argument("-d", "--datatype", dest="datatype", default="fend",
        help="Which corrections (if any) to apply to counts. [default: %(default)s]",
        choices=["raw", "fend", "distance", "enrichment"])
    add_silent_argument(parser)
    parser.add_argument(dest="project", type=str,
        help="The name of a HiFive HiC project file to pull data from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiC multi-resolution heatmaps to.")
    return

def add_hic_interval_subparser(subparsers):
    """Add command 'hic-interval' arguments to parser."""
    parser = subparsers.add_parser("hic-interval",
        help="HiFive Binning Function: Create a tabular interaction file containing data from a HiFive HiC project. Data are a genomic-interval format (chr1 start1 stop1 chr2 start2 stop2).")
    parser.add_argument("-c", "--chromosome", dest="chrom", default=None, required=True, type=str,
        help="The chromosome from which to pull interaction data from.")
    parser.add_argument("--chromosome2", dest="chrom2", default=None, required=False, type=str,
        help="The second chromosome from which to pull interaction data from if pulling trans data.")
    parser.add_argument("-s", "--start", dest="start", default=None, required=False, type=int,
        help="The start coordinate of the pulled region to return. (None indicates the first valid bin on the chromosome) [default %(default)s]")
    parser.add_argument("-e", "--stop", dest="stop", default=None, required=False, type=int,
        help="The stop coordinate of the pulled region to return. (None indicates the last valid bin on the chromosome) [default %(default)s]")
    parser.add_argument("--start2", dest="start2", default=None, required=False, type=int,
        help="The start coordinate of the second chromosome pulled region to return. (None indicates the first valid bin on the chromosome) [default %(default)s]")
    parser.add_argument("--stop2", dest="stop2", default=None, required=False, type=int,
        help="The stop coordinate of the second chromosome pulled region to return. (None indicates the last valid bin on the chromosome) [default %(default)s]")
    parser.add_argument("-b", "--binsize", dest="binsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into. [default: %(default)s]")
    parser.add_argument("-m", "--max-distance", dest="maxdist", default=None, type=int,
        help="The maximum interaction distance to return (None indicates no maximum). [default: %(default)s]")
    parser.add_argument("-d", "--data-type", dest="datatype", default="fend",
        help="Which corrections (if any) to apply to counts. [default: %(default)s]",
        choices=["raw", "fend", "distance", "enrichment", "expected"])
    parser.add_argument("-M", "--matrix", dest="matrix", default=False, action="store_true",
        help="Store output as a tab-separated matrix of values.")
    parser.add_argument("-y", "--dynamically-bin", dest="dynamic", default=False, action="store_true",
        help="Dynamically bin heatmap.")
    parser.add_argument("-x", "--expansion-binsize", dest="expbinsize", default=10000, type=int,
        help="The size of bins, in base pairs, to group data into for expanding under-populated bins. [default: %(default)s]")
    parser.add_argument("-f", "--minobservations", dest="minobs", default=20, type=int,
        help="The minimum number of observed reads in a bin for it to be considered valid. [default: %(default)s]")
    parser.add_argument("-a", "--search-distance", dest="search", default=0, type=int,
        help="The furthest distance from the bin minpoint to expand bounds. If set to zero, there is no limit on expansion distance. [default: %(default)s]")
    parser.add_argument("-v", "--remove-failed", dest="remove", default=False, action="store_true",
        help="If a non-zero 'search-distance' is given, it is possible for a bin not to meet the 'minobservations' criteria before stopping looking. If this occurs and 'remove-failed' is set, the observed and expected values for that bin are zero.")
    parser.add_argument("-i", "--image-file", dest="image", default=None, type=str,
        help="Save the data as an image to this file.")
    parser.add_argument("-p", "--pdf", dest="pdf", default=False, action="store_true",
        help="Format the image in PDF format. [default: %(default)s]")
    parser.add_argument("-r", "--rotate", dest="rotate", default=False, action="store_true",
        help="Rotate the plot 45 degrees (cis binned only). [default: %(default)s]")
    parser.add_argument("-t", "--ticks", dest="ticks", default=False, action="store_true",
        help="Add tick marks and labels to the plot (pdf format and binned only). [default: %(default)s]")
    parser.add_argument("-l", "--legend", dest="legend", default=False, action="store_true",
        help="Add color scale to the plot (pdf format only). [default: %(default)s]")
    parser.add_argument("-k", "--keyword", dest="keywords", default=[], type=str, action='append',
        help="Additional keyword arguments to pass to plotting function.")
    add_silent_argument(parser)
    parser.add_argument(dest="project", type=str,
        help="The name of a HiFive HiC project file to pull data from.")
    parser.add_argument(dest="output", type=str,
        help="The name of the file to write HiC interval to.")
    return

def add_quasar_subparser(subparsers):
    """Add command 'quasar' arguments to parser. This command is MPI-compatible."""
    parser = subparsers.add_parser("quasar", help="HiFive HiC QuASAR scoring function: Create a new QuASAR transformation file or use an existing one and find quality and/or replicate scores from it.")
    parser.add_argument("-p", "--hic", dest="hic", type=str, required=False, default=None,
        help="The name of the HiFive HiC project to use for QuASAR transformation.")
    parser.add_argument("-P", "--hic2", dest="hic2", type=str, required=False, default=None,
        help="The name of the second HiFive HiC project to use for QuASAR replicate scoring.")
    parser.add_argument("-Q", "--quasar2", dest="quasar2", type=str, required=False, default=None,
        help="The name of the second HiFive QuASAR file to write to and pull data from for QuASAR replicate scoring.")
    parser.add_argument("-o", "--report", dest="report", type=str, required=False, default=None,
        help="The name of the file to write all scoring results to. This can be either a txt or PDF file and format will be determined by the file suffix.")
    parser.add_argument("-c", "--chromosomes", dest="chroms", required=False, type=str, default=None,
        action='store', help="A comma-separated list of chromosomes to find QuASAR transformations for (None indicates all chromosomes). [default: %(default)s]")
    parser.add_argument("-r", "--resolutions", dest="resolutions", required=False, type=str,
        default="10000,40000,200000,1000000", action='store', help="A comma-separated list of resolutions to find QuASAR transformations for. [default: %(default)s]")
    parser.add_argument("-d", "--coverages", dest="coverages", required=False, type=str,
        default="0,1000000,2000000,4000000,8000000,16000000,32000000", action='store', help="A comma-separated list of coverages to find QuASAR transformations for. [default: %(default)s]")
    parser.add_argument("--scores_only", dest="scores_only", action='store_true',
        help="Report only scores in report, no additional analyses. [default: %(default)s]")
    parser.add_argument("--seed", dest="seed", default=None, type=int,
        help="The seed value for the random number generating function. [default: %(default)s]")
    add_silent_argument(parser)
    parser.add_argument(dest="quasar", type=str,
        help="The name of the HiFive QuASAR file to write to and pull data from.")
    return

def add_common_heatmap_arguments(parser):
    """Add common heatmap arguments to parse."""
    parser.add_argument("-t", "--trans", dest="trans", default=False, action="store_true",
        help="Calculate and include trans interactions in heatmaps. [default: %(default)s]")
    parser.add_argument("-c", "--chromosomes", dest="chroms", default=None, type=str,
        help="A comma-separated list of chromosomes to include in the heatmaps (None indicates all chromosomes). [default: %(default)s]")
    parser.add_argument("-f", "--minobservations", dest="minobs", default=20, type=int,
        help="The minimum number of observed reads in a bin for it to be considered valid. [default: %(default)s]")
    return

def add_silent_argument(parser):
    """Add silent argmuent to parser."""
    parser.add_argument("-q", "--quiet", dest="silent", required=False, default=False,
        action='store_true', help="Silence output messages. [default: %(default)s]")
    return


if __name__ == "__main__":
    if "mpi4py" in sys.modules.keys():
        try:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            num_procs = comm.Get_size()
        except:
            comm = None
            rank = 0
            num_procs = 1
    else:
        comm = None
        rank = 0
        num_procs = 1
    main()

