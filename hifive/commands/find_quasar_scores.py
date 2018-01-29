#!/usr/bin/env python

import sys

from ..quasar import Quasar
from ..hic import HiC

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_state()
    num_procs = comm.Get_size()
except:
    comm = None
    rank = 0
    num_procs = 1
    pass


def run(args):
    try:
        q1 = Quasar(args.quasar, 'a', silent=args.silent)
        if args.quasar2 is not None:
            q2 = Quasar(args.quasar2, 'a', silent=args.silent)
    except:
        print >> sys.stderr, ("Quasar file(s) are not of the correct filetype.\n"),
        return None
    if args.hic is None:
        if rank == 0:
            try:
                missing = False
                chroms = q1.storage['chromosomes'][...]
                coverages = q1.storage['coverages'][...]
                resolutions = q1.storage['resolutions'][...]
                for chrom in chroms:
                    for cov in coverages:
                        for res in resolutions:
                            if ("%s.%iC.%iR.invalid" % (chrom, cov, res) not in q1.storage.attrs and
                                ("valid.%s.%iC.%iR" % (chrom, cov, res) not in q1.storage or
                                "dist.%s.%iC.%iR" % (chrom, cov, res) not in q1.storage or
                                "corr.%s.%iC.%iR" % (chrom, cov, res) not in q1.storage)):
                                missing = True
                if missing:
                    print >> sys.stderr, ("Quasar file appears incomplete. Rerun with HiC project argument.\n"),
            except:
                print >> sys.stderr, ("Quasar file appears incomplete. Rerun with HiC project argument.\n"),
                missing = True
        else:
            missing = False
        if comm is not None:
            missing = comm.bcast(missing, root=0)
        if missing:
            return None

    if args.quasar2 is not None and args.hic2 is None:
        if rank == 0:
            try:
                missing = False
                chroms = q2.storage['chromosomes'][...]
                coverages = q2.storage['coverages'][...]
                resolutions = q2.storage['resolutions'][...]
                for chrom in chroms:
                    for cov in coverages:
                        for res in resolutions:
                            if ("%s.%iC.%iR.invalid" % (chrom, cov, res) not in q1.storage.attrs and
                                ("valid.%s.%iC.%iR" % (chrom, cov, res) not in q2.storage or
                                "dist.%s.%iC.%iR" % (chrom, cov, res) not in q2.storage or
                                "corr.%s.%iC.%iR" % (chrom, cov, res) not in q2.storage)):
                                missing = True
                if missing:
                    print >> sys.stderr, ("Second Quasar file appears incomplete. Rerun with HiC project argument.\n"),
            except:
                print >> sys.stderr, ("Second Quasar file appears incomplete. Rerun with HiC project argument.\n"),
                missing = True
        else:
            missing = False
        if comm is not None:
            missing = comm.bcast(missing, root=0)
        if missing:
            return None

    if args.hic2 is not None and args.quasar2 is None:
        if rank == 0:
            print >> sys.stderr, ("Specify a second Quasar file or don't specify a second HiC file.\n"),
        return None
    try:
        args.coverages = args.coverages.split(',')
        for i in range(len(args.coverages)):
            args.coverages[i] = int(args.coverages[i])
    except:
        if rank == 0:
            print >> sys.stderr, ("The coverages argument is not properly formattted. It should be a comma-separated list of integers.\n"),
        return None
    try:
        args.resolutions = args.resolutions.split(',')
        for i in range(len(args.resolutions)):
            args.resolutions[i] = int(args.resolutions[i])
    except:
        if rank == 0:
            print >> sys.stderr, ("The resolutions argument is not properly formattted. It should be a comma-separated list of integers.\n"),
        return None
    if args.chroms is None:
        args.chroms = []
    else:
        args.chroms = args.chroms.split(',')
    if args.hic is not None:
        try:
            hic1 = HiC(args.hic)
        except:
            if rank == 0:
                print >> sys.stderr, ("HiC file appears incomplete. Rerun with HiC project argument.\n"),
            return None
        q1.find_transformation(hic1, chroms=args.chroms, resolutions=args.resolutions,
                         coverages=args.coverages, seed=args.seed)
    qscores = q1.find_quality_scores(chroms=args.chroms)
    rscores = None
    q1.save()
    if args.quasar2 is not None:
        if args.hic2 is not None:
            try:
                hic2 = HiC(args.hic2)
            except:
                if rank == 0:
                    print >> sys.stderr, ("Second HiC file appears incomplete. Rerun with HiC project argument.\n"),
                return None
            q2.find_transformation(hic2, chroms=args.chroms, resolutions=args.resolutions,
                             coverages=args.coverages, seed=args.seed)
        q2.save()
        rscores = q1.find_replicate_scores(q2, chroms=args.chroms)
        q2.close()
    if args.report is not None:
        q1.print_report(args.report, qscores=qscores, rscores=rscores, scores_only=args.scores_only)
    q1.close()
