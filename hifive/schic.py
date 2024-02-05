#!/usr/bin/env python

"""A class for scHiC model preparation."""

import os
import sys

import numpy
import h5py


class scHiC(object):

    """
    This is the class for handling scHiC analysis.

    This class relies on :class:`Fend <hifive.fend.Fend>` and :class:`HiCData <hifive.hic_data.HiCData>` for genomic position and interaction count data. Use this class to perform filtering of interactions and 3D model preparation.

        .. note::
          This class is also available as hifive.scHiC

    When initialized, this class creates an h5dict in which to store all data associated with this object.

    :param filename: The file name of the h5dict. This should end with the suffix '.hdf5'
    :type filename: str.
    :param mode: The mode to open the h5dict with. This should be 'w' for creating or overwriting an h5dict with name given in filename.
    :type mode: str.
    :param silent: Indicates whether to print information about function execution for this object.
    :type silent: bool.
    :returns: :class:`HiC <hifive.hic.HiC>` class object.

    :attributes: * **file** (*str.*) - A string containing the name of the file passed during object creation for saving the object to.
                 * **silent** (*bool.*) - A boolean indicating whether to suppress all of the output messages.
                 * **history** (*str.*) - A string containing all of the commands executed on this object and their outcome.
                 * **normalization** (*str.*) - A string stating which type of normalization has been performed on this object. This starts with the value 'none'.
                 * **comm** (*class*) - A link to the MPI.COMM_WORLD class from the mpi4py package. If this package isn't present, this is set to 'None'.
                 * **rank** (*int.*) - The rank integer of this process, if running with mpi, otherwise set to zero.
                 * **num_procs** (*int.*) - The number of processes being executed in parallel. If mpi4py package is not present, this is set to one.

    In addition, many other attributes are initialized to the 'None' state.
    """

    def __init__(self, filename, mode='r', silent=False):
        """Create a HiC object."""
        self.file = os.path.abspath(filename)
        self.filetype = 'schic_project'
        self.silent = silent
        self.history = ''
        if mode != 'w':
            self.load()
        return None

    def __getitem__(self, key):
        """Dictionary-like lookup."""
        if key in self.__dict__:
            return self.__dict__[key]
        else:
            return None

    def __setitem__(self, key, value):
        """Dictionary-like value setting."""
        self.__dict__[key] = value
        return None

    def load_data(self, filename):
        """
        Load fend-pair counts and fend object from :class:`HiCData <hifive.hic_data.HiCData>` object.

        :param filename: Specifies the file name of the :class:`HiCData <hifive.hic_data.HiCData>` object to associate with this analysis.
        :type filename: str.
        :returns: None

        :Attributes: * **datafilename** (*str.*) - A string containing the relative path of the HiCData file.
                     * **fendfilename** (*str.*) - A string containing the relative path of the Fend file associated with the HiCData file.
                     * **fends** (*filestream*) - A filestream to the hdf5 Fragment file such that all saved Fend attributes can be accessed through this class attribute.
                     * **data** (*filestream*) - A filestream to the hdf5 FiveCData file such that all saved HiCData attributes can be accessed through this class attribute.
                     * **chr2int** (*dict.*) - A dictionary that converts chromosome names to chromosome indices.
                     * **filter** (*ndarray*) - A numpy array of type int32 and size N where N is the number of fends. This contains the inclusion status of each fend with a one indicating included and zero indicating excluded and is initialized with all fends included.

        When a HiCData object is associated with the project file, the 'history' attribute is updated with the history of the HiCData object.
        """
        self.history += "HiC.load_data(filename='%s') - " % (filename)
        filename = os.path.abspath(filename)
        # ensure data h5dict exists
        if not os.path.exists(filename):
            if not self.silent:
                print("Could not find %s. No data loaded." % (filename), file=sys.stderr)
            self.history += "Error: '%s' not found\n" % filename
            return None
        self.datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(filename)),
                                       os.path.dirname(self.file)), os.path.basename(filename))
        self.data = h5py.File(filename, 'r')
        self.history = self.data['/'].attrs['history'] + self.history
        fendfilename = self.data['/'].attrs['fendfilename']
        if fendfilename[:2] == './':
            fendfilename = fendfilename[2:]
        parent_count = fendfilename.count('../')
        fendfilename = '/'.join(os.path.abspath(filename).split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
        self.fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(fendfilename),
                                       os.path.dirname(self.file)), os.path.basename(fendfilename))
        # ensure fend h5dict exists
        if not os.path.exists(fendfilename):
            if not self.silent:
                print("Could not find %s." % (fendfilename), file=sys.stderr)
            self.history += "Error: '%s' not found\n" % fendfilename
            return None
        self.fends = h5py.File(fendfilename, 'r')
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.fends['chromosomes']):
            self.chr2int[chrom] = i
        # create arrays
        self.filter = numpy.ones(self.fends['fends'].shape[0], dtype=numpy.int32)
        self.cis_filter = numpy.ones(self.data['cis_data'].shape[0], dtype=numpy.int32)
        self.trans_filter = numpy.ones(self.data['trans_data'].shape[0], dtype=numpy.int32)
        self.history += "Succcess\n"
        return None

    def save(self, out_fname=None):
        """
        Save analysis parameters to h5dict.

        :param filename: Specifies the file name of the :class:`HiC <hifive.hic.HiC>` object to save this analysis to.
        :type filename: str.
        :returns: None
        """
        self.history.replace("'None'", "None")
        if not out_fname is None:
            original_file = os.path.abspath(self.file)
            out_fname = os.path.abspath(out_fname)
            if 'datafilename' in self.__dict__:
                datafilename = self.datafilename
                if datafilename[:2] == './':
                    datafilename = datafilename[2:]
                parent_count = datafilename.count('../')
                datafilename = '/'.join(original_file.split('/')[:-(1 + parent_count)] +
                                        datafilename.lstrip('/').split('/')[parent_count:])
                datafilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(datafilename)),
                                               os.path.dirname(out_fname)), os.path.basename(datafilename))
            if 'fendfilename' in self.__dict__:
                fendfilename = self.fendfilename
                if fendfilename[:2] == './':
                    fendfilename = fendfilename[2:]
                parent_count = fendfilename.count('../')
                fendfilename = '/'.join(original_file.split('/')[:-(1 + parent_count)] +
                                        fendfilename.lstrip('/').split('/')[parent_count:])
                fendfilename = "%s/%s" % (os.path.relpath(os.path.dirname(os.path.abspath(fendfilename)),
                                               os.path.dirname(out_fname)), os.path.basename(fendfilename))
        else:
            out_fname = self.file
            if 'datafilename' in self.__dict__:
                datafilename = self.datafilename
            if 'fendfilename' in self.__dict__:
                fendfilename = self.fendfilename
        datafile = h5py.File(out_fname, 'w')
        for key in self.__dict__.keys():
            if key in ['data', 'fends', 'file', 'chr2int', 'comm', 'rank', 'num_procs', 'silent']:
                continue
            elif self[key] is None:
                continue
            elif key == 'fendfilename':
                datafile.attrs[key] = fendfilename
            elif key == 'datafilename':
                datafile.attrs[key] = datafilename
            elif isinstance(self[key], numpy.ndarray):
                datafile.create_dataset(key, data=self[key])
            elif not isinstance(self[key], dict):
                datafile.attrs[key] = self[key]
        datafile.close()
        return None

    def load(self):
        """
        Load analysis parameters from h5dict specified at object creation and open h5dicts for associated :class:`HiCData <hifive.hic_data.HiCData>` and :class:`Fend <hifive.fend.Fend>` objects.

        Any call of this function will overwrite current object data with values from the last :func:`save` call.

        :returns: None
        """
        # set parameters to init state
        self.history = ''
        # load data hdf5 dict
        datafile = h5py.File(self.file, 'r')
        for key in datafile.keys():
            self[key] = numpy.copy(datafile[key])
        for key in datafile['/'].attrs.keys():
            self[key] = datafile['/'].attrs[key]
        # ensure data h5dict exists
        if 'datafilename' in self.__dict__:
            datafilename = self.datafilename
            if datafilename[:2] == './':
                datafilename = datafilename[2:]
            parent_count = datafilename.count('../')
            datafilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                datafilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(datafilename):
                if not self.silent:
                    print("Could not find %s. No data loaded." % (datafilename), file=sys.stderr)
            else:
                self.data = h5py.File(datafilename, 'r')
        # ensure fend h5dict exists
        if 'fendfilename' in self.__dict__:
            fendfilename = self.fendfilename
            if fendfilename[:2] == './':
                fendfilename = fendfilename[2:]
            parent_count = fendfilename.count('../')
            fendfilename = '/'.join(self.file.split('/')[:-(1 + parent_count)] +
                                fendfilename.lstrip('/').split('/')[parent_count:])
            if not os.path.exists(fendfilename):
                if not self.silent:
                    print("Could not find %s. No fends loaded." % (fendfilename), file=sys.stderr)
            else:
                self.fends = h5py.File(fendfilename, 'r')
        # create dictionary for converting chromosome names to indices
        self.chr2int = {}
        for i, chrom in enumerate(self.fends['chromosomes']):
            self.chr2int[chrom] = i
        datafile.close()
        return None

    def filter_interactions(self, minsize=21, maxsize=10000, excluded_chroms=['Y'], haploid=True, maxsep=40, maxsteps=3, resolution=1000000, invalid_fends=None):
        # remove user-specified invalid fends
        if invalid_fends is not None:
            for i in invalid_fends:
                self.filter[i] = 0

        # remove small fragments
        fends = self.fends['fends'][...]
        sizes = fends['stop'] - fends['start']
        prev_filt = numpy.sum(self.filter)
        self.filter[numpy.where((sizes < minsize) | (sizes > maxsize))[0]] = 0

        # remove excluded chromosomes
        chr_indices = self.fends['chr_indices'][...]
        for chrom in excluded_chroms:
            chrint = self.chr2int[chrom]
            self.filter[chr_indices[chrint]:chr_indices[chrint + 1]] = 0

        # load data
        cis_data = self.data['cis_data'][...].astype(numpy.int64)
        trans_data = self.data['trans_data'][...].astype(numpy.int64)

        # filter fends without interactions
        contacts = numpy.zeros(chr_indices[-1], dtype=numpy.int32)
        contacts += numpy.bincount(cis_data[:, 0], minlength=chr_indices[-1])
        contacts += numpy.bincount(cis_data[:, 1], minlength=chr_indices[-1])
        contacts += numpy.bincount(trans_data[:, 0], minlength=chr_indices[-1])
        contacts += numpy.bincount(trans_data[:, 1], minlength=chr_indices[-1])
        self.filter[numpy.where(contacts == 0)[0]] = 0

        # if haploid, make sure that fragments don't have multiple contacts
        if haploid:
            where = numpy.where(contacts > 2)[0]
            self.filter[where] = 0
            where = numpy.where(contacts == 2)[0]
            for i in where:
                where1 = numpy.where((cis_data[:, 0] == i) | (cis_data[:, 1] == i))[0]
                if where1.shape[0] == 2:
                    f1 = cis_data[where1[0], numpy.where(cis_data[where1[0], :2] != i)[0]]
                    f2 = cis_data[where1[1], numpy.where(cis_data[where1[1], :2] != i)[0]]
                    if abs(f2 - f1) > maxsep:
                        self.filter[i] = 0
                    else:
                        if abs(f1 - i) < abs(f2 - i):
                            self.cis_filter[where1[0]] = 0
                        else:
                            self.cis_filter[where1[1]] = 0
                elif where1.shape[0] == 0:
                    where1 = numpy.where((trans_data[:, 0] == i) | (trans_data[:, 1] == i))[0]
                    f1 = trans_data[where1[0], numpy.where(trans_data[where1[0], :2] != i)[0]]
                    f2 = trans_data[where1[1], numpy.where(trans_data[where1[1], :2] != i)[0]]
                    if abs(f2 - f1) > maxsep or fends['chr'][f1] != fends['chr'][f2]:
                        self.filter[i] = 0
                    else:
                        self.trans_filter[where1[numpy.random.randint(0, 2)]] = 0
                else:
                    self.filter[i] = 0

        if not self.silent:
            print("Filtering fends: %i of %i kept" % (numpy.sum(self.filter), self.filter.shape[0]), file=sys.stderr)

        # filter by connectivity graph
        self.cis_filter[numpy.where((self.filter[cis_data[:, 0]] == 0) | (self.filter[cis_data[:, 1]] == 0))] = 0
        self.trans_filter[numpy.where((self.filter[trans_data[:, 0]] == 0) | (self.filter[trans_data[:, 1]] == 0))] = 0
        mapping = fends['mid'] / resolution
        N = 0
        for i in range(chr_indices.shape[0] - 1):
            mapping[chr_indices[i]:chr_indices[i + 1]] += N
            N = mapping[chr_indices[i + 1] - 1] + 1
        cis_data[:, 0] = mapping[cis_data[:, 0]]
        cis_data[:, 1] = mapping[cis_data[:, 1]]
        trans_data[:, 0] = mapping[trans_data[:, 0]]
        trans_data[:, 1] = mapping[trans_data[:, 1]]
        cindices = cis_data[:, 0] * N + cis_data[:, 1]
        tindices = trans_data[:, 0] * N + trans_data[:, 1]
        connections = numpy.unique(numpy.r_[cindices, tindices])
        connections0 = connections // N
        connections1 = connections % N
        counts = numpy.bincount(numpy.searchsorted(connections, cindices), weights=(self.cis_filter * cis_data[:, 2]), minlength=connections.shape[0])
        counts += numpy.bincount(numpy.searchsorted(connections, tindices), weights=(self.trans_filter * trans_data[:, 2]), minlength=connections.shape[0])
        indices = numpy.triu_indices(N, 1)
        steps = numpy.zeros((N, N), dtype=numpy.int32)
        steps[connections1, connections0] = 1
        steps[connections0, connections1] = 1
        steps[numpy.arange(steps.shape[0]), numpy.arange(steps.shape[0])] = 0
        for i in range(1, maxsteps - 1):
            for j in range(indices[0].shape[0]):
                X = indices[0][j]
                Y = indices[1][j]
                if steps[X, Y] != 0:
                    continue
                temp = steps[X, :] + steps[Y, :]
                valid = numpy.where((steps[X, :] > 0) & (steps[Y, :] > 0))[0]
                if valid.shape[0] > 0:
                    steps[X, Y] = numpy.amin(temp[valid])
                    steps[Y, X] = steps[X, Y]
        valid = numpy.ones(connections1.shape[0], dtype=numpy.bool)
        for i in range(valid.shape[0]):
            if counts[i] > 1 or connections0[i] == connections1[i]:
                continue
            temp = steps[X, :] + steps[Y, :]
            tvalid = numpy.where((steps[X, :] > 0) & (steps[Y, :] > 0))[0]
            if tvalid.shape[0] == 0 or numpy.amin(temp[tvalid]) > maxsteps:
                fends1 = numpy.where(mapping == connections0[i])[0]
                fends2 = numpy.where(mapping == connections1[i])[0]
                where = numpy.where((cis_data[:, 0] == connections0[i]) & (cis_data[:, 1] == connections1[i]))[0]
                self.cis_filter[where] = 0
                where = numpy.where((trans_data[:, 0] == connections0[i]) & (trans_data[:, 1] == connections1[i]))[0]
                self.trans_filter[where] = 0
        if not self.silent:
            print("Filtering interactions: %i of %i cis kept, %i of %i trans kept" % (numpy.sum(self.cis_filter), self.cis_filter.shape[0], numpy.sum(self.trans_filter), self.trans_filter.shape[0]), file=sys.stderr)
        return None

    def run_simulation(self, resolution=1000000, tmpdir='./', seed=2001):
        RNG = numpy.random.RandomState(seed=seed)
        # bin data at selected resolution
        fends = self.fends['fends'][...]
        cis_data = self.data['cis_data'][...].astype(numpy.int64)
        cis_data[numpy.where((self.filter[cis_data[:, 0]] == 0) | (self.filter[cis_data[:, 1]] == 0))[0], 2] = 0
        cis_data[numpy.where(numpy.logical_not(self.cis_filter))[0], 2] = 0
        trans_data = self.data['trans_data'][...].astype(numpy.int64)
        trans_data[numpy.where((self.filter[trans_data[:, 0]] == 0) | (self.filter[trans_data[:, 1]] == 0))[0], 2] = 0
        trans_data[numpy.where(numpy.logical_not(self.trans_filter))[0], 2] = 0
        chr_indices = self.fends['chr_indices'][...]
        mapping = fends['mid'] // resolution
        N = 0
        bin_indices = numpy.zeros(chr_indices.shape[0], dtype=numpy.int32)
        for i in range(chr_indices.shape[0] - 1):
            mapping[chr_indices[i]:chr_indices[i + 1]] += N
            N = mapping[chr_indices[i + 1] - 1] + 1
            bin_indices[i + 1] = N
        cis_data[:, 0] = mapping[cis_data[:, 0]]
        cis_data[:, 1] = mapping[cis_data[:, 1]]
        trans_data[:, 0] = mapping[trans_data[:, 0]]
        trans_data[:, 1] = mapping[trans_data[:, 1]]
        cindices = cis_data[:, 0] * N + cis_data[:, 1]
        tindices = trans_data[:, 0] * N + trans_data[:, 1]
        connections = numpy.unique(numpy.r_[cindices, tindices])
        connections0 = connections // N
        connections1 = connections % N
        counts = numpy.bincount(numpy.searchsorted(connections, cindices), weights=cis_data[:, 2], minlength=connections.shape[0])
        counts += numpy.bincount(numpy.searchsorted(connections, tindices), weights=trans_data[:, 2], minlength=connections.shape[0])
        where = numpy.where(counts > 0)[0]
        connections0 = connections0[where]
        connections1 = connections1[where]
        counts = counts[where]

        # determine if we need to specifically include non-interacting contraints (no interactions within 250 Kb)
        constraints = numpy.zeros((N, N), dtype=numpy.int32)
        width = int(round(500000. / resolution) - 1) / 2
        if width >= 1:
            indices = numpy.triu_indices(N, 1)
            constraints[indices] = -1

        # add interaction constraints
        constraints[connections0, connections1] = 2
        constraints[connections1, connections0] = counts

        # filter out chromosomes with no trans contacts
        filt = numpy.ones(bin_indices.shape[0] - 1, dtype=numpy.bool)
        filt1 = numpy.ones(N, dtype=numpy.bool)
        M = N
        for i in range(bin_indices.shape[0] - 1):
            if (numpy.sum(constraints[:bin_indices[i], bin_indices[i]:bin_indices[i + 1]] == 2) +
                numpy.sum(constraints[bin_indices[i]:bin_indices[i + 1], bin_indices[i + 1]:] == 2) == 0):
                filt[i] = False
                filt1[bin_indices[i]:bin_indices[i + 1]] = False
                print("%i removed" % i, file=sys.stderr)
                M -= bin_indices[i + 1] - bin_indices[i]


        # if needed, add specific non-interacting constraints
        if width >= 1:
            w2 = width * 2
            for i in range(bin_indices.shape[0]):
                start = bin_indices[i]
                stop = bin_indices[i + 1]
                for j in range(start, stop - 2 * width - 1):
                    for k in range(j + 2 * width + 1, stop):
                        if numpy.sum(constraints[(k - w2):min(stop, k + w2 + 1), max(0, j - w2):(j + w2 + 1)]) == 0:
                            constraints[j, k] = 0

        # add adjacent constraints
        for i in range(bin_indices.shape[0] - 1):
            constraints[numpy.arange(bin_indices[i], bin_indices[i + 1] - 1), numpy.arange(bin_indices[i] + 1, bin_indices[i + 1])] = 1

        # intialize coordinates
        coords = numpy.zeros((N, 3), dtype=numpy.float32)
        angles = (2.0 * numpy.pi) / M * numpy.arange(M)
        offsets = RNG.uniform(0, 300.0, M)
        angle2 = RNG.uniform(0, 2.0 * numpy.pi, M)
        distances = 3000.0 + offsets * numpy.cos(angle2)
        where = numpy.where(filt1)[0]
        coords[where, 2] = offsets * numpy.sin(angle2)
        coords[where, 0] = distances * numpy.cos(angles)
        coords[where, 1] = distances * numpy.sin(angles)
        scale = 100.0
        size = 99999.0 / scale
        coords /= scale
        coords += size / 2

        # write PBD
        output = open('%s/out.pdb' % tmpdir, 'w')
        print("HEADER", file=output)
        print("TITLE  Chromatin", file=output)
        print("MODEL %s" % ("1".rjust(4)), file=output)
        for i in range(bin_indices.shape[0] - 1):
            if not filt[i]:
                continue
            for j in range(bin_indices[i], bin_indices[i + 1]):
                if j == bin_indices[i]:
                    mol = "DBB"
                elif j == bin_indices[i + 1] - 1:
                    mol = "DBE"
                else:
                    mol = "DBM"
                print("%s%s %s%s%s %s%s%s   %s%s%s%s%s          %s%s%s" % (
                    "ATOM".ljust(6),                       # [ 1 - 4 ] Entry type
                    str(j + 1).rjust(5),                   # [ 7 - 11] Atom serial number
                    str("ZZ").ljust(4),                    # [13 - 16] Atom name
                    " ",                                   # [  17   ] Alternate location indicator
                    str(mol).rjust(3),                     # [18 - 20] Residue name
                    chr(i + 97),                           # [  22   ] Chain identifier
                    str(j - bin_indices[i]).rjust(4),      # [23 - 26] Residue sequence number
                    " ",                                   # [  27   ] Code for insertion of residues
                    str("%0.3f" % coords[j, 0]).rjust(8),  # [31 - 38] Orthoginal X coordinate (in angstroms)
                    str("%0.3f" % coords[j, 1]).rjust(8),  # [39 - 46] Orthoginal X coordinate (in angstroms)
                    str("%0.3f" % coords[j, 2]).rjust(8),  # [47 - 54] Orthoginal X coordinate (in angstroms)
                    str("1.0").rjust(6),                   # [55 - 60] Occupancy
                    str("1.0").rjust(6),                   # [61 - 66] Temperature factor
                    str(" ").ljust(4),                     # [73 - 76] Segment identifier
                    str("ZZ").rjust(2),                    # [77 - 78] Element symbol
                    str(" ").rjust(2)),                    # [79 - 80] Charge on atom
                    file=output)
        print("ENDMDL", file=output)
        for i in range(bin_indices.shape[0] - 1):
            if not filt[i]:
                continue
            print("CONECT%s%s" % (str(bin_indices[i] + 1).rjust(5), str(bin_indices[i] + 2).rjust(5)), file=output)
            for j in range(bin_indices[i] + 1, bin_indices[i + 1] - 1):
                print("CONECT%s%s%s" % (str(j + 1).rjust(5), str(j).rjust(5), str(j + 2).rjust(5)), file=output)
            print("CONECT%s%s" % (str(bin_indices[i + 1]).rjust(5), str(bin_indices[i + 1] - 1).rjust(5)), file=output)
        print("END", file=output)
        output.close()

        # write GRO
        output = open('%s/out.gro' % tmpdir, 'w')
        vel = None#RNG.uniform(-0.1, 0.1, coords.shape)
        print("Chromatin", file=output)
        print("%i" % (M), file=output)
        for i in range(bin_indices.shape[0] - 1):
            if not filt[i]:
                continue
            for j in range(bin_indices[i], bin_indices[i + 1]):
                if j == bin_indices[i]:
                    mol = "DBB"
                elif j == bin_indices[i + 1] - 1:
                    mol = "DBE"
                else:
                    mol = "DBM"
                if vel is not None:
                    print("%s%s%s%s%s%s%s%s%s%s" % (
                        str(j - bin_indices[i]).rjust(5),     # [ 1 - 5 ] Residue sequence number
                        str(mol).ljust(5),                    # [ 6 - 10] Residue name
                        "ZZ".rjust(5),                        # [11 - 15] Atom name
                        str(j + 1).rjust(5),                  # [16 - 20] Atom serial number
                        ("%0.3f" % coords[j, 0]).rjust(8),    # [21 - 28] Orthoginal X coordinate(in angstroms)
                        ("%0.3f" % coords[j, 1]).rjust(8),    # [29 - 36] Orthoginal Y coordinate (in angstroms)
                        ("%0.3f" % coords[j, 2]).rjust(8),    # [37 - 44] Orthoginal Z coordinate (in angstroms)
                        ("%0.4f" % vel[j, 0]).rjust(8),       # [45 - 52] Velocity of X
                        ("%0.4f" % vel[j, 1]).rjust(8),       # [53 - 60] Velocity of Y
                        ("%0.4f" % vel[j, 2]).rjust(8)),      # [61 - 68] Velocity of Z
                        file=output)
                else:
                    print("%s%s%s%s%s%s%s" % (
                        str(j - bin_indices[i]).rjust(5),     # [ 1 - 5 ] Residue sequence number
                        str(mol).ljust(5),                    # [ 6 - 10] Residue name
                        "ZZ".rjust(5),                        # [11 - 15] Atom name
                        str(j + 1).rjust(5),                  # [16 - 20] Atom serial number
                        ("%0.3f" % coords[j, 0]).rjust(8),    # [21 - 28] Orthoginal X coordinate (in angstroms)
                        ("%0.3f" % coords[j, 1]).rjust(8),    # [29 - 36] Orthoginal Y coordinate (in angstroms)
                        ("%0.3f" % coords[j, 2]).rjust(8)),   # [37 - 44] Orthoginal Z coordinate (in angstroms)
                        file=output)
        print("%s%s%s" % (("%0.5f" % size).rjust(10), ("%0.5f" % size).rjust(10), ("%0.5f" % size).rjust(10)), file=output)
        output.close()

        # write TOP
        output = open('%s/out.top' % tmpdir, 'w')
        print('#include "charmm27.ff/forcefield.itp"', file=output)
        print("", file=output)
        print("[ moleculetype ]", file=output)
        print("Chromatin  0", file=output)
        print("", file=output)
        print("[ atoms ]", file=output)
        for i in range(bin_indices.shape[0] - 1):
            if not filt[i]:
                continue
            for j in range(bin_indices[i], bin_indices[i + 1]):
                if j == bin_indices[i]:
                    mol = "DBB"
                elif j == bin_indices[i + 1] - 1:
                    mol = "DBE"
                else:
                    mol = "DBM"
                print("%i  ZZ  %i  %s  ZZ  1  0.0 1.0" % (j + 1, j - bin_indices[i], mol), file=output)
        print("", file=output)
        params = numpy.array([[6, 10000, 10000, 1.0 / 25.0], [0, 150.0, 180.0, 1.0], [18.0, 120.0, 150.0, 1.0]])
        params /= scale
        print("[ distance_restraints ]", file=output)
        pos = 1
        for i in range(N - 1):
            if not filt1[i]:
                continue
            for j in range(i + 1, N):
                if not filt1[j]:
                    continue
                bond = constraints[i, j]
                if bond < 1: #< 0:
                    continue
                elif bond < 2:
                    d_low = params[bond, 0]
                    d_up1 = params[bond, 1]
                    d_up2 = params[bond, 2]
                    f = params[bond, 3]
                else:
                    if resolution < 250000:
                        d_low = params[bond, 0]
                        d_up1 = params[bond, 1]
                        d_up2 = params[bond, 2]
                        f = params[bond, 3]
                    else:
                        target = 150.0  / scale / (constraints[j, i] ** 2)
                        d_low = target * 0.8
                        d_up1 = target * 1.2
                        d_up2 = d_up1 + 30.0 / scale
                        f = params[bond, 3]
                print("%i  %i  %i  %i  %i  %0.2f  %0.2f  %0.2f  %0.2f" % (
                    i + 1,  # First atom j + 1,  # Second atom
                    1,      # type
                    pos,    # index
                    1,      # type'
                    d_low,  # lower distance bound
                    d_up1,  # upper distance bound 1
                    d_up2,  # upper distance bound 2
                    f),     # restraint scale factor
                    file=output)
                pos += 1
        print("", file=output)
        print("[ system ]", file=output)
        print("DNA as beads", file=output)
        print("", file=output)
        print("[ molecules ]", file=output)
        print("Chromatin  1", file=output)
        output.close()

        # write MDP
        temps = [[1000, 500]]
        tsteps = [4]
        msteps = [1000]
        temp_str = []
        time_str = []
        t = 0
        for i in range(len(temps)):
            step = numpy.linspace(temps[i][0], temps[i][1], tsteps[i] + 1)[:-1]
            for j in range(len(step)):
                temp_str += ["%i" % int(round(step[j])),"%i" % int(round(step[j]))]
                time_str += ["%i" % t, "%i" % (t + msteps[i] - 1)]
                t += msteps[i]
        output = open('%s/out.mdp' % tmpdir, 'w')
        print("%s = %s" % ("integrator".ljust(24), "steep"), file=output)
        print("%s = %i" % ("nsteps".ljust(24), 5000), file=output)
        print("%s = %s" % ("comm-mode".ljust(24), "Linear"), file=output)
        print("%s = %s" % ("annealing".ljust(24), "single"), file=output)
        print("%s = %i" % ("annealing-npoints".ljust(24), len(temp_str)), file=output)
        print("%s = %s" % ("annealing-time".ljust(24), ' '.join(time_str)), file=output)
        print("%s = %s" % ("annealing-temp".ljust(24), ' '.join(temp_str)), file=output)
        print("%s = %s" % ("constraints".ljust(24), "none"), file=output)
        print("%s = %s" % ("cutoff-scheme".ljust(24), "Verlet"), file=output)
        print("%s = %s" % ("disre".ljust(24), "simple"), file=output)
        print("%s = %i" % ("nstdisreout".ljust(24), 0), file=output)
        output.close()


