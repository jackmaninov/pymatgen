# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

import os
import re
import itertools
import warnings
import logging
import math

import six
from six import string_types
import fortranformat as ff
from copy import deepcopy

from monty.io import zopen
from monty.os.path import zpath
from monty.json import MontyDecoder

from enum import Enum
from tabulate import tabulate

import scipy.constants as const

from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.util.string import str_delimited
from pymatgen.util.io_utils import clean_lines
from monty.json import MSONable

"""
Classes for reading/manipulating/writing WIEN2k input files. Initial focus is
on case.struct, case.inc, case.inm, case.inso, case.int, case.inorb, case.klist
"""

logger = logging.getLogger(__name__)

BOHR = const.physical_constants['atomic unit of length'][0] * 10 ** 10

class Struct(MSONable):
    """
    case.struct reader/writer
    """

    def __init__(self, structure, comment=None, relativistic=True, symops=None, lattice=None):
        if structure.is_ordered:
            site_properties = {}
            self.comment = structure.formula if comment is None else comment
            self.relativistic = relativistic
            self.structure = structure.copy(site_properties=site_properties)
            self.symops = symops
            self.lattice = lattice
        else:
            raise ValueError("Structure with partial occupancies cannot be "
                             "converted into a WIEN2k struct!")

    @property
    def natoms(self):
        """
        Number of ATOM sites #TODO understand this
        :return:
        """
        syms = [site.specie.symbol for site in self.structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    @property
    def site_symbols(self):
        """
        Sequence of symbols associated with the Struct. Similar to 6th line in
        vasp 5+ POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [a[0] for a in itertools.groupby(syms)]

    def __setattr__(self, name, value):
        super(Struct, self).__setattr__(name, value)  # how is this not recursive?

    @staticmethod
    def from_file(filename):
        """
        read from a case.struct file
        :param filename: case.struct name (string)
        :return: Struct
        """

        with zopen(filename, "rt") as f:
            return Struct.from_string(f.read())

    @staticmethod
    def from_string(data):
        """

        :param data: string containing case.struct
        :return: Struct
        """

        lattype = {'P': 1, 'F': 2, 'B': 3, 'CXY': 4, 'CYZ': 5, 'CXZ': 6, 'R': 7, 'H': 8}
        rela = {'RELA': True, 'NREL': False}

        chunks = re.split("\n\s*\n", data.rstrip(), flags=re.MULTILINE)
        try:
            if chunks[0] == "":
                chunks.pop(0)
                chunks[0] = "\n" + chunks[0]
        except IndexError:
            raise ValueError("Empty STRUCT")

        lines = tuple(chunks[0].split("\n"))

        # line 1:
        comment = lines[0].rstrip()

        # line 2:
        reader2 = ff.FortranRecordReader('(A4,23X,I3)')
        try:
            structtype = lattype.get(reader2.read(lines[1])[0].rstrip())
        except IndexError:
            raise ValueError("Unknown lattice type")
        try:
            natoms = reader2.read(lines[1])[1]
        except ValueError:
            raise ValueError("Unknown NAT")

        # line 3:
        reader3 = ff.FortranRecordReader('(13X,A4)')
        try:
            relativistic = rela.get(reader3.read(lines[2])[0])
        except ValueError:
            raise ValueError("Unclear relativistic mode in STRUCT")

        # line 4:
        reader4 = ff.FortranRecordReader('(6F10.6)')
        try:
            [a, b, c, alpha, beta, gamma] = reader4.read(lines[3])
            assert (a is not None or b is not None or c is not None), 'Missing unit cell parameter'
            if alpha is None:
                alpha = 90.0
            if beta is None:
                beta = 90.0
            if gamma is None:
                gamma = 90.0
        except ValueError:
            raise ValueError("Error reading unit cell parameters")

        a *= BOHR
        b *= BOHR
        c *= BOHR

        # construct lattice
        if structtype is lattype['P']:  # primative (not hexagonal)
            lattice = Lattice.from_parameters(a, b, c,
                                              alpha, beta, gamma)
        elif structtype is lattype['F']:  # face-centered
            lattice = Lattice([[a / 2, b / 2, 0],
                               [a / 2, 0, c / 2],
                               [0, b / 2, c / 2]])
        elif structtype is lattype['B']:  # body-centered
            lattice = Lattice([[a / 2, -b / 2, c / 2],
                               [a / 2, b / 2, -c / 2],
                               [-a / 2, b / 2, c / 2]])
        elif structtype is lattype['CXY']:  # C-base-centered (orthorhombic only)
            lattice = Lattice([[a / 2, -b / 2, 0],
                               [a / 2, b / 2, 0],
                               [0, 0, c]])
        elif structtype is lattype['CYZ']:  # A-base-centered (orthorhombic only)
            lattice = Lattice([[a, 0, 0],
                               [0, -b / 2, c / 2],
                               [0, b / 2, c / 2]])  # B-base-centered (orthorh. and monoclinic symmetry)
        elif structtype is lattype['CXZ']:
            gamma_r = math.radians(gamma)
            lattice = Lattice([[a * math.sin(gamma_r) / 2, a * math.cos(gamma_r) / 2, -c],
                               [0, b, 0],
                               [a * math.sin(gamma_r) / 2, a * math.cos(gamma_r) / 2, c / 2]])
        # TODO implement hex/rhomb coordinates
        elif structtype is lattype['R']:  # WIEN2k rhombohedral setting
            lattice = Lattice([[a / (math.sqrt(3) / 2), -a / 2, c / 3],
                               [a / (math.sqrt(3) / 2), a / 2, c / 3],
                               [-a / math.sqrt(3), 0, c / 3]])
        elif structtype is lattype['H']:  # WIEN2k hexagonal setting
            lattice = Lattice([[math.sqrt(3) * a / 2, -a / 2, 0],
                               [0, a, 0],
                               [0, 0, c]])
        else:
            raise ValueError('Unknown STRUCT type')

        # line 5 (repeats for natoms)

        reader5 = ff.FortranRecordReader('(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8)')
        reader6 = ff.FortranRecordReader('(15X,I2,17X,I2)')
        reader7 = ff.FortranRecordReader('(A10,5X,I5,5X,F10.8,5X,F10.5,5X,F10.5)')
        reader8 = ff.FortranRecordReader('(20X,3F10.7)')
        reader11 = ff.FortranRecordReader('(I4)')
        reader12 = ff.FortranRecordReader('(3I2,F10.7)')
        # reader15 = ff.FortranRecordReader('(I8)')  #not used

        species = []
        coords = []
        sitenames = []
        isplits = []
        npts = []
        r0s = []
        rmts = []
        zs = []
        altcoordss = []
        rotlocs = []

        currentline = 4

        for i in range(natoms):
            # site line 1 (ex: ATOM    1: X=...)
            try:
                [sitenumber, x, y, z] = reader5.read(lines[currentline])
                assert sitenumber == i + 1
            except ValueError:
                raise ValueError('Error reading coordinates of Atom ', i + 1)
            currentline += 1
            # site line 2 (ex: MULT= 1 ISPLIT = 2)
            try:
                [multiplicity, isplit] = reader6.read(lines[currentline])
            except ValueError:
                raise ValueError('Bad multiplicity/ISPLIT on Atom ', i + 1)
            currentline += 1
            altcoords = []
            for j in range(multiplicity - 1):  # Read alternate site coordinates
                try:
                    altcoords.append(reader5.read(lines[currentline])[1:])
                except ValueError:
                    raise ValueError('Improper multiplicity on Atom ', i + 1)
                currentline += 1
            try:
                [atomname, npt, r0, rmt, atomic_z] = reader7.read(lines[currentline])
            except ValueError:
                raise ValueError('Bad site parameters on Atom ', i + 1)

            currentline += 1
            rotloc = []
            for rotlocline in lines[currentline:currentline + 3]:
                try:
                    rotloc.append(reader8.read(rotlocline))
                except ValueError:
                    raise ValueError('Error reading ROTLOC on site ', i + 1)
            currentline += 3

            for site in [[float(x), float(y), float(z)]] + altcoords:
                species.append(Element.from_Z(int(atomic_z)))  # NOTE: WIEN2k supports non-integer Z!
                coords.append(site)
                sitenames.append(atomname.strip())
                isplits.append(isplit)
                npts.append(npt)
                r0s.append(r0)
                rmts.append(rmt)
                zs.append(atomic_z)
                altcoordss.append(altcoords)
                rotlocs.append(rotloc)

        [nsym] = reader11.read(lines[currentline])
        currentline += 1

        symops = []

        for i in range(nsym):
            symop = []
            for symopline in lines[currentline:currentline + 3]:
                symop.append(reader12.read(symopline))
            currentline += 3
            symops.append(symop)
            # symopindex = reader15.read(lines[currentline]) #not necessary to store this
            currentline += 1
            # assert symopindex = i+1,'Misordered symmetry operations'

        # build structure

        site_properties = {'sitename': sitenames, 'ISPLIT': isplits, 'NPT': npts, 'R0': r0s,
                           'RMT': rmts, 'Z': zs, 'ROTLOC': rotlocs}
        struct = Structure(lattice, species, coords, to_unit_cell=False,
                           validate_proximity=False, coords_are_cartesian=False,
                           site_properties=site_properties)

        return Struct(struct, comment=comment, relativistic=relativistic, symops=symops,
                      lattice=lattice)

    def __repr__(self):
        """
        :return: String representation of Struct file
        """
        return self.get_string()

    def __str__(self):
        """
        :return: String representation of Struct file
        """
        return self.get_string()

    def get_string(self):
        """

        :return: string representation of the struct
        """

        latt = self.structure.lattice
        lines = [self.comment]

        # line 2
        # for now, output all lattices as P
        writer2 = ff.FortranRecordWriter('(A4,A23,I3)')
        lines.append(writer2.write(['P   ', "LATTICE,NONEQUIV.ATOMS:",
                                    self.structure.num_sites]))

        # line 3
        lines.append("MODE OF CALC=" + "RELA" if self.relativistic else "NREL")

        # line 4
        writer4 = ff.FortranRecordWriter('(6F10.6)')
        lines.append(writer4.write(list(x / BOHR for x in latt.abc) + list(latt.angles)))

        # line 5 - 7
        writer5 = ff.FortranRecordWriter('(A4,I4,A4,F10.8,A3,F10.8,A3,F10.8)')
        writer6 = ff.FortranRecordWriter('(A15,I2,A17,I2)')
        writer7 = ff.FortranRecordWriter('(A10,A5,I5,A5,F10.8,A5,F10.5,A5,F10.5)')
        writer8 = ff.FortranRecordWriter('(A20,3F10.7)')

        for i, site in enumerate(self.structure):
            coords = site.frac_coords
            lines.append(writer5.write(['ATOM', i, ': X=', coords[0], ' Y=',
                                        coords[1], ' Z=', coords[2]]))
            lines.append(writer6.write(['MULT=', 1,
                                        'ISPLIT=', self.structure.site_properties['ISPLIT'][i]]))
            lines.append(writer7.write([str(self.structure.site_properties['sitename'][i]).ljust(10),
                                        ' NPT=', self.structure.site_properties['NPT'][i],
                                        '  R0=', self.structure.site_properties['R0'][i],
                                        ' RMT=', self.structure.site_properties['RMT'][i],
                                        '   Z:', site.Z]))
            rotloc = self.structure.site_properties['ROTLOC'][i]
            lines.append(writer8.write(['LOCAL ROT MATRIX:'.ljust(20)] + rotloc[0]))
            lines.append(writer8.write([''] + rotloc[1]))
            lines.append(writer8.write([''] + rotloc[2]))

        # lines 11 - 15

        writer11 = ff.FortranRecordWriter('(I4,A35)')
        writer12 = ff.FortranRecordWriter('(3I2,F10.7)')
        writer15 = ff.FortranRecordWriter('(I8)')


        lines.append(writer11.write([len(self.symops) if self.symops else 0, 'NUMBER OF SYMMETRY OPERATIONS'.rjust(35)]))
        for (i, symop) in enumerate(self.symops or []):
            for line in symop:
                lines.append(writer12.write(line))
            lines.append(writer15.write([i]))

        return "\n".join(lines)

    def write_file(self, filename):
        """
        Write Struct to a file.
        :param filename:  filename to write Struct to
        """
        with zopen(filename, "wt", newline="\n") as f:
            f.write(self.get_string())


class KListSupportedModes(Enum):
    Automatic = 0
    Gamma = 1
    Monkhorst = 2

    # TODO factor this out, since the paradigm doesn't apply

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s):
        c = s.lower()[0]
        for m in KListSupportedModes:
            if m.name.lower()[0] == c:
                return m
        raise ValueError("Can't interprete Kpoint mode %s" % s)


class Kpoints(MSONable):
    """
    case.klist reader/writer
    """
    supported_modes = KListSupportedModes

    def __init__(self, comment="Default gamma", num_kpts=0,
                 style=supported_modes.Gamma,
                 kpts=((1, 1, 1),), kpts_shift=(0, 0, 0),
                 kpts_weights=None, coord_type=None, labels=None,
                 tet_number=0, tet_weight=0):
        """
        WIEN2k Kpoint class constructor. Note that the file should be saved as case.klist!

        Args:
            comment (str): String comment for Kpoints
            num_kpts: number of k-points to automatically generate. If
                set to zero, explicit cell divisions necessary.
            style: Modes of generating k-points. Use one of the
                Kpoints.supported_modes enum types.
            kpts (2D array): 2D array of kpoints.
            kpts_shift (3x1 aray): Shift for Kpoints.
            coord_type: in line-mode, this variable specifies whether the
                kpoints were given in Cartesian or Reciprocal coordinates.
            labels: In line-mode, this should provide a list of labels for¬
                each kpt. It is optional in explicit kpoint mode as comments for¬
                k-points.¬
            tet_number: For explicit kpoints, specifies the number of¬
                tetrahedrons for the tetrahedron method.¬
            tet_weight: For explicit kpoints, specifies the weight for each¬
                tetrahedron for the tetrahedron method.¬

        The default behavior of the construtor is for a Gamma centered,
        1x1x1 k-point grid with no shift.
        """

        self.comment = comment
        self.num_kpts = num_kpts
        self.kpts = kpts
        self.style = style
        self._style = None
        self.coord_type = coord_type
        self.kpts_weights = kpts_weights
        self.kpts_shift = kpts_shift
        self.labels = labels
        self.tet_number = tet_number
        self.tet_weight = tet_weight

    @property
    def style(self):
        return self._style

    @style.setter
    def style(self, style):
        if isinstance(style, six.string_types):
            style = Kpoints.supported_modes.from_string(style)

        if style in (Kpoints.supported_modes.Automatic,
                     Kpoints.supported_modes.Gamma,
                     Kpoints.supported_modes.Monkhorst) and len(self.kpts) > 1:
            raise ValueError("For fully automatic or automatic gamma or monk "
                             "kpoints, only a single line for the number of "
                             "divisions is allowed.")

        self._style = style

    @staticmethod
    def automatic(subdivisions):
        """
        Convenient static constructor for a fully automatic Kpoint grid, with
        gamma centered Monkhorst-Pack grids and the number of subdivisions
        along each reciprocal lattice vector determined by the default WIEN2k
        kgen behaviour.

        Args:
            subdivisions: Parameter determining number of subdivisions along
                each reciprocal lattice vector.

        Returns:
            Kpoints object
        """
        return Kpoints("Fully automatic kpoint scheme", 0,
                       style=Kpoints.supported_modes.Automatic,
                       kpts=[[subdivisions]])

    @staticmethod
    def gamma_automatic(kpts=(1, 1, 1), shift=(0, 0, 0)):
        """
        Convenient static constructor for an automatic Gamma centered Kpoint
        grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice
                vectors. Defaults to (1,1,1)
            shift: Shift to be applied to the kpoints. Defaults to (0,0,0).

        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0,
                       Kpoints.supported_modes.Gamma, kpts=[kpts],
                       kpts_shift=shift)

    @staticmethod
    def monkhorst_automatic(kpts=(2, 2, 2), shift=(0, 0, 0)):
        """
        Convenient static constructor for an automatic Monkhorst pack Kpoint
        grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice
                vectors. Defaults to (2,2,2)
            shift: Shift to be applied to the kpoints. Defaults to (0,0,0).

        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0,
                       Kpoints.supported_modes.Monkhorst, kpts=[kpts],
                       kpts_shift=shift)

    @staticmethod
    def monkhorst_automatic(kpts=(2, 2, 2), shift=(0, 0, 0)):
        """
        Convenient static constructor for an automatic Monkhorst pack Kpoint
        grid.

        Args:
            kpts: Subdivisions N_1, N_2 and N_3 along reciprocal lattice
                vectors. Defaults to (2,2,2)
            shift: Shift to be applied to the kpoints. Defaults to (0,0,0).

        Returns:
            Kpoints object
        """
        return Kpoints("Automatic kpoint scheme", 0,
                       Kpoints.supported_modes.Monkhorst, kpts=[kpts],
                       kpts_shift=shift)

    @staticmethod
    def automatic_density(structure, kppa, force_gamma=False):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density. Uses Gamma centered meshes for hexagonal cells and
        Monkhorst-Pack grids otherwise.

        Algorithm:
            Uses a simple approach scaling the number of divisions along each
            reciprocal lattice vector proportional to its length.

        Args:
            structure (Structure): Input structure
            kppa (int): Grid density
            force_gamma (bool): Force a gamma centered mesh (default is to
                use gamma only for hexagonal cells or odd meshes)

        Returns:
            Kpoints
        """
        comment = "pymatgen 4.7.6+ generated KPOINTS with grid density = " + \
                  "%.0f / atom" % kppa
        if math.fabs((math.floor(kppa ** (1 / 3) + 0.5)) ** 3 - kppa) < 1:
            kppa += kppa * 0.01
        latt = structure.lattice
        lengths = latt.abc
        ngrid = kppa / structure.num_sites
        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)

        num_div = [int(math.floor(max(mult / l, 1))) for l in lengths]

        is_hexagonal = latt.is_hexagonal()

        has_odd = any([i % 2 == 1 for i in num_div])
        if has_odd or is_hexagonal or force_gamma:
            style = Kpoints.supported_modes.Gamma
        else:
            style = Kpoints.supported_modes.Monkhorst

        return Kpoints(comment, 0, style, [num_div], [0, 0, 0])

    @staticmethod
    def automatic_density_by_vol(structure, kppvol, force_gamma=False):
        """
        Returns an automatic Kpoint object based on a structure and a kpoint
        density per inverse Angstrom of reciprocal cell.

        Algorithm:
            Same as automatic_density()

        Args:
            structure (Structure): Input structure
            kppvol (int): Grid density per Angstrom^(-3) of reciprocal cell
            force_gamma (bool): Force a gamma centered mesh

        Returns:
            Kpoints
        """
        vol = structure.lattice.reciprocal_lattice.volume
        kppa = kppvol * vol * structure.num_sites
        return Kpoints.automatic_density(structure, kppa,
                                         force_gamma=force_gamma)

    @staticmethod
    def from_file(filename):
        """
        Reads a Kpoints object from a KPOINTS file.

        Args:
            filename (str): filename to read from.

        Returns:
            Kpoints object
        """
        with zopen(filename, "rt") as f:
            return Kpoints.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads a Kpoints object from a case.klist string.

        Args:
            string (str): case.klist string.

        Returns:
            Kpoints object
        """

        reader = ff.FortranRecordReader('(I10,4I10,f5.1)')
        reader1 = ff.FortranRecordReader('(I10,4I10,3f5.1,2x,i8,\' k, div: (\',3i3,\')\')')
        kpts = []
        kpts_weights = []

        lines = [line.strip() for line in string.splitlines()]

        for line_number, line in enumerate(lines):
            if line == "END":
                break
            if line_number == 1:
                try:
                    [kpt, k_x, k_y, k_z, k_weight, i_div, junk1, junk2, num_kpts0, n_divx, n_divy,
                     ndivz] = reader1.read(
                        line)
                    style = Kpoints.supported_modes.Gamma if k_x == k_y == k_z == 0 else Kpoints.supported_modes.Monkhorst
                except ValueError:
                    raise ValueError('Error reading first k-point', line)
            else:
                try:
                    [kpt, k_x, k_y, k_z, k_weight, i_div] = reader.read(line)
                except ValueError:
                    raise ValueError('Error reading k-point:', line)

            kpts.append([float(k_x / i_div), float(k_y / i_div), float(k_z / i_div)])
            kpts_weights.append(float(k_weight))

        return Kpoints.gamma_automatic(kpts) if style == Kpoints.supported_modes.Gamma \
            else Kpoints.monkhorst_automatic(kpts)


# class In0(MSONable):
#
# class In1(MSONable):
#
# class In2(MSONable):
#
# class Inc(MSONable):
#
# class Inm(MSONable):
#
# class Inso(MSONable):
#
# class Inst(MSONable):
#
# class Kgen(MSONable):


class WIEN2kInput(dict, MSONable):
    """
    Class to contain a set of WIEN2k input objects, corresponding to a 'case'.
    """

    def __init__(self, struct, kpoints, optional_files=None,
                 **kwargs):
        super
        self.update({'struct': struct,
                     'kpoints': kpoints})
        if optional_files is not None:
            self.update(optional_files)

    def __str__(self):
        output = []
        for k, v in self.items():
            output.append(k)
            output.append(str(v))
            output.append("")
        return "\n".join(output)

    def as_dict(self):
        d = {k: v.as_dict() for k, v in self.items()}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        dec = MontyDecoder()
        sub_d = {"optional_files": {}}
        for k, v in d.items():
            if k in ["struct", "kpoints"]:
                sub_d[k.lower()] = dec.process_decoded(v)
            elif k not in ["@module", "@class"]:
                sub_d["optional_files"][k] = dec.process_decoded(v)
        return cls(**sub_d)

    def write_input(self, output_dir=".", make_dir_if_not_present=True):
        """
        Write a WIEN2k case

        :param output_dir (str): Director to write to. Defaults to current
            directory (".").
        :param make_dir_if_not_present (bool): Creat the directory if not present.
            Defaults to True.
        :return: none
        """

        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.items():
            with zopen(os.path.join(output_dir, k), "wt") as f:
                f.write(v.__str__())

    @staticmethod
    def from_directory(input_dir, casename=None, optional_files=None):
        """
        Read in a WIEN2k case from a directory. Note that only case.struct, .in0,
            .in1(c), .in2(c), .inc, .inm, .inso, .inst, .kgen, .klist, are read
            unless optional_filenames is specified.

        TODO: figure out what to read and progress through this
        :param input_dir (str): case directory to read
        :param optional_files (dict): Optional files to read in as well as a dict
            of {filename: Object type}. Object type must have a static method from_file
        :return: WIEN2kInput object
        """

        sub_d = {}
        if casename is None:
            case = os.path.basename(input_dir)
        else:
            case = casename
        for fname, ftype in [(case + ".struct", Struct), (case + ".in0", In0),
                             (case + ".in1", In1), (case + ".in2", In2),
                             (case + ".inc", Inc), (case + ".inm", Inm),
                             (case + ".inso", Inso), (case + ".inst", Inst),
                             (case + ".kgen", Kgen), (case + ".klist", Kpoints)]:
            fullzpath = zpath(os.path.join(input_dir, fname))
            sub_d[fname.lower()] = ftype.from_file(fullzpath)
        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = \
                    ftype.from_file(os.path.join(input_dir, fname))
        return WIEN2kInput(**sub_d)


VALID_TELNES3_TAGS = ("VERBOSITY", "ATOMS", "DETECTOR POSITION", "MODUS", "SPLIT",
                      "BRANCHING RATIO", "NONRELATIVISTIC", "RELATIVISTIC", "INITIALIZATION",
                      "QGRID", "ORIENTATION SENSITIVE", "SELECTION RULE", "LSELECTION RULE",
                      "EXTEND POTENTIALS", "FERMI ENERGY", "CORE WAVEFUNCTION",
                      "FINAL STATE WAVEFUNCTION", "NOHEADERS", "DOSONLY", "NBTOT", "END")


class TelnesTags(dict):
    """
    TELNES3 optional control parameters

    """

    def __init__(self, params=None):
        """
        :param params: a dictionary of input parameters
        """
        super(TelnesTags, self).__init__()
        if params:
            self.update(params)

    def __setitem__(self, key, val):
        """
        Add parameter-val pair. Warns is parameter is not in the list of valid TELNES tags.
        Also strips parameter and val.
        :param key: dict key value
        :param val: value to be associated with key in dict
        :return:
        """
        if key.strip().upper() not in VALID_TELNES3_TAGS:
            warnings.warn(key.strip() + " not in VALID_TENLES3_TAGS list")
        super(TelnesTags, self).__setitem__(key.strip(),
                                            TelnesTags.proc_val(key.strip(), val.strip())
                                            if isinstance(val, string_types) else val)

    def as_dict(self):
        """
        Dict representation

        :return: dict of paremeters from telnestags object
        """

        tags_dict = dict(self)
        tags_dict['@module'] = self.__class__.__module__
        tags_dict['@class'] = self.__class__.__name__
        return tags_dict

    @staticmethod
    def from_dict(d):
        """
        Creates Tags object from a dictionary.

        Args:
            d: Dict of feff parameters and values.

        Returns:
            Tags object
        """
        i = TelnesTags()
        for k, v in d.items():
            if k not in ("@module", "@class"):
                i[k] = v
        return i

    def get_string(self, sort_keys=False, pretty=False):
        """
        Returns a string representation of the TelnesTags.  The reason why this
        method is different from the __str__ method is to provide options
        for pretty printing.

        Args:
            sort_keys: Set to True to sort the Telnes parameters alphabetically.
                Defaults to False.
            pretty: Set to True for pretty aligned output. Defaults to False.

        Returns:
            String representation of TelnesTags.
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)
        lines = []
        for k in keys:
            lines.append([k, self._stringify_val(self[k])])
        if pretty:
            return tabulate(lines)
        else:
            return str_delimited(lines, None, " ")

    @staticmethod
    def _stringify_val(val):
        """
        Convert the given value to string.
        """
        if isinstance(val, list):
            if isinstance(val[0], list): #recurse into lists of lists
                return "\n".join([TelnesTags._stringify_val(i) for i in val])
            elif isinstance(val[0], float):
                return " ".join(['{:.3f}'.format(i) for i in val])
            else:
                return " ".join([str(i) for i in val])
        if isinstance(val, float):
            return '{:.3f}'.format(val)
        else:
            return str(val)

    def __str__(self):
        return self.get_string()

    @staticmethod
    def fromstring(input):
        """
        Read the optional tag components of a case.innes file
        :param input: input string
        :return: TelnesTags object
        """
        params = {}

        lines = list(clean_lines(input.split('\n')))
        for i, line in enumerate(lines):
            key = line.strip()
            if key in VALID_TELNES3_TAGS:
                if key == "END":
                    break
                elif key == "INITIALIZATION":
                    params[key] = [lines[i + 1].strip().split(), lines[i + 2].strip().split()]
                elif key == "QGRID":
                    params[key] = [lines[i+1]]
                    if params[key][0].upper() == "L":
                        params[key] = ["L", float(lines[i + 2].strip())]
                elif key not in ("DOSONLY", "NONRELATIVISTIC", "NOHEADERS"):
                    params[key] = lines[i + 1]
                elif key in ("DOSONLY", "NONRELATIVISTIC", "NOHEADERS"):
                    params[key] = True
                else:
                    warnings.warn(key + " parameter invalid")

        return TelnesTags(params)

    @staticmethod
    def proc_val(key, val):
        """
        Static helper method to convert Telnes parameters to proper types
        :param key: Telnes parameter key
        :param val: Actual value of Telnes parameter
        :return: properly formatted parameter
        """

        list_type_keys = ("ATOMS", "DETECTOR POSITION", "ORIENTATION SENSITIVE", "EXTEND POTENTIALS",
                          "QGRID")
        boolean_type_keys = ()
        float_type_keys = ("SPLIT", "BRANCHING RATIO", "FERMI ENERGY")
        int_type_keys = ("VERBOSITY", "RELATIVISTIC")

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key.lower() == 'cif':
                m = re.search(r"\w+.cif", val)
                return m.group(0)

            if key in list_type_keys:
                output = list()
                toks = re.split(r"\s+", val)

                for tok in toks:
                    m = re.match(r"(\d+)\*([\d\.\-\+]+)", tok)
                    if m:
                        output.extend([smart_int_or_float(m.group(2))] *
                                      int(m.group(1)))
                    else:
                        output.append(smart_int_or_float(tok))
                return output
            if key in boolean_type_keys:
                m = re.search(r"^\W+([TtFf])", val)
                if m:
                    if m.group(1) == "T" or m.group(1) == "t":
                        return True
                    else:
                        return False
                raise ValueError(key + " should be a boolean type!")

            if key in float_type_keys:
                return float(val)
            if key in int_type_keys:
                return int(val)

        except ValueError:
            return val

        return val


class Innes(MSONable):
    pass

    def __init__(self, config_dict=None, comment="EELS spectrum", absorbing_atom=1, edge='K',
                 e_onset=535.0, e_beam=200.0, e_start=0.0, e_final=35.0, e_step=0.05,
                 collection_sa=5.0, convergence_sa=1.87, qmesh_nr=5, qmesh_nt=2,
                 broadening=0.30, edge_n=None, edge_l=None):
        """
        Initialize a case.innes file to be written
        :param config_dict: set of optional tags to be included (defined in VALID_TELNES3_TAGS)
        :param comment: human-readable comment
        :param absorbing_atom: index to absorbing atom in this calculation's case.struct
        :param edge: absorption edge to calculate. Converted into edge_n, edge_l if they are not
        specified
        :param e_onset: absorption onset energy, in eV
        :param e_beam: exciting electron beam energy, in keV
        :param e_start: starting energy of energy grid
        :param e_final: final energy of energy grid
        :param e_step: step size of energy grid
        :param collection_sa: collection semiangle, in mrad
        :param convergence_sa: convergence semiangle, in mrad
        :param qmesh_nr: NR parameter of Q-mesh
        :param qmesh_nt: NT parameter of Q-mesh
        :param broadening: spectrometer broadening, in eV
        :param edge_n: edge n-quantum number, overrides edge
        :param edge_l: edge l-quantum number, overrides edge
        """
        self.absorbing_atom = absorbing_atom
        self.edge = edge
        self.e_onset = e_onset
        self.e_beam = e_beam
        self.comment = comment
        self.config_dict = deepcopy(config_dict)
        self.e_start = e_start
        self.e_final = e_final
        self.e_step = e_step
        self.collection_sa = collection_sa
        self.convergence_sa = convergence_sa
        self.qmesh_nr = qmesh_nr
        self.qmesh_nt = qmesh_nt
        self.broadening = broadening

        EDGE_DEFAULTS = {'K': [1, 0],
                         'L1': [2, 0],
                         'L23': [2, 1],
                         'M1': [3, 0],
                         'M23': [3, 1],
                         'M45': [3, 2]}

        if (edge_n or edge_l):
            self.edge_n = edge_n
            self.edge_l = edge_l
        else:
            self.edge_n, self.edge_l = EDGE_DEFAULTS[edge]

    @property
    def tags(self):
        """
        TELNES3 optional parameters
        :return: the tags in config_dict
        """
        return TelnesTags(self.config_dict)

    def __str__(self):
        """
        String representation of case.innes
        :return: string
        """
        output = ['{:<60}'.format(self.comment),
                  str(self.absorbing_atom),
                  ' '.join([str(self.edge_n), str(self.edge_l)]),
                  '{:.2f}'.format(self.e_onset),
                  '{:.0f}'.format(self.e_beam),
                  ' '.join(["{:.4f} {:.4f} {:.4f}".format(self.e_start, self.e_final, self.e_step)]),
                  ' '.join(['{:.2f}'.format(self.collection_sa), '{:.2f}'.format(self.convergence_sa)]),
                  ' '.join([str(self.qmesh_nr), str(self.qmesh_nt)]),
                  '{:.2f}'.format(self.broadening)]
        output.extend(["%s\n%s" % (k, TelnesTags._stringify_val(v))
                       for k, v in six.iteritems(self.config_dict)])

        output.append("END")
        return "\n".join(output)

    def write_file(self, filename):
        """
        Write the case.innes file
        :param filename: filename and path to write to
        :return:
        """
        with zopen(filename, "wt") as f:
            f.write(self.__str__() + "\n")

    @staticmethod
    def from_string(input):
        """
        Reads a case.innes from input
        :param input: input string
        :return:
        """

        lines = list(clean_lines(input.split('\n')))
        try:
            comment = lines[0]
            absorbing_atom = int(lines[1])
            edge_n, edge_l = [int(x) for x in lines[2].split()]
            e_onset = float(lines[3])
            e_beam = float(lines[4])
            e_start, e_final, e_step = [float(x) for x in lines[5].split()]
            collection_sa, convergence_sa = [float(x) for x in lines[6].split()]
            qmesh_nr, qmesh_nt = [int(x) for x in lines[7].split()]
            broadening = float(lines[8])
        except:
            warnings.warn("Formatting error in mandatory lines of " + filename)
            return None

        config_dict = TelnesTags.fromstring('\n'.join(lines[9:]))

        return Innes(config_dict, comment, absorbing_atom, '', e_onset, e_beam, e_start, e_final,
                     e_step, collection_sa, convergence_sa, qmesh_nr, qmesh_nt, broadening, edge_n,
                     edge_l)

    @staticmethod
    def from_file(filename):
        """
        Reads case.innes, pulling mandatory parameters and optional tags
        :param filename: filename to read
        :return: Innes object
        """

        with zopen(filename, "rt") as f:
            return Innes.from_string(f.read())

class Control(dict, MSONable):
    """
    Class to create my wien2k.in shell script. Results in a bash-script to be executed
    by WIEN2KCMD
    """

    def __init__(self, params=None):
        """
        Creates a control object
        :param params: Dict of input environment variables to set.
        """
        super(Control, self).__init__()
        if params:
            self.update(params)

    def __setitem__(self, key, val):
        """
        Add parameter-val pair to Control. Cleans the parameter and val by stripping
        leading and trailing white space.
        :param key: key to set
        :param value: value to assign to key
        """
        super(Control, self).__setitem__(
            key.strip(), Control.proc_val(key.strip(), val.strip())
            if isinstance(val, six.string_types) else val)

    def as_dict(self):
        d = dict(self)
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        return d

    @classmethod
    def from_dict(cls, d):
        return Control({k: v for k, v in d.items() if k not in ("@module",
                                                                "@class")})

    def get_string(self, sort_keys=False, pretty=False):
        """
        Returns a string representation of the Control file. Provides options for
        pretty printing
        :param sort_keys: (bool) set to True to sort the parameters alphabetically. Defaults
            to False
        :param pretty: Set to True for pretty aligned output. Defaults to False.
        :return: string representation
        """
        keys = self.keys()
        if sort_keys:
            keys = sorted(keys)

        lines = []
        for k in keys:
            lines.append(f'export {k}="{self[k]}"')

        lines.append(f'eval {self["COMMAND"]}') #execute the command last

        if pretty:
            return str(tabulate([[l.split('=')[0], "=", l.split('=')[-1]] for l in lines[0:-1]], tablefmt = "plain"))
        else:
            return "\n".join(lines)

    def __str__(self):
        return self.get_string(sort_keys=True, pretty=False)

    def write_file(self, filename):
        """
        Write Control to a file.
        :param filename: (str) filename to write to.
        """
        with zopen(filename, "w", newline='\n') as f:  #Write in binary to force UNIX CR
            f.write(self.__str__())

    @staticmethod
    def from_file(filename):
        """
        Reads a Control object from a file
        :param filename: (str) filename of the file to be read
        :return: Incar object
        """
        with zopen(filename, "rt") as f:
            return Control.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads a Control object from a string. Does not import the final execution command.
        :param string: (str) to read
        :return: Control object
        """
        lines = list(clean_lines(string.splitlines()))
        params = {}
        for line in lines:
            m = re.match(r'export (\w+)\s*=\s*(.*)', line)
            if m:
                key = m.group(1).strip()
                val = m.group(2).strip()
                val = Control.proc_val(key,val)
                params[key] = val
        return Control(params)

    @staticmethod
    def proc_val(key, val):
        """
        Static helper method to convert Control parameters to proper types, e.g.,
        integers, floats, lists, etc.
        :param key: Control parameter key
        :param val: Actual value of Control parameter
        :return:
        """
        float_keys = ("WIENEC", "WIENCC", )
        int_keys = ("WIENFC")

        def smart_int_or_float(numstr):
            if numstr.find(".") != -1 or numstr.lower().find("e") != -1:
                return float(numstr)
            else:
                return int(numstr)

        try:
            if key in float_keys:
                return float(re.search(r"^-?\d*\.?\d*[e|E]?-?\d*", val).group(0))

            if key in int_keys:
                return int(re.match(r"^-?[0-9]+", val).group(0))

        except ValueError:
            pass

        # Not in standard keys. We will try a hierarchy of conversions.
        try:
            val = int(val)
            return val
        except ValueError:
            pass

        try:
            val = float(val)
            return val
        except ValueError:
            pass

        if "true" in val.lower():
            return True

        if "false" in val.lower():
            return False

        return val.strip()

    def diff(self, other):
        """
        Diff function for Control. Compares two Controls and indicates which
        parameters are the same and which are not. Useful for checking whether
        two runs were done using the same parameters.
        :param other: Other Control object to compare to.
        :return: Dict of the following format:
            {"Same" : parameters_that_are_the_same,
            "Different": parameters_that_are_different}
            Note that the parameters are return as full dictionaries of values.
            E.g. {"ISIF":3}
        """
        similar_param = {}
        different_param = {}
        for k1, v1 in self.items():
            if k1 not in other:
                different_param[k1] = {"CONTROL1": v1, "CONTROL2": None}
            elif v1 != other[k1]:
                different_param[k1] = {"CONTROL1": v1, "CONTROL2": other[k1]}
            else:
                similar_param[k1] = v1
        for k2, v2 in other.items():
            if k2 not in similar_param and k2 not in different_param:
                if k2 not in self:
                    different_param[k2] = {"CONTROL1": None, "CONTROL2": v2}
        return {"Same": similar_param, "Different": different_param}

    def __add__(self, other):
        """
        Add all the values of another Control object to this object.
        Facilitates the use of "standard" Control.
        """
        params = {k: v for k, v in self.items()}
        for k, v in other.items():
            if k in self and v != self[k]:
                raise ValueError("Controls have conflicting values!")
            else:
                params[k] = v
        return Control(params)

class In0(MSONable):
    """
    Case.in0 class
    """
    pass

class In1(MSONable):
    """
    Case.in1 class
    """

class In2(MSONable):
    """
    Case.in2 class
    """