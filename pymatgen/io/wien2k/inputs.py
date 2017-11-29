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
import numpy as np
from numpy.linalg import det
from collections import OrderedDict, namedtuple
from hashlib import md5
import fortranformat as ff

from monty.io import zopen
from monty.os.path import zpath
from monty.json import MontyDecoder

from enum import Enum
from tabulate import tabulate

import scipy.constants as const

from pymatgen import SETTINGS
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element, get_el_sp
from monty.design_patterns import cached_class
from pymatgen.util.string_utils import str_delimited
from pymatgen.util.io_utils import clean_lines
from monty.json import MSONable

"""
Classes for reading/manipulating/writing WIEN2k input files. Initial focus is
on case.struct, case.inc, case.inm, case.inso, case.int, case.inorb, case.klist
"""



logger = logging.getLogger(__name__)

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
        super(Struct, self).__setattr__(name, value) # how is this not recursive?

    @staticmethod
    def from_file(filename, debug=0):
        """
        read from a case.struct file
        :param filename: case.struct name (string)
        :return: Struct
        """

        with zopen(filename, "rt") as f:
            return Struct.from_string(f.read(), debug)

    @staticmethod
    def from_string(data, debug=0):
        """

        :param data: string containing case.struct
        :return: Struct
        """

        LATTYPE = {'P': 1, 'F': 2, 'B': 3, 'CXY': 4, 'CYZ': 5, 'CXZ': 6, 'R': 7, 'H': 8}
        RELA = {'RELA': True, 'NREL': False}
        UNITS  = {'ang': 1, 'bohr': 2}



        chunks = re.split("\n\s*\n", data.rstrip(), flags=re.MULTILINE)
        try:                                            #TODO: why?
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
            structtype = LATTYPE.get(reader2.read(lines[1])[0].rstrip())
        except IndexError:
            raise ValueError("Unknown lattice type")
        try:
            natoms = reader2.read(lines[1])[1]
        except:
            raise ValueError("Unknown NAT")
        if debug: print("LATTYPE=", structtype)
        if debug: print("number atoms=", natoms)

        #line 3:
        reader3 = ff.FortranRecordReader('(13X,A4)')
        try:
            relativistic = RELA.get(reader3.read(lines[2])[0])
        except:
            raise ValueError("Unclear relativistic mode in STRUCT")

        if debug: print("RELA", relativistic)
        if debug: print("in Bohrs?: ", inbohr)

        #line 4:
        reader4 = ff.FortranRecordReader('(6F10.6)')
        try:
            [a, b, c, alpha, beta, gamma] = reader4.read(lines[3])
            assert (a is not None or b is not None or c is not None), 'Missing unit cell parameter'
            if alpha is None: alpha = 90.0
            if beta is None: beta = 90.0
            if gamma is None: gamma = 90.0
        except:
            raise ValueError("Error reading unit cell parameters")

        #convert bohrs to angstroms
        bohr = const.physical_constants['atomic unit of length'][0] * 10**(10)
        a *= bohr
        b *= bohr
        c *= bohr
        #construct lattice
        if structtype is LATTYPE['P']:              #primative (not hexagonal)
            lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        elif structtype is LATTYPE['F']:            #face-centered
            lattice = Lattice([[ a/2,  b/2,    0],
                               [ a/2,   0,   c/2],
                               [   0,  b/2,  c/2]])
        elif structtype is LATTYPE['B']:            #body-centered
            lattice = Lattice([[ a/2, -b/2,  c/2],
                               [ a/2,  b/2, -c/2],
                               [-a/2,  b/2,  c/2]])
        elif structtype is LATTYPE['CXY']:          #C-base-centered (orthorhombic only)
            lattice = Lattice([[ a/2, -b/2,    0],
                               [ a/2,  b/2,    0],
                               [   0,    0,    c]])
        elif structtype is LATTYPE['CYZ']:          #A-base-centered (orthorhombic only)
            lattice = Lattice([[   a,    0,    0],
                               [   0, -b/2,  c/2],
                               [   0,  b/2,  c/2]]) #B-base-centered (orthorh. and monoclinic symmetry)
        elif structtype is LATTYPE['CXZ']:
            gamma_r = math.radians(gamma)
            lattice = Lattice([[ a * math.sin(gamma_r)/2, a * math.cos(gamma_r)/2,  -c],
                               [ 0, b, 0],
                               [ a * math.sin(gamma_r)/2, a * math.cos(gamma_r)/2, c/2]])
        elif structtype is LATTYPE['R']:            #WIEN2k rhombohedral setting
            lattice = Lattice([[ a / (math.sqrt(3)/2), -a/2, c/3],
                               [ a / (math.sqrt(3)/2),  a/2, c/3],
                               [-a / math.sqrt(3),        0, c/3]])
        elif structtype is LATTYPE['H']:            #WIEN2k hexagonal setting
            lattice = Lattice([[math.sqrt(3) * a / 2, -a / 2, 0],
                               [0, a, 0],
                               [0, 0, c]])
        else:
            raise ValueError('Unknown STRUCT type')

        #TODO implement hex/rhomb coordinates

        if debug: print("Lattice: ", lattice)

        #line 5 (repeats for natoms)

        reader5 = ff.FortranRecordReader('(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8)')
        reader6 = ff.FortranRecordReader('(15X,I2,17X,I2)')
        reader7 = ff.FortranRecordReader('(A10,5X,I5,5X,F10.8,5X,F10.5,5X,F10.5)')
        reader8 = ff.FortranRecordReader('(20X,3F10.7)')
        reader11 = ff.FortranRecordReader('(I4)')
        reader12 = ff.FortranRecordReader('(3I2,F10.7)')
        reader15 = ff.FortranRecordReader('(I8)')

        species = []
        coords = []
        sitenames = []
        NPTs = []
        R0s = []
        RMTs = []
        Zs = []
        altcoordss = []
        ROTLOCs = []


        currentline = 4

        for i in range(natoms):
            # site line 1 (ex: ATOM    1: X=...)
            try:
                [atomindex, x, y, z] = reader5.read(lines[currentline])
            except:
                raise ValueError('Error reading coordinates of Atom ', i)
            if debug: print("Atomindex, x, y, z: ", atomindex, x, y, z)
            currentline += 1
            # site line 2 (ex: MULT= 1 ISPLIT = 2)
            try:
                [multiplicity, isplit] = reader6.read(lines[currentline])
            except:
                raise ValueError('Bad multiplicity/ISPLIT on Atom ', i)
            currentline += 1
            if debug: print("multiplicity, isplit: ", multiplicity,isplit)
            altcoords=[]
            for j in range(multiplicity - 1):           #Read alternate site coordinates
                if debug: print("loop j: ",j)
                try:
                    if debug: print("on line ", currentline, ": ", lines[currentline])
                    if debug: print (reader5.read(lines[currentline]))
                    altcoords.append(reader5.read(lines[currentline])[1:])
                except:
                    raise ValueError('Improper multiplicity on Atom ', i)
                currentline += 1
            try:
                [atomname, NPT, R0, RMT, Z] = reader7.read(lines[currentline])
            except:
                raise ValueError('Bad site parameters on Atom ', i)


            currentline += 1
            ROTLOC=[]
            for rotlocline in lines[currentline:currentline + 3]:
                try:
                    ROTLOC.append(reader8.read(rotlocline))
                    if debug: print ("ROTLOC:", ROTLOC[-1])
                except:
                    raise ValueError('Error reading ROTLOC on site ', i)
            currentline += 3

            for site in [[float(x), float(y), float(z)]] + altcoords:
                species.append(Element.from_Z(int(Z)))  #NOTE: WIEN2k supports non-integer Z!
                coords.append(site)
                sitenames.append(atomname.strip())
                NPTs.append(NPT)
                R0s.append(R0)
                RMTs.append(RMT)
                Zs.append(Z)
                altcoordss.append(altcoords)
                ROTLOCs.append(ROTLOC)

            if debug:
                print(atomindex, x, y, z, ROTLOC)

        [nsym] = reader11.read(lines[currentline])
        if debug: print("NSYM: ", nsym)
        currentline += 1

        symops=[]

        for i in range(nsym):
            symop=[]
            for symopline in lines[currentline:currentline + 3]:
                symop.append(reader12.read(symopline))
            currentline += 3
            symops.append(symop)
            symopindex = reader15.read(lines[currentline])
            currentline += 1
            # assert symopindex = i+1,'Misordered symmetry operations'
            if debug: print ("Symop ", i, ": ", symops[-1])

        #build structure

        site_properties={'sitename': sitenames, 'NPT': NPTs, 'R0': R0s,
                         'RMT': RMTs, 'Z': Zs,
                         'ROTLOC': ROTLOCs}

        if debug: print(site_properties)
        if debug: print (species)
        if debug: print(coords)
        if debug: print(site_properties)

        struct = Structure(lattice, species, coords, to_unit_cell=False,
                           validate_proximity=False,coords_are_cartesian=False,
                           site_properties=site_properties)

        return Struct(struct, comment=comment, relativistic=relativistic, symops=symops,
                      lattice=lattice)

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
        lines.append(writer4.write(list(latt.abc)+list(latt.angles)))

        # line 5 - 7
        writer5 = ff.FortranRecordWriter('(A4,I4,A4,F10.8,A3,F10.8,A3,F10.8)')
        writer6 = ff.FortranRecordWriter('(A15,I2,A17,I2)')
        writer7 = ff.FortranRecordWriter('(A10,A5,I5,A5,F10.8,A5,F10.5,A5,F10.5)')
        writer8 = ff.FortranRecordWriter('(A20,3F10.7)')

        for (i,site) in enumerate(self.structure):
            coords = site.frac_coords
            lines.append(writer5.write(['ATOM', i, ': X=', coords[0], ' Y=',
                                       coords[1], ' Z=', coords[2]]))
            lines.append(writer6.write(['MULT=', 1, 'ISPLIT=', 0]))
            lines.append(writer7.write([str(site.specie).ljust(10),
                                        ' NPT=', self.structure.site_properties['NPT'][i],
                                        '  R0=', self.structure.site_properties['R0'][i],
                                        ' RMT=', self.structure.site_properties['RMT'][i],
                                        '   Z:', site.Z]))
            lines.append(writer8.write(['LOCAL ROT MATRIX:'.ljust(20), 1, 0, 0]))
            lines.append(writer8.write(['', 0, 1, 0]))
            lines.append(writer8.write(['', 0, 0, 1]))
            #TODO extract LOCROT, but otherwise let x symmetry generate it

        # lines 11 - 15

        writer11 = ff.FortranRecordWriter('(I4,A35)')
        writer12 = ff.FortranRecordWriter('(3I2,F10.7)')
        writer15 = ff.FortranRecordWriter('(I8)')

        lines.append(writer11.write([len(self.symops), 'NUMBER OF SYMMETRY OPERATIONS'.rjust(35)]))
        for (i,symop) in enumerate(self.symops):
            for line in symop:
                lines.append(writer12.write(line))
            lines.append(writer15.write([i]))

        return "\n".join(lines)


class Klist_supported_modes(Enum):
    Automatic = 0
    Gamma = 1
    Monkhorst = 2

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s):
        c = s.lower()[0]
        for m in Klist_supported_modes:
            if m.name.lower()[0] == c:
                return m
        raise ValueError("Can't interprete Kpoint mode %s" % s)

class Klist(MSONable):
    """
    case.klist reader/writer
    """
    supported_modes = Klist_supported_modes

    def __init__(self, comment="Default gamma", num_kpts=0,
                 style=supported_modes.Gamma,
                 kpts=((1, 1, 1),), kpts_shift=(0, 0, 0),
                 kpts_weights=None, coord_type=None, labels=None,
                 tet_number=0, tet_weight=0, tet_connctions=None):
        """
        WIEN2k Kpoint class constructor.

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
            tet_connections: For explicit kpoints, specifies the connections¬
                of the tetrahedrons for the tetrahedron method.¬
                format is a list of tuples, [ (sym_weight, [tet_vertices]),¬
                ...]¬
        The default behavior of the construtor is for a Gamma centered,
        1x1x1 k-point grid with no shift.
        """

        self.comment = comment
        self.num_kpts = num_kpts
        self.kpts = kpts
        self.style = style
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
            style = Klist.supported_modes.from_string(style)

        if style in (Klist.supported_modes.Automatic,
                     Klist.supported_modes.Gamma,
                     Klist.supported_modes.Monkhorst) and len(self.kpts) > 1:
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
        return Klist("Fully automatic kpoint scheme", 0,
                     style=Klist.supported_modes.Automatic,
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
        return Klist("Automatic kpoint scheme", 0,
                     Klist.supported_modes.Gamma, kpts=[kpts],
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
        return Klist("Automatic kpoint scheme", 0,
                     Klist.supported_modes.Monkhorst, kpts=[kpts],
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
        return Klist("Automatic kpoint scheme", 0,
                     Klist.supported_modes.Monkhorst, kpts=[kpts],
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
            Klist
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
            style = Klist.supported_modes.Gamma
        else:
            style = Klist.supported_modes.Monkhorst

        return Klist(comment, 0, style, [num_div], [0, 0, 0])

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
            Klist
        """
        vol = structure.lattice.reciprocal_lattice.volume
        kppa = kppvol * vol * structure.num_sites
        return Klist.automatic_density(structure, kppa,
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
            return Klist.from_string(f.read())

    @staticmethod
    def from_string(string):
        """
        Reads a Kpoints object from a case.klist string.

        Args:
            string (str): case.klist string.

        Returns:
            Kpoints object
        """

        reader=ff.FortranRecordReader('(I10,4I10,f5.1)')
        reader1=ff.FortranRecordReader('(I10,4I10,3f5.1,2x,i8,\' k, div: (\',3i3,\')\')')
        kpts = []
        kpts_weights = []


        lines = [line.strip() for line in string.splitlines()]

        for line_number, line in enumerate(lines):
            if line == "END": break
            if line_number == 1:
                try:
                    [KPT, KX, KY, KZ, K_WEIGHT, IDIV, junk1, junk2, num_kpts0, n_divx, n_divy, ndivz ] = reader1.read (line)
                    if KX == KY == KZ == 0: style = Klist.supported_modes.Gamma
                except:
                    raise ValueError('Error reading first k-point', line)
            else:
                try:
                    [KPT, KX, KY, KZ, K_WEIGHT, IDIV] = reader.read (line)
                except:
                    raise ValueError('Error reading k-point:', line)

            kpts.append([float(KX/IDIV), float(KY/IDIV), float(KZ/IDIV)])
            kpts_weights.append(float(K_WEIGHT))


        return Klist.gamma_automatic(kpts) if style == Klist.supported_modes.Gamma \
            else Klist.monkhorst_automatic(kpts)

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
    def from_directory(input_dir, casename = None, optional_files = None):
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
                             (case + ".kgen", Kgen), (case + ".klist", Klist)]:
            fullzpath = zpath(os.path.join(input_dir, fname))
            sub_d[fname.lower()] = ftype.from_file(fullzpath)
        sub_d["optional_files"] = {}
        if optional_files is not None:
            for fname, ftype in optional_files.items():
                sub_d["optional_files"][fname] = \
                    ftype.from_file(os.path.join(input_dir, fname))
        return WIEN2kInput(**sub_d)


