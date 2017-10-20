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
on case.struct, case.inc, case.inm, case.inso, case.int, case.inorb
"""

__author__ = "Eamon McDermott"
__version__ = "0.1"
__email__ = "eamon.mcdermott@cea.fr"
__status__ = "Development"
__date__ = "Oct 19, 2017"

logger = logging.getLogger(__name__)

class Site(MSONable):
    """
    WIEN2k site within a structure

    """

    def __init__(self, sitenum, coords, multi, isplit, NPT=781, R0, RMT, Z):
        if structure.is_ordered:



class Struct(MSONable):
    """
    case.struct reader/writer
    """

    def __init__ (self, structure, comment=None, relativistic=True):
        if structure.is_ordered:
            site_properties = {}
            self.comment = structure.formula if comment is None else comment
            self.relativistic = relativistic
        else:
            raise ValueError("Structure with partial occupancies cannot be "
                             "converted into a WIEN2k struct!")

    @property
    def natoms(self):
        """
        Number of ATOM sites #TODO understand this
        :return:

        syms = [site.specie.symbol for site in self.structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]
        """
        return self.natoms
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
        reader2=ff.FortranRecordReader('(A4,23X,I3)')
        try:
            structtype=LATTYPE.get(reader2.read(lines[1])[0].rstrip())
        except IndexError:
            raise ValueError("Unknown lattice type")
        try:
            natoms = reader2.read(lines[1])[1]
        except:
            raise ValueError("Unknown NAT")

        #line 3:
        reader3=ff.FortranRecordReader('(13X,A4)')
        try:
            relativistic = RELA.get(reader3.read(lines[2])[0])
        except
            raise ValueError("Unclear relativistic mode in STRUCT")

        #line 4:
        reader4=ff.FortranRecordReader('(6F10.6)')
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


        #line 5 (repeats for natoms)

        reader5 = ff.FortranRecordReader('(4X,I4,4X,F10.8,3X,F10.8,3X,F10.8)')
        reader6 = ff.FortranRecordReader('(15X,I2,17X,I2)')
        reader7 = ff.FortranRecordReader('(A10,5X,I5,5X,F10.8,5X,F10.5,5X,F10.5)')
        reader8 = ff.FortranRecordReader('(20X,3F10.7)')
        reader11 = ff.FortranRecordReader('(I4)')
        reader12 = ff.FortranRecordReader('(3I2,F10.7)')
        reader15 = ff.FortranRecordReader('(I8)')

        sites = []
        currentline = 4

        for i in range(natoms):
            try:
                [atomindex, x, y, z] = reader5.read(lines[currentline])
            except:
                raise ValueError('Error reading coordinates of Atom ', i)
            currentline += 1
            try:
                [multiplicity, isplit] = reader6.read(lines[currentline])
            except:
                raise ValueError('Bad multiplicity/ISPLIT on Atom ', i)
            currentline += 1
            altsites=[]
            for j in range(multiplicity):           #Read alternate site coordinates
                try:
                    altsites.append(reader5.read(lines[currentline]))
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
                except:
                    raise ValueError('Error reading ROTLOC on site ', i)
            currentline += 3

            #TODO append a Site object to sites() here
            if debug:
                print(atomindex, x, y, z, ROTLOC)

        nsym = reader11.read(lines[currentline])
        currenline += 1

        symops=[]

        for i in range(nsym):
            symop=[]
            for symopline in lines[currentline:currentline + 3]:
                symop.append(reader12.read(symopline))
            currentline += 3
            symops.append(symop)
            symopindex = reader15(lines[currentline])
            currentline += 1
            assert (symopindex = i+1,'Misordered symmetry operations'

    return natoms

class Kpoints_supported_modes(Enum):
    Automatic = 0
    Gamma = 1
    Monkhorst = 2

    def __str__(self):
        return self.name

    @staticmethod
    def from_string(s):
        c = s.lower()[0]
        for m in Kpoints_supported_modes:
            if m.name.lower()[0] == c:
                return m
        raise ValueError("Can't interprete Kpoint mode %s" % s)

class Kpoints(MSONable):
    """
    case.klist reader/writer
    """

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


