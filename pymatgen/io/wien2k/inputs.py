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

class Kpoints(MSONable):
    """
    case.klist reader/writer
    """

