# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import glob
import itertools
import logging
import math
import os
import re
import warnings
import xml.etree.cElementTree as ET
from collections import defaultdict
from io import StringIO
import fortranformat as ff

import numpy as np
from monty.io import zopen, reverse_readfile
from monty.json import MSONable
from monty.json import jsanitize
from monty.re import regrep
from six import string_types
from six.moves import map, zip

from pymatgen.analysis.nmr import NMRChemicalShiftNotation
from pymatgen.core.composition import Composition
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Structure
from pymatgen.core.units import unitized
from pymatgen.electronic_structure.bandstructure import BandStructure, \
    BandStructureSymmLine, get_reconstructed_band_structure
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType, Magmom
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.entries.computed_entries import \
    ComputedEntry, ComputedStructureEntry
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar
from pymatgen.util.io_utils import clean_lines, micro_pyawk

"""
Classes for reading/manipulating/writing WIEN2k ouput files.
"""


logger = logging.getLogger(__name__)


def _wien2krun_float(f):
    """
    Large numbers are often represented as ********* in WIEN output.
    This function parses these values as np.nan
    """
    try:
        return float(f)
    except ValueError as e:
        f = f.strip()
        if f == '*' * len(f):
            warnings.warn('Float overflow (*******) encountered in vasprun')
            return np.nan
        raise e

class dayfile(MSONable)
    """
    Parse the WIEN2k case.dayfile for important SCF results
    """

    def __init__(self,filename):
        self.filename = filename
        self.is_stopped = False

        run_stats = {}
        header = []

    time_patt = re.compile(r"\(()")

    all_lines = []
    for line in reverse_readfile(self.filename):
        clean = line.strip()
        all_lines.append(clean)
        if clean.find("stop due to .stop file") != -1:
            self.is_stopped = True
        else:


class Scffile(MSONable):
    """
    Parser for a WIEN2k case.scf file, where most output is collected
    """

    def __init__(self, filename):
        self.filename = filename
        self.is_stopped = False

        header = []
        iteration = None
        spinpolarized = None
        potentials = []
        ifft = []
        ifft_enhancement = None
        density = None
        lattice_c

        all_lines = []
        for line in reverse_readfile(self.filename):
            clean = line.strip()
            all_lines.append(clean)
            if clean.find