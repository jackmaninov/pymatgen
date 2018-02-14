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
            warnings.warn('Float overflow (*******) encountered in WIEN2krun')
            return np.nan
        raise e


# def keep_last(seq):
#     seen = {}
#     result = []
#     for item in seq:
#         if item in seen: continue
#         seen[marker] = 1
#         result.append(item)
#     return result


class Dayfile(MSONable):
    """
    Parse the WIEN2k case.dayfile for important SCF results (convergence, early termination, etc)
    """

    def __init__(self, filename):
        self.filename = filename
        self.is_stopped = False

        ext_run_stats = []
        run_stats = {}
        energy_convergence = []
        charge_convergence = []
        header = []
        cycle_start = []
        running_prog = None
        current_cycle = None
        last_utime = None
        last_stime = None
        processors = 1

        prog_patt = re.compile(r"^\>")
        stop_patt = re.compile(r"crashed|\.stop|stop error")
        utime_patt = re.compile(r"(^|\W)([0-9]*\.[0-9]{3})u")
        stime_patt = re.compile(r"(^|\W)([0-9]*\.[0-9]{3})s")

        all_lines = []
        for line in zopen(self.filename):
            clean = line.strip()
            all_lines.append(clean)

            if stop_patt.search(clean):
                self.is_stopped = True
            else:
                if clean.find("cycle") != -1:
                    tok = clean.strip().split()
                    current_cycle = tok[1]
                    cycle_start.append([current_cycle, ' '.join(tok[2:8])])
                    continue
                if current_cycle is None:
                    header.append(clean)
                if prog_patt.search(clean):
                    tok = clean.strip().split()
                    if running_prog is not None:
                        ext_run_stats.append([current_cycle, running_prog, last_utime, last_stime])
                    running_prog = tok[1]
                if utime_patt.search(clean):
                    # only keep the last (summary) time for parallel programs
                    last_utime = utime_patt.search(clean).group(2)
                    last_stime = stime_patt.search(clean).group(2)
                if clean.find(":ENERGY") != -1:
                    energy_convergence.append([current_cycle, clean.strip().split()[-1]])
                if clean.find(":CHARGE") != -1:
                    charge_convergence.append([current_cycle, clean.strip().split()[-1]])
                if clean.find("processors") != -1:
                    tok = clean.strip().split()
                    processors = tok[-2]

        if running_prog is not None and last_utime is not None:
            ext_run_stats.append([current_cycle, running_prog, last_utime, last_stime])

        self.ext_run_stats = ext_run_stats
        self.user_time = round(sum([float(prog[2]) for prog in ext_run_stats]), 3)  # approximate
        self.sys_time = round(sum([float(prog[3]) for prog in ext_run_stats]), 3)

        run_stats['User time (sec)'] = self.user_time
        run_stats['System time (sec)'] = self.sys_time
        run_stats['cores'] = processors

        self.run_stats = run_stats
        self.energy_convergence = energy_convergence
        self.charge_convergence = charge_convergence
        self.all_lines = all_lines
        self.header = header

class Scffile(MSONable):
    """
    Parser for a WIEN2k case.scf file, where most output is collected
    """

    nsites = None

    def __init__(self, filename):
        self.filename = filename

        header = []
        last_iteration = None
        spin = None
        potentials = []
        ifft = []
        ifft_enhancement = None
        density = None
        lattice_constant = None
        total_energy = None
        total_energy_warning = None
        efermi = None
        nelect = None
        total_mag = None
        magnetization = []
        charge = []
        chargeline = {}
        gap = None

        for line in reverse_readfile(self.filename):
            clean = line.strip()
            if clean.find(":ITE") != -1:
                tok = re.split('[.:]', clean)
                last_iteration = int(tok[2])
            if ":NATO" in clean:
                tok = clean.split()
                self.nsites, self.natoms = int(tok[2]), int(tok[5])
            if clean.find(":ENE") != -1:
                tok = clean.split('=')
                total_energy = float(tok[-1])
                total_energy_warning = (clean.find("WARNING") != -1)
            if clean.find(":FER") != -1:
                tok = clean.split('=')
                efermi = float(tok[-1])
            if clean.find(":NOE") != -1:
                tok = clean.split('=')
                nelect = float(tok[-1])
            if clean.find("SPINPOLARIZED") != -1:
                spin = (clean.find("NON-SPINPOLARIZED") == -1)
            if clean.find(":MMTOT") != -1:
                spin = True
                tok = clean.split('=')
                total_mag = float(tok[-1])
            if clean.find(":MMI") != -1:
                tok = clean.split()
                magnetization.append(float(tok[-1]))
            if clean.find(":GAP") != -1:
                tok = clean.split()
                gap = float(tok[6])
            if clean.find(":CTO") != -1:
                tok = clean.split()
                if spin:
                    chargeline['up'], chargeline['dn'], chargeline['tot'] = float(tok[-7]), \
                        float(tok[-4]), float(tok[-1])
                else:
                    chargeline['tot'] = float(tok[-1])
                if clean.find("INTERSTITIAL") != -1:
                    charge.insert(0, dict(chargeline))  #put INTERSTITAL at the beginning so when we reverse it is at end
                else:
                    charge.append(dict(chargeline))

            if all([nelect is not None, efermi is not None,
                    spin is not None,
                    last_iteration is not None]):
                break

        self.magnetization = tuple(magnetization[::-1])
        self.charge = tuple(charge[::-1])
        self.spin = spin
        self.total_energy = total_energy
        self.total_energy_warning = total_energy_warning
        self.efermi = efermi
        self.nelect = nelect
        self.total_mag = total_mag
        self.last_iteration = last_iteration
        self.data = {}
        self.gap = gap

        """
        TODO other items read in pymatgen.io.vasp.outputs.py.Outcar can be extracted
        from BerryPI calculations (Piezoelectric epsilon, Born effective charge)
        
        WIENncm exists to deal with non-collinear spin structures, but is not included
        with the main distribution.
        """

    def read_core_state_eigen(self):
        """
        Read out the core states
        TODO: better error handling
        :return:
            A list of dict over the atoms [{"AO":[core state eig]}]
            The core state eigenenergy list for each AO is over each SCF iteration
        """

        cl = [defaultdict(list) for i in range(self.nsites)]
        with zopen(self.filename, "rt") as fscf:
            line = fscf.readline()
            while line != "":
                line = fscf.readline()
                if ":ITE" in line:
                    tok = re.split('[.:]', line)
                    current_iteration = int(tok[2])
                if ":1S" in line:
                    while (line != "\n" and line.find("Ry") != -1):
                        tok = re.split('[ :]', line)
                        try:
                            cl[int(line[4:7])-1][str.lower(line[9:12]).strip()].append(float(tok[-2]))
                        except:
                            continue
                        line = fscf.readline()

        return cl


    @property
    def is_spin(self):
        """
        :return: returns true if run is spin-polarized
        """
        return self.spin

class OutputNMR(MSONable):
    """
    TODO: Parse the NMR chemical shifts of an x_nmr_lapw run

    :param self:
    :return:
    """

    def __init__(self, filename):
        """
        Read and parse case.
        :param filename: case.outputnmr to parse
        """
        warnings.warn('TODO: WIEN2k NMR parsing not implemented')


class Eels(MSONable):
    """
    Parse case.broadspec for EELS spectrum.
    """

    def __init__(self, data, filename):
        self.data = np.array(data)
        self.filename = filename

    @property
    def energies(self):
        """
        Returns the spectrum energies in eV.
        :return: spectrum energies
        """

        return self.data[:, 0]

    @property
    def total_spectrum(self):
        """
        Returns the total EELS spectrum
        :return: total spectrum
        """

        return self.data[:, 1]

    @property
    def first_edge(self):
        """
        Returns the first edge of the EELS spectrum
        :return: first edge
        """

        if len(self.data[0])>2:
            print("Moo")
            return self.data[:, 2]
        else:
            warnings.warn("No first edge in " + self.filename)

    @property
    def second_edge(self):
        """
        Returns the second edge of the EELS spectrum
        :return: second edge
        """

        if len(self.data[0])>3:
            return self.data[:, 3]
        else:
            warnings.warn("No second edge in " + self.filename)
            return None


    @staticmethod
    def from_file(spectrum):
        """
        Load and parse a broadened/unbroadened EELS spectrum
        :param spectrum: filename and path for case.broadspec
        :return: Eels object
        """
        data = np.loadtxt(spectrum)
        return Eels(data, spectrum)

    def as_dict(self):
        """
        Returns dict representation of Eels object
        :return: dict representation
        """
        d = MSONable.as_dict(self)
        d["data"] = self.data.tolist()
        return d