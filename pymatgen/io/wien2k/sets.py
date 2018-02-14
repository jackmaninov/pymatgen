# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals, print_function

import abc
import re
import os
import glob
import shutil
import warnings
from itertools import chain
from copy import deepcopy

import six
import numpy as np

from monty.serialization import loadfn
from monty.io import zopen

from pymatgen.core.structure import Structure
from pymatgen.io.wien2k.inputs import Struct, Kpoints, In0, In1, In2, Innes, Control
from pymatgen.io.wien2k.outputs import Scffile, Dayfile
from monty.json import MSONable
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.sites import PeriodicSite

"""
This module defines the Wien2kInputSet abstract base class concrete
implementations for the parameters to use in calculations.  The basic 
concept behind an input set is to specify a scheme to generate a consistent 
set of WIEN2k inputs from a structure without further user intervention. 
This ensures comparability across runs.

Read the following carefully before implementing new input sets:

1. 99% of what needs to be done can be done by specifying user_settings
   to override some of the defaults of various input sets. Unless there is an
   extremely good reason to add a new set, DO NOT add one. E.g., if you want
   to turn the hubbard U off, just set "LDAU": False as a user_setting.
2. All derivative input sets should inherit from one of the usual MPRelaxSet or
   MITRelaxSet, and proper superclass delegation should be used where possible.
   In particular, you are not supposed to implement your own as_dict or
   from_dict for derivative sets unless you know what you are doing.
   Improper overriding the as_dict and from_dict protocols is the major
   cause of implementation headaches. If you need an example, look at how the
   MPStaticSet or MPNonSCFSets are constructed.

The above are recommendations. The following are UNBREAKABLE rules:
1. All input sets must take in a structure or list of structures as the first
   argument.
2. user_settings and user_kpoints_settings are absolute. Any new sets you
   implement must obey this. If a user wants to override your settings,
   you assume he knows what he is doing. Do not magically override user
   supplied settings. You can issue a warning if you think the user is wrong.
3. All input sets must save all supplied args and kwargs as instance variables.
   E.g., self.my_arg = my_arg and self.kwargs = kwargs in the __init__. This
   ensures the as_dict and from_dict work correctly.
"""

__author__ = "Eamon McDermott"
__version__ = "1.0"
__maintainer__ = "Eamon McDermott"
__email__ = "eamon.mcdermott@gmail.com"
__date__ = "Feb 5 2018"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

class Wien2kInputSet(six.with_metaclass(abc.ABCMeta, MSONable)):
    """
    Base class representing a set of WIEN2k input parameters (a case) with a structure
    supplied as init parameters. Typically, you should not inherit from this class. Start
    from DictSet or MPRelaxSet or MITRelaxSet.
    """

    @property
    @abc.abstractmethod
    def control(self):
        """Command-and-control object"""
        pass

    @property
    @abc.abstractmethod
    def struct(self):
        """Struct object"""
        pass

    @property
    @abc.abstractmethod
    def kpoints(self):
        """Kpoints object"""
        pass

    @property
    @abc.abstractmethod
    def in0(self):
        """case.in0 file"""
        pass

    @property
    @abc.abstractmethod
    def in1(self):
        """case.in1 file"""
        pass

    @property
    @abc.abstractmethod
    def in2(self):
        """case.in2 file"""
        pass

    @property
    @abc.abstractmethod
    def innes(self):
        """Innes object"""
        pass

    #TODO more input files

    @property
    def all_input(self):
        """
        Returns all input files as a dict of {filename: wien2k object}
        """

        #case = self.case

        return {'control': self.control,
                'struct': self.struct,
                'innes': self.innes,
                'klist_pymatgen': self.kpoints,
                'in0': self.in0,
                'in1': self.in1,
                'in2': self.in2}


    def write_input(self, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        """
        Writes a set of WIEN2k input to a directory
        :param output_dir: Directory to output the WIEN2k input files
        :param make_dir_if_not_present: Set to True if you want the
                directory (and the whole path) to be created if it is not
                present.
        :param include_cif: Whether to write a CIF file in the output
                directory for easier opening by VESTA.
        :return:
        """

        # we MUST write to a directory with the case name
        #if os.path.basename(os.path.normpath(output_dir)) != self.case:
        #    output_dir = os.path.join(output_dir, self.case)
        #    warnings.warn(f"Output directory MUST match case. Writing to {output_dir}")

        case = os.path.basename(os.path.normpath(output_dir)) # files must match dir name

        if make_dir_if_not_present and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for k, v in self.all_input.items():
            v.write_file(os.path.join(output_dir, f"{case}.{k}"))
        if include_cif:
            s = self.all_input['struct'].structure
            fname = os.path.join(output_dir, "%s.cif" % re.sub(r'\s', "",
                                                               s.formula))
            s.to(filename=fname)

    def as_dict(self, verbosity=2):
        """
        output set as a dict
        :param verbosity: verbosity level
        :return: dict representation
        """
        d=MSONable.as_dict(self)
        if verbosity == 1:
            d.pop("structure", None)
        return d

def _load_yaml_config(fname):
    #TODO are there actually any parameters I would want to load by default into wien2k.in?

    config = loadfn(os.path.join(MODULE_DIR, "%s.yaml" % fname))
    config["Control"].update(loadfn(os.path.join(MODULE_DIR, "WIENControlBase.yaml")))
    return config

class DictSet(Wien2kInputSet):
    """
    Concrete implementation of Wien2kInputSet that is initialized from a dict
    settings. This allows arbitrary settings to be input. In general,
    this is rarely used directly unless there is a source of settings in yaml
    format (e.g., from a REST interface). It is typically used by other
    Wien2kInputSets for initialization.

    Args:
        structure (Structure): The structure to create inputs for.
        case (string): optional case name (defaults to pymatgen)
        config_dict (dict): The config dictionary to use.
        files_to_transfer (dict): A dict of {filename: filepath}. This allows the
            transfer of files from a previous calculation.
        user_control_settings (string): TODO decide how to pass in control object
        user_kpoints_settings (dict or Kpoints): Allow user to overrid kpoints setting
            by supplying a dict e.g. ("reciprocal_density": 1000).
            User can also supply Kpoints object. Default is None
            TODO update this with WIEN2k kpoints definition
        sort_structure (bool): Whether to sort the structure (using the
            default sort order of electronegativity) before generating input
            files. Defaults to True, the behavior you would want most of the
            time. This ensures that similar atomic species are grouped
            together.
        functional (str): Functional to use in case.in0. Valid vaules: "XC_LDA",
            "XC_PBE", "XC_PBESOL", "XC_TPSS", "XC_SCAN", "XC_WC", "XC_MBJ", etc
            TODO make a translation table object
        force_gamma (bool): Force gamma centered kpoint generation (e.g. shift k-mesh
            is false). Defaults to False for Monkhorst-Pack
        reduce_structure (None/str): efore generating the input files,
            generate the reduced structure. Default (None), does not
            alter the structure. Valid values: None, "niggli", "LLL".
        TODO vdw: Add default van-der-Waals correction.
    """

    def __init__(self, structure, config_dict, case="pymatgen", files_to_transfer = None,
                 user_control_settings = None, user_kpoints_settings = None,
                 sort_structure = True, functional="XC_PBE", force_gamma=False,
                 reduce_structure=None, vdw=None):
        if reduce_structure:
            structure = structure.get_reduced_structure(reduce_structure)
        if sort_structure:
            structure = structure.get_sorted_structure()
        self.structure = structure
        self._config_dict = deepcopy(config_dict)
        self.case = case
        self.files_to_transfer = files_to_transfer or {}
        self.sort_structure = sort_structure
        self.functional = functional #TODO roll this into an override?
        self.force_gamma = force_gamma
        self.reduce_structure = reduce_structure
        self.user_control_settings = user_control_settings or {}
        self.user_kpoints_settings = user_kpoints_settings
        self.vdw = vdw.lower() if vdw is not None else None #TODO probably will detect off functional

        """if self.vdw:
                set up DFTD3 params
        """

    @property
    def control(self):
        settings = dict(self._config_dict["Control"])
        settings.update(self.user_control_settings)
        structure = self.structure
        control = Control() #TODO implement in inputs.py
        comp = structure.composition
        """
        elements = sorted([el for el in comp.elements if comp[el] > 0],
                          key=lambda e: e.X)
        most_electroneg = elements[-1].symbol
        struct = Struct(structure)
        hubbard_u = settings.get("LDAU", False)
        """

        for k, v in settings.items():
            """
            TODO add structure-dependant control settings here
            e.g. VASP DictSet sets LDAU parameters, so I could turn on -orb and
            set up Indm/Inorb objects. MPRelax set also sets the E convergence
            per atom"""
            control[k] = v


        return control

    @property
    def struct(self):
        return Struct(self.structure)

    @property
    def nelect(self):
        """
        Gets the default number of electrons for a given structure.
        :return: number of electrons
        """
        return int(round(
            sum([self.structure.composition.element_composition[ps.element]
                 * ps.ZVAL
                 for ps in self.struct])))

    @property
    def in0(self):
        return In0()

    @property
    def in1(self):
        return In1()

    @property
    def in2(self):
        return In2()

    @property
    def innes(self):
        settings = self.user_control_settings or self._config_dict["Innes"]
        #TODO break out configuration from the user_control_settings.
        #TODO Need to think of a good place to put the "Atom of interest" into the control set...
        return Innes(settings)

    @property
    def kpoints(self):
        """
        Writes out a case.klist file using the fully automated grid method. TODO Should run x kgen
        through some other function!
        :return:
        """
        settings = self.user_kpoints_settings or self._config_dict["KList"]

        if isinstance(settings, Kpoints):
            return settings

        if settings.get('grid_density'):
            kpoints = Kpoints.automatic_density(
                self.structure, int(settings['grid_density']),
                self.force_gamma)
            self.user_control_settings["WIEN_KPOINTS"] = kpoints.num_kpts

        elif settings.get('reciprocal_density'):
            return Kpoints.automatic_density_by_vol(
                self.structure, int(settings['reciprocal_density']),
                self.force_gamma)

        # Raise error. Unsure of which kpoint generation to use
        else:
            raise ValueError(
                "Invalid KPoint Generation algo : Supported Keys are "
                "grid_density: for Kpoints.automatic_density generation, "
                "reciprocal_density: for KPoints.automatic_density_by_vol "
                "generation, and length  : for Kpoints.automatic generation")

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return self.__class__.__name__

    def write_input(self, output_dir,
                    make_dir_if_not_present=True, include_cif=False):
        super(DictSet, self).write_input(output_dir=output_dir,
                                         make_dir_if_not_present=make_dir_if_not_present,
                                         include_cif=include_cif)
        for k,v in self.files_to_transfer.items():
            with zopen(v, "rb") as fin, \
                    zopen(os.path.join(output_dir,k), "wb") as fout:
                shutil.copyfileobj(fin, fout)
        #TODO need to see how well this copies big files...



class LETIRelaxSet(DictSet):
    """
    TODO implement a set of default WIEN2k settings for my project. Should consist of
    things like:
    RKMAX
    Functional
    executing commands
    """
    CONFIG = _load_yaml_config("LETIRelaxSet")

    def __init__(self, structure, **kwargs):
        super(LETIRelaxSet, self).__init__(
            structure, LETIRelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs




class LETIStaticSet(DictSet):
    """
    TODO same as Relax set. Important to implement from_prev_calc to do my imports
    """

    def __init__(self, structure, prev_control=None, prev_kpoints=None,
                 reciprocal_density=100, **kwargs):
        super(LETIStaticSet, self).__init__(structure, **kwargs)
        if isinstance(prev_control, six.string_types):
            prev_control = Control.from_file(prev_control)
        if isinstance(prev_kpoints, six.string_types):
            prev_kpoints = Kpoints.from_file(prev_kpoints)

        self.prev_control = prev_control
        self.prev_kpoints = prev_kpoints
        self.reciprocal_density = reciprocal_density
        self.structure = structure
        self.kwargs = kwargs

    @property
    def control(self):
        parent_control = super(LETIStaticSet, self).control
        control = Control(self.prev_control) if self.prev_control is not None else \
            Control(parent_control)

        # Control.update() #set some reset defaults

        control["WIENEC"] = min(control.get("WIENEC", 1), parent_control["WIENEC"])
        return control

    @property
    def kpoints(self):
        return kpoints

    @classmethod
    def from_prev_calc(cls, prev_calc_dir, standardize=False, sym_prec=0.1,
                       international_monoclinic=True, reciprocal_density=100,
                       small_gap_multiply=None, **kwargs):
        """
        Generate a set of WIEN2k input files for static calculations from a directory of
        a previous WIEN2k run
        :param prev_calc_dir:
        :param standardize:
        :param sym_prec:
        :param international_monoclinic:
        :param reciprocal_density:
        :param small_gap_multiply:
        :param kwargs:
        :return:
        """

        case = os.path.basename(os.path.normpath(prev_calc_dir))

        prev_control = Control.from_file(os.path.join(prev_calc_dir, 'wien2k.in'))
        prev_struct = Struct.from_file(os.path.join(prev_calc_dir, f"{case}.struct"))
        prev_kpoints = Kpoints.from_file(os.path.join(prev_calc_dir, f"{case}.klist"))
        prev_scffile = Scffile(os.path.join(prev_calc_dir, f"{case}.scf"))

        # multiply the reciprocal density if needed:
        if small_gap_multiply:
            gap = prev_scffile.gap
            if gap <= small_gap_multiply[0]:
                reciprocal_density *= small_gap_multiply[1]

        return LETIStaticSet(structure=prev_struct.structure, prev_control=prev_control,
                             prev_kpoints=prev_kpoints, reciprocal_density=reciprocal_density,
                             **kwargs)

