# coding: utf-8

from __future__ import unicode_literals

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks import Firework, Workflow, ScriptTask, FileTransferTask
import pymatgen.io.wien2k.outputs
from pymatgen.io.wien2k.inputs import Innes

__author__ = 'Eamon McDermott'
__copyright__ = 'Copyright 2019'
__version__ = '0.1'
__maintainer__ = 'Shyue Ping Ong'
__email__ = 'ongsp@ucsd.edu'
__date__ = 'Apr 10, 2019'


class ArchiveWIEN2kOutputTask(FiretaskBase):
    """
    Load a WIEN2k output file (using classes from outputs.py) and store
    it in stored_data

    Args:
        output_file (str): Name of the file to create, including the path,
            minus any format-specific extension.
        format (str): Optional. one of "zip", "tar", "bztar" or "gztar".
        Defaults to gztar.
    """

    _fw_name = 'ArchiveWIEN2kOutputTask'
    required_params = ["output_file"]
    optional_params = ["format"]

    def run_task(self, fw_spec):
        toSave = eval(self.get("format")).from_file(self.get("output_file"))
        return FWAction(stored_data=toSave.as_dict())

class DeployInnesInputTask(FiretaskBase):
    """
    Write a WIEN2k Innes input file from spec

    Args:
        case_name (str): Name of file to create
        Innes (dict): Innes to create
    """

    _fw_name = 'DeployInnesInputTask'
    required_params = ["output_file"]
    optional_params = ["Innes"]   #TODO do I pass a dict or just what I want to change?

    def run_task(self, fw_spec):
        a = self.get("Innes")
        a.write_file(self.get("output_file"))

class TelnesRunTask(FiretaskBase):
    """
    Run TELNES given an Innes spec

    Args:
        case_name (str): Name of case (no extension)
        origin (str): originating directory to copy
        Innes (dict): Innes spec
    """

    _fw_name = 'TelnesRunTask'
    required_params = ["case_name", "origin"]
    optional_params = ["Innes"]

    @staticmethod
    def TelnesList(case='case', path=''):
        TelnesExtensions = ['kgen', 'innes', 'inc', 'inb', 'corewavef', 'final', 'vsp', 'vtotal',
                            'struct', 'qtl', 'rotij', 'ortho', 'ctr', 'sdlm', 'matrix', 'cdos', 'dos',
                            'xdos', 'sp2', 'angular', 'eelstable']
        return ([path + (case + '.' + ext) for ext in TelnesExtensions])

    def run_task(self, fw_spec):
        case = self.get("case_name")
        mkdirTask = ScriptTask.from_str('mkdir ' + case)
        deployTask = FileTransferTask(
            {'files': self.TelnesList(case, self.get("origin")), 'dest': case,
             'mode': 'copy'})  # , 'ignore_errors':'True'
        deployTask2 = DeployInnesInputTask({'output_file': case + '/' + case + ".innes", 'Innes': self.get("Innes")})
        runTask = ScriptTask.from_str('cd '+ case +' && x telnes3 && x broadening >> STDOUT')
        archiveTask = ArchiveWIEN2kOutputTask( {'output_file': case + '/' + case + '.broadspec',
                                                'format': 'pymatgen.io.wien2k.outputs.Eels'})
        cleanTask = ScriptTask.from_str('yes|rm -rf '+case)
        return FWAction(detours=Firework([mkdirTask, deployTask, deployTask2, runTask, archiveTask]))  #, cleanTask]))

class EelsAngleSweepTask(FiretaskBase):
    """
    Generate all Eels output over a range of angles. Run a TelnesRunTasks for each

    Args:
        case_name (str): Name of the WIEN2k case (no extension)
        origin (str): originating directory to copy
        Innes (dict): Base Innes spec (angles will be overwritten)
        step (int): step size
    """

    _fw_name = 'EelsAngleSweepTask'
    required_params = ["case_name", "origin", "Innes"]
    optional_params = ["step", "phase", "component", "method", "optimized"]

    @staticmethod
    def frange(start, stop, step):
        i = start
        while i <= stop:
            yield i
            i += step

    def run_task(self, fw_spec):
        step = float(self.get("step"))
        if not step:
            step = 15.0

        tasklist=[]


        baseInnes = self.get("Innes")

        for alpha in self.frange(0, 90, step):
            for beta in self.frange(0, 90, step):
                for gamma in self.frange(0, 90, step):
                    caseInnes=Innes.from_dict(baseInnes.as_dict())
                    caseInnes.config_dict["ORIENTATION SENSITIVE"] = [alpha, beta, gamma]
                    tasklist.append(Firework(TelnesRunTask({"case_name": self.get("case_name"),
                                                         "origin": self.get("origin"),
                                                         "Innes": caseInnes.as_dict()})))
        return FWAction(additions=tasklist)