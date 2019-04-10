# coding: utf-8

from __future__ import unicode_literals

from fireworks.core.firework import FiretaskBase, FWAction
import pymatgen.io.wien2k.outputs

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
        base_name (str): Name of the file to create, including the path,
            minus any format-specific extension.
        format (str): Optional. one of "zip", "tar", "bztar" or "gztar".
        Defaults to gztar.
    """

    _fw_name = 'ArchiveWIEN2kOutputTask'
    required_params = ["output_file"]
    optional_params = ["format"]

    def run_task(self, fw_spec):
        print(self.get("format"))
        toSave = self.get("format").from_file(self.get("output_file"))
        return FWAction(stored_data=toSave.as_dict())
