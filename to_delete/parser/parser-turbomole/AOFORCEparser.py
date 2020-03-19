# Copyright 2016-2018 Arvid Conrad Ihrig, Aliaksei Mazheika
#                     Fritz-Haber-Institut der Max-Planck-Gesellschaft
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""This module constructs the parser for the AOFORCE module from TurboMole"""

import logging
import re
import os
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.simple_parser import AncillaryParser
from turbomoleparser.TurbomoleCommon import RE_FLOAT, RE_FIXED_FLOAT
import turbomoleparser.TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class AOFORCEparser(object):

    mode_units = {
        # FIXME: add other units for mode energies supported by turbomole
        "cm**(-1)": "inversecm"
    }

    def __init__(self, context, key="aoforce"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__hessian = None
        self.__modes = None
        self.__energies = None
        self.__energy_unit = None
        self.__activity_infrared = None
        self.__activity_raman = None
        self.__intensities = None
        self.__current_columns = None
        self.__current_row = None
        self.__prefix_width = None
        self.__vibration_files = dict()

    def purge_data(self):
        self.__hessian = None
        self.__modes = None
        self.__energies = None
        self.__activity_infrared = None
        self.__activity_raman = None
        self.__intensities = None
        self.__current_columns = None
        self.__current_row = None
        self.__prefix_width = None
        self.__vibration_files = dict()

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{5,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(r"\s*a o f o r c e - program\s*$",
                    name="header",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        sub_matchers = [
            header,
            self.__context["geo"].build_qm_geometry_matcher(),
            self.__context["geo"].build_orbital_basis_matcher(),
            self.__context["orbitals"].build_ir_rep_matcher(),
            self.__context["method"].build_dft_functional_matcher(),
            self.__context["method"].build_dftd3_vdw_matcher(),
            self.__build_hessian_matcher(),
            self.__build_normal_mode_file_matcher(),
            self.__build_vibration_spectrum_file_matcher(),
            self.__build_normal_modes_matcher(),
            self.__build_total_energy_matcher(),
        ]

        return self.__context.build_module_matcher("force", sub_matchers, "AOFORCE",
                                                   self.__context["method"].add_default_functional)

    def __build_hessian_file_matcher(self):
        def store_file_name(backend, groups):
            self.__vibration_files["hessian"] = groups[0]
        return SM(r"\s*\*{3}\s*projected\s+hessian\s+written\s+onto\s+\$hessian\s+\(projected\),"
                  r"\s+file=<([^>]+)>\s*\*{3}\s*$",
                  name="hessian file",
                  startReAction=store_file_name
                  )

    def __build_normal_mode_file_matcher(self):
        def store_file_name(backend, groups):
            self.__vibration_files["modes"] = groups[0]
        return SM(r"\s*\*{3}\s*normal\s+modes\s+written\s+onto\s+\$vibrational\s+normal\s+modes,"
                  r"\s+file=<([^>]+)>\s*\*{3}\s*$",
                  name="normal modes file",
                  startReAction=store_file_name
                  )

    def __build_vibration_spectrum_file_matcher(self):
        def store_file_name(backend, groups):
            self.__vibration_files["spectrum"] = groups[0]

        sub = SM("\s*file\s*=\s*<(.+)>\s*\*{3}\s*$",
                 name="spectrum filename",
                 required=True,
                 startReAction=store_file_name)
        return SM(r"\s*\*{3}\s*vibrational\s+spectroscopic\s+data\s+written\s+"
                  r"onto\s*\$vibrational\s+spectrum\s*$",
                  name="spectrum filename",
                  subMatchers=[sub]
                  )

    def __build_hessian_matcher(self):
        re_header = re.compile(r"\s*([0-9]+)\s*[A-z]+")
        re_element = re.compile(r"\s*("+RE_FIXED_FLOAT+")")
        indices = {"x": 0, "y": 1, "z": 2}

        def prepare_data(backend, groups):
            n_atoms = self.__context["geo"].num_atoms()
            self.__hessian = np.ndarray(shape=(n_atoms, n_atoms, 3, 3), dtype=float)
            self.__hessian[:, :, :, :] = 0.0
            self.__current_columns = list()
            self.__current_row = -1

        def process_header(backend, groups):
            match = re_header.match(groups[0])
            del self.__current_columns[:]
            while match:
                self.__current_columns.append(int(match.group(1)) - 1)
                match = re_header.match(groups[0], pos=match.end(0))

        def process_row(backend, groups):
            if groups[0]:
                self.__current_row = int(groups[0]) - 1
            axis1 = indices[groups[1]]
            match = re_element.match(groups[2])
            index = 0
            while match:
                axis2 = index % 3
                atom2 = self.__current_columns[index // 3]
                value = float(match.group(0))
                self.__hessian[self.__current_row, atom2, axis1, axis2] = value
                self.__hessian[atom2, self.__current_row, axis2, axis1] = value
                match = re_element.match(groups[2], pos=match.end(0))
                index += 1

        def write_data(backend, gIndex, section):
            if "hessian" in self.__vibration_files:
                hessian_parser = self.__build_hessian_file_parser()
                dir_name = os.path.dirname(os.path.abspath(self.__context.fName))
                file_name = os.path.normpath(os.path.join(dir_name,
                                                          self.__vibration_files["hessian"]))
                try:
                    with open(file_name) as file_in:
                        hessian_parser.parseFile(file_in)
                except IOError:
                    logger.warning("Could not find auxiliary file '%s' in directory '%s'." % (
                        self.__vibration_files["modes"], dir_name))
            backend.addArrayValues("hessian_matrix", self.__hessian,
                                   self.__context.index_configuration(), unit="hartree * bohr**-2")

        block_row = SM(r"\s*(?:([0-9]+)+\s*[A-z]+\s*)?d([xyz])\s+((?:"+RE_FIXED_FLOAT+"\s*)+)\s*$",
                       name="hessian row",
                       repeats=True,
                       startReAction=process_row
                       )

        block_header = SM(r"\s*ATOM((?:\s+[0-9]+\s*[A-z]+)+)\s*$",
                          name="column header",
                          startReAction=process_header,
                          repeats=True,
                          subMatchers=[
                              SM(r"\s*dx\s+dy\s+dz\s+dx\s+dy\s+dz\s*$", name="column header"),
                              block_row
                          ]
                          )

        return SM(r"\s*CARTESIAN\s+FORCE\s+CONSTANT\s+MATRIX\s+\(hartree/bohr\*\*2\)\s*$",
                  name="Cartesian Hessian",
                  startReAction=prepare_data,
                  subMatchers=[
                      block_header,
                      self.__build_hessian_file_matcher()
                  ],
                  onClose={None: write_data}
                  )

    def __build_total_energy_matcher(self):

        def set_total_energy(backend, groups):
            backend.addRealValue("energy_total", float(groups[0]), unit="hartree")

        scf_etot = SM(r"\s*\*\s*SCF-energy\s+:\s+("+RE_FLOAT+")\s*\*\s*$",
                      name="SCF total energy")
        scf_zp = SM(r"\s*\*\s*SCF\s+\+\s+E\(vib0\)\s+:\s+"
                    r"(?P<energy_current__hartree>"+RE_FLOAT+")\s*\*\s*$",
                    name="SCF total energy",
                    startReAction=set_total_energy)

        return SM(r"\s*\*\s*zero\s+point\s+VIBRATIONAL\s+energy\s+:\s+"
                  r"(?P<energy_zero_point__hartree>"+RE_FLOAT+")\s*Hartree\s*\*\s*$",
                  name="zero point energy",
                  # TODO: define allowed values for ZP-corrections and add output
                  subMatchers=[
                      scf_etot,
                      scf_zp
                  ])

    def __prepare_vibrations__data(self, backend, groups):
        if self.__modes is None:
            n_atoms = self.__context["geo"].num_atoms()
            self.__modes = np.ndarray(shape=(3 * n_atoms, n_atoms, 3), dtype=float)
            self.__modes[:, :, :] = 0.0
            self.__energies = np.ndarray(shape=(3 * n_atoms), dtype=float)
            self.__energies[:] = 0.0
            try:
                self.__energy_unit = self.mode_units[groups[0]]
            except KeyError:
                logger.error("unknown vibrational mode energy unit found: " + groups[0])
                self.__energy_unit = self.mode_units["cm**(-1)"]
            self.__intensities = np.ndarray(shape=(3 * n_atoms), dtype=float)
            self.__intensities[:] = 0.0
            self.__activity_infrared = np.ndarray(shape=(3 * n_atoms), dtype=bool)
            self.__activity_infrared[:] = False
            self.__activity_raman = np.ndarray(shape=(3 * n_atoms), dtype=bool)
            self.__activity_raman[:] = False
            backend.addValue("x_turbomole_vibrations_num_modes", 3 * n_atoms)
            self.__current_columns = list()
            self.__current_row = -1

    def __build_normal_modes_matcher(self):
        re_header = re.compile(r"\s*([0-9]+)")
        re_element = re.compile(r"\s*("+RE_FLOAT+")")
        re_activity = re.compile(r"\s*(-|YES)")
        indices = {"x": 0, "y": 1, "z": 2}

        def process_header(backend, groups):
            match = re_header.match(groups[0])
            del self.__current_columns[:]
            while match:
                self.__current_columns.append(int(match.group(1)) - 1)
                match = re_header.match(groups[0], pos=match.end(0))

        def process_row(backend, groups):
            if groups[0]:
                self.__current_row = int(groups[0]) - 1
            axis1 = indices[groups[1]]
            match = re_element.match(groups[2])
            index = 0
            while match:
                value = float(match.group(0))
                self.__modes[self.__current_columns[index], self.__current_row, axis1] = value
                match = re_element.match(groups[2], pos=match.end(0))
                index += 1

        def extract_frequencies(backend, groups):
            match = re_element.match(groups[0])
            index = 0
            while match:
                value = float(match.group(0))
                self.__energies[self.__current_columns[index]] = value
                match = re_element.match(groups[0], pos=match.end(0))
                index += 1

        def extract_infrared_activity(backend, groups):
            match = re_activity.match(groups[0])
            index = 0
            while match:
                value = match.group(0) == "YES"
                self.__activity_infrared[self.__current_columns[index]] = value
                match = re_activity.match(groups[0], pos=match.end(0))
                index += 1

        def extract_raman_activity(backend, groups):
            match = re_activity.match(groups[0])
            index = 0
            while match:
                value = match.group(0) == "YES"
                self.__activity_raman[self.__current_columns[index]] = value
                match = re_activity.match(groups[0], pos=match.end(0))
                index += 1

        def extract_intensities(backend, groups):
            match = re_element.match(groups[0])
            index = 0
            while match:
                value = float(match.group(0))
                self.__intensities[self.__current_columns[index]] = value
                match = re_element.match(groups[0], pos=match.end(0))
                index += 1

        def write_data(backend, gIndex, section):
            # check if the external files are present and read them if possible (higher precision)
            if "modes" in self.__vibration_files:
                modes_parser = self.__build_normal_modes_file_parser()
                dir_name = os.path.dirname(os.path.abspath(self.__context.fName))
                file_name = os.path.normpath(os.path.join(dir_name,
                                                          self.__vibration_files["modes"]))
                try:
                    with open(file_name) as file_in:
                        modes_parser.parseFile(file_in)
                except IOError:
                    logger.warning("Could not find auxiliary file '%s' in directory '%s'." % (
                        self.__vibration_files["modes"], dir_name))
            if "spectrum" in self.__vibration_files:
                spectrum_parser = self.__build_vibrational_spectrum_file_parser()
                dir_name = os.path.dirname(os.path.abspath(self.__context.fName))
                file_name = os.path.normpath(os.path.join(dir_name,
                                                          self.__vibration_files["spectrum"]))
                try:
                    with open(file_name) as file_in:
                        spectrum_parser.parseFile(file_in)
                except IOError:
                    logger.warning("Could not find auxiliary file '%s' in directory '%s'." % (
                        self.__vibration_files["spectrum"], dir_name))

            backend.addArrayValues("x_turbomole_vibrations_normal_modes", self.__modes,
                                   self.__context.index_configuration(), unit="bohr")
            backend.addArrayValues("x_turbomole_vibrations_mode_energies", self.__energies,
                                   self.__context.index_configuration(), unit=self.__energy_unit)
            backend.addArrayValues("x_turbomole_vibrations_intensities", self.__intensities,
                                   self.__context.index_configuration(), unit="kcal * mole ** -1")
            backend.addArrayValues("x_turbomole_vibrations_infrared_activity",
                                   self.__activity_infrared, self.__context.index_configuration())
            backend.addArrayValues("x_turbomole_vibrations_raman_activity", self.__activity_raman,
                                   self.__context.index_configuration())

        frequencies = SM(r"\s*frequency((?:\s+"+RE_FLOAT+")+)\s*$",
                         name="frequencies",
                         required=True,
                         startReAction=extract_frequencies
                         )
        intensities = SM(r"\s*intensity\s*\(km/mol\)((?:\s+"+RE_FLOAT+")+)\s*$",
                         name="intensities",
                         required=True,
                         startReAction=extract_intensities
                         )
        activity_infrared = SM(r"\s*IR((?:\s+(?:-|YES))+)\s*$",
                               name="IR activity",
                               required=True,
                               startReAction=extract_infrared_activity
                               )
        activity_raman = SM(r"\s*RAMAN((?:\s+(?:-|YES))+)\s*$",
                            name="RAMAN activity",
                            required=True,
                            startReAction=extract_raman_activity
                            )

        block_row = SM(r"\s*(?:([0-9]+)\s+[A-z]+\s+)?([xyz])((?:\s+"+RE_FLOAT+")+)\s*$",
                       name="mode displacements",
                       repeats=True,
                       startReAction=process_row
                       )

        block_header = SM(r"\s*mode((?:\s+[0-9]+)+)\s*$",
                          name="column header",
                          startReAction=process_header,
                          repeats=True,
                          subMatchers=[
                              frequencies,
                              activity_infrared,
                              intensities,
                              activity_raman,
                              block_row
                          ]
                          )

        return SM(r"\s*NORMAL\s+MODES\s+and\s+VIBRATIONAL\s+FREQUENCIES\s+\((.+)\)\s*$",
                  name="vibrational spectrum",
                  startReAction=self.__prepare_vibrations__data,
                  subMatchers=[
                      block_header
                  ],
                  onClose={None: write_data}
                  )

    def __build_hessian_file_parser(self):
        re_element = re.compile(r"\s*("+RE_FLOAT+")")

        def prepare_data(backend, gIndex, section):
            self.__prefix_width = len(str(self.__context["geo"].num_atoms()))+1
            self.__current_columns = 0

        def process_line(backend, groups):
            row_index = int(groups[0][0:self.__prefix_width])
            if self.__current_row != row_index:
                self.__current_columns = 0
                self.__current_row = row_index
            atom_row = (row_index - 1) // 3
            axis_row = (row_index - 1) % 3
            match = re_element.match(groups[1])
            while match:
                atom_col = self.__current_columns // 3
                axis_col = self.__current_columns % 3
                value = float(match.group(0))
                self.__hessian[atom_row, atom_col, axis_row, axis_col] = value
                match = re_element.match(groups[1], pos=match.end(0))
                self.__current_columns += 1

        data = SM(r"(\s*[0-9]+\s*[0-9]+)((?:\s+"+RE_FLOAT+")+)\s*$",
                  name="hessian matrix",
                  startReAction=process_line,
                  repeats=True)
        normal_modes = SM(r"",
                          name="hessian matrix file root",
                          onOpen={None: prepare_data},
                          subMatchers=[
                              SM(r"\$hessian \(projected\)\s*$", name="start", required=True),
                              data,
                              SM(r"\$end\s*$", name="end", required=True)
                          ]
                          )

        return AncillaryParser(fileDescription=normal_modes,
                               parser=self.__context.parser,
                               cachingLevelForMetaName=dict(),
                               superContext=self.__context
                               )

    def __build_normal_modes_file_parser(self):
        re_element = re.compile(r"\s*("+RE_FLOAT+")")

        def prepare_data(backend, gIndex, section):
            self.__prefix_width = len(str(self.__context["geo"].num_atoms()))+1
            self.__current_columns = 0

        def process_line(backend, groups):
            row_index = int(groups[0][0:self.__prefix_width])
            if self.__current_row != row_index:
                self.__current_columns = 0
                self.__current_row = row_index
            atom = (row_index - 1) // 3
            axis = (row_index - 1) % 3
            match = re_element.match(groups[1])
            while match:
                value = float(match.group(0))
                self.__modes[self.__current_columns, atom, axis] = value
                match = re_element.match(groups[1], pos=match.end(0))
                self.__current_columns += 1

        data = SM(r"(\s*[0-9]+\s*[0-9]+)((?:\s+"+RE_FLOAT+")+)\s*$",
                  name="displacements",
                  startReAction=process_line,
                  repeats=True)
        normal_modes = SM(r"",
                          name="normal modes file root",
                          onOpen={None: prepare_data},
                          subMatchers=[
                              SM(r"\$vibrational normal modes\s*$", name="start", required=True),
                              data,
                              SM(r"\$end\s*$", name="end", required=True)
                          ]
                          )

        return AncillaryParser(fileDescription=normal_modes,
                               parser=self.__context.parser,
                               cachingLevelForMetaName=dict(),
                               superContext=self.__context
                               )

    def __build_vibrational_spectrum_file_parser(self):

        def process_header(backend, groups):
            try:
                self.__energy_unit = self.mode_units[groups[0]]
            except KeyError:
                logger.error("unknown vibrational mode energy unit found: " + groups[0])
                self.__energy_unit = self.mode_units["cm**(-1)"]

        def process_line(backend, groups):
            mode_index = int(groups[0]) - 1
            self.__energies[mode_index] = float(groups[1])
            self.__intensities[mode_index] = float(groups[2])
            self.__activity_infrared[mode_index] = groups[3] == "YES"
            self.__activity_raman[mode_index] = groups[4] == "YES"

        header = SM(r"#\s+mode\s+symmetry\s+wave\s+number\s+IR\s+intensity\s+selection\s+rules\s*$",
                    name="header",
                    required=True
                    )
        header_units = SM(r"#\s+(\S+)\s+km/mol\s+IR\s+RAMAN",
                          name="units",
                          required=True,
                          startReAction=process_header
                          )
        data = SM(r"\s*([0-9]+)\s+[a-z'\"]*\s+(" + RE_FLOAT + ")\s+("
                  + RE_FLOAT + ")\s+(-|YES)\s+(-|YES)\s*$",
                  name="vibrational mode",
                  startReAction=process_line,
                  repeats=True)
        normal_modes = SM(r"",
                          name="vibrational spectrum file root",
                          subMatchers=[
                              SM(r"\$vibrational spectrum\s*$", name="start", required=True),
                              header,
                              header_units,
                              data,
                              SM(r"\$end\s*$", name="end", required=True)
                          ]
                          )

        return AncillaryParser(fileDescription=normal_modes,
                               parser=self.__context.parser,
                               cachingLevelForMetaName=dict(),
                               superContext=self.__context
                               )
