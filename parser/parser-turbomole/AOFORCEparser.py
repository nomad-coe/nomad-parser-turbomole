"""This module constructs the parser for the AOFORCE module from TurboMole"""

import logging
import re
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT, RE_FIXED_FLOAT
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class AOFORCEparser(object):

    def __init__(self, context, key="aoforce"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__hessian = None
        self.__modes = None
        self.__energies = None
        self.__intensities = None
        self.__current_columns = None
        self.__current_row = None

    def purge_data(self):
        pass

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
            self.__context.build_start_time_matcher(),
            header,
            self.__context["geo"].build_qm_geometry_matcher(),
            self.__context["geo"].build_orbital_basis_matcher(),
            self.__context["orbitals"].build_ir_rep_matcher(),
            self.__context["method"].build_dft_functional_matcher(),
            self.__context["method"].build_dftd3_vdw_matcher(),
            self.build_hessian_matcher(),
            self.build_normal_modes_matcher(),
            self.build_total_energy_matcher(),
            self.__context.build_end_time_matcher("grad")
        ]

        return self.__context.build_module_matcher("force", sub_matchers, "AOFORCE",
                                                   self.__context["method"].add_default_functional)

    def build_hessian_matcher(self):
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
                      block_header
                  ],
                  onClose={None: write_data}
                  )

    def build_total_energy_matcher(self):

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

    def build_normal_modes_matcher(self):
        re_header = re.compile(r"\s*([0-9]+)")
        re_element = re.compile(r"\s*("+RE_FLOAT+")")
        indices = {"x": 0, "y": 1, "z": 2}
        units = {"cm**(-1)": "???"}

        def prepare_data(backend, groups):
            n_atoms = self.__context["geo"].num_atoms()
            self.__modes = np.ndarray(shape=(3 * n_atoms, n_atoms, 3), dtype=float)
            self.__modes[:, :, :] = 0.0
            self.__energies = np.ndarray(shape=(3 * n_atoms), dtype=float)
            self.__energies[:] = 0.0
            self.__intensities = np.ndarray(shape=(3 * n_atoms), dtype=float)
            self.__intensities[:] = 0.0
            backend.addValue("x_turbomole_vibrations_num_modes", 3 * n_atoms)
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

        def extract_intensities(backend, groups):
            match = re_element.match(groups[0])
            index = 0
            while match:
                value = float(match.group(0))
                self.__intensities[self.__current_columns[index]] = value
                match = re_element.match(groups[0], pos=match.end(0))
                index += 1

        def write_data(backend, gIndex, section):
            backend.addArrayValues("x_turbomole_vibrations_normal_modes", self.__modes,
                                   self.__context.index_configuration(), unit="bohr")
            backend.addArrayValues("x_turbomole_vibrations_mode_energies", self.__energies,
                                   self.__context.index_configuration(), unit="inversecm")
            backend.addArrayValues("x_turbomole_vibrations_intensities", self.__intensities,
                                   self.__context.index_configuration(), unit="kcal * mole ** -1")

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
                              intensities,
                              block_row
                          ]
                          )

        return SM(r"\s*NORMAL\s+MODES\s+and\s+VIBRATIONAL\s+FREQUENCIES\s+\(cm\*\*\(-1\)\)\s*$",
                  name="vibrational spectrum",
                  startReAction=prepare_data,
                  subMatchers=[
                      block_header
                  ],
                  onClose={None: write_data}
                  )
