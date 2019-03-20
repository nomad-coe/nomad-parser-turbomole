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

"""This module constructs the parser for the CCSDF12 module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from turbomoleparser.TurbomoleCommon import RE_FLOAT
import turbomoleparser.TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class CCSDF12parser(object):

    def __init__(self, context, key="ccsdf12"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__previous_energy = None
        self.__auxiliary_indices = {"method": -1, "config": -1}
        self.__correlation_indices = dict()

    def purge_data(self):
        self.__previous_energy = None
        self.__auxiliary_indices = {"method": -1, "config": -1}
        self.__correlation_indices = dict()

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{5,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(r"\s*\*\s*C\s*C\s*S\s*D\s*F\s*1\s*2\s+P\s*R\s*O\s*G\s*R\s*A\s*M\s*\*\s*$",
                    name="header",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        sub_matchers = [
            header,
            self.__context["geo"].build_qm_geometry_matcher(),
            self.__context["geo"].build_orbital_basis_matcher(),
            self.__context["method"].build_wave_function_model_matcher(),
            self.__context["geo"].build_auxiliary_basis_matcher(),
            self.__build_mp2_starting_point_matcher(),
            self.__build_cc_scf_matcher(),
            self.__build_ccsdt_correction_matcher(),
        ]

        return self.__context.build_module_matcher("ccsdf12", sub_matchers, "CCSDF12")

    def __create_auxiliary_sections(self, backend, method):
        if self.__context["method"].get_main_method() != method:
            # intermediate calculation, create extra method/config sections for results
            self.__auxiliary_indices["config"] = backend.openSection(
                "section_single_configuration_calculation")
            self.__auxiliary_indices["method"] = backend.openSection("section_method")
            backend.addValue("single_configuration_to_calculation_method_ref",
                             self.__auxiliary_indices["method"], self.__auxiliary_indices["config"])
            backend.addValue("single_configuration_calculation_to_system_ref",
                             self.__context.index_system(), self.__auxiliary_indices["config"])
        else:
            self.__auxiliary_indices["config"] = self.__context.index_configuration()
            self.__auxiliary_indices["method"] = self.__context.index_method()
        self.__context.set_index_overwrites(index_config=self.__auxiliary_indices["config"],
                                            index_method=self.__auxiliary_indices["method"])

    def __close_auxiliary_sections(self, backend, gIndex, section):
        if self.__auxiliary_indices["config"] != self.__context.index_configuration():
            backend.closeSection("section_single_configuration_calculation",
                                 self.__auxiliary_indices["config"])
        if self.__auxiliary_indices["method"] != self.__context.index_method():
            backend.closeSection("section_method", self.__auxiliary_indices["method"])
        self.__auxiliary_indices["config"] = -1
        self.__auxiliary_indices["method"] = -1
        self.__context.set_index_overwrites()

    def __build_mp2_starting_point_matcher(self):

        def setup(backend, gIndex, section):
            self.__create_auxiliary_sections(backend, "MP2")
            self.__correlation_indices["MP2"] = self.__auxiliary_indices["config"]

        return SM(r"\s*Calculate\s+integrals\s+\(ia\|jb\)\s+for\s+MP2\s+start\s+guess\s*$",
                  name="MP2 starting point",
                  subMatchers=[
                      self.__context["method"].build_correlation_energy_matcher()
                  ],
                  onOpen={None: setup},
                  onClose={None: self.__close_auxiliary_sections}
                  )

    def __build_cc_scf_matcher(self):
        references = {"section_single_configuration_calculation": None}

        def extract_iteration_data(backend, groups):
            index_scf = backend.openSection("section_scf_iteration")
            self.__backend.setSectionInfo("section_scf_iteration", index_scf, references)
            if self.__previous_energy:
                backend.addRealValue("energy_change_scf_iteration",
                                     float(groups[0]) - self.__previous_energy, index_scf,
                                     unit="hartree")
            backend.addRealValue("energy_total_scf_iteration", float(groups[0]), index_scf,
                                 unit="hartree")
            self.__previous_energy = float(groups[0])
            backend.closeSection("section_scf_iteration", index_scf)

        def convergence(backend, groups):
            backend.addValue("number_of_scf_iterations", int(groups[0]),
                             self.__auxiliary_indices["config"])

        scf_iteration = SM(r"\s*[0-9]+\s+("+RE_FLOAT+")" + 5 * (r"\s+"+RE_FLOAT) + "\s*$",
                           name="iteration",
                           repeats=True,
                           startReAction=extract_iteration_data
                           )
        scf_converged = SM(r"\s*CC\s+equations\s+converged\s+in\s+([0-9]+)\s+iterations\.\s*$",
                           name="CC SCF converged",
                           startReAction=convergence)

        scf_cycle = SM(r"\s*Iter\.\s+CCSD\s+energy\s+Norm\(Omega\)\s+"
                       r"Norm\(t1\)\s+Norm\(t2\)\s+cpu\s+wall",
                       name="CCSD SCF iterations",
                       subMatchers=[
                           scf_iteration,
                           scf_converged
                       ],
                       )

        def setup(backend, gIndex, section):
            self.__create_auxiliary_sections(backend, "CCSD")
            references["section_single_configuration_calculation"] = \
                self.__auxiliary_indices["config"]
            self.__correlation_indices["CCSD"] = self.__auxiliary_indices["config"]
            if "MP2" in self.__correlation_indices:
                index_link = backend.openSection("section_calculation_to_calculation_refs")
                self.__backend.setSectionInfo("section_calculation_to_calculation_refs",
                                              index_link, references)
                backend.addValue("calculation_to_calculation_ref",
                                 self.__correlation_indices["MP2"], index_link)
                backend.addValue("calculation_to_calculation_kind", "starting_point", index_link)
                backend.closeSection("section_calculation_to_calculation_refs", index_link)

        return SM(r"\s*\*\s*OPTIMIZATION\s+OF\s+THE\s+GROUND\s+"
                  r"STATE\s+CLUSTER\s+AMPLITUDES\s*\*\s*$",
                  name="CC ground state",
                  subMatchers=[
                      scf_cycle,
                      self.__context["method"].build_correlation_energy_matcher()
                  ],
                  onOpen={None: setup},
                  onClose={None: self.__close_auxiliary_sections}
                  )

    def __build_ccsdt_correction_matcher(self):
        references = {"section_single_configuration_calculation": None}

        def setup(backend, gIndex, section):
            self.__create_auxiliary_sections(backend, "CCSD(T)")
            references["section_single_configuration_calculation"] = \
                self.__auxiliary_indices["config"]
            self.__correlation_indices["CCSD(T)"] = self.__auxiliary_indices["config"]
            if "CCSD" in self.__correlation_indices:
                index_link = backend.openSection("section_calculation_to_calculation_refs")
                self.__backend.setSectionInfo("section_calculation_to_calculation_refs",
                                              index_link, references)
                backend.addValue("calculation_to_calculation_ref",
                                 self.__correlation_indices["CCSD"], index_link)
                backend.addValue("calculation_to_calculation_kind", "starting_point", index_link)
                backend.closeSection("section_calculation_to_calculation_refs", index_link)

        return SM(r"\s*time in triples corr",
                  name="CCSD(T) correction",
                  subMatchers=[
                      self.__context["method"].build_correlation_energy_matcher()
                  ],
                  onOpen={None: setup},
                  onClose={None: self.__close_auxiliary_sections}
                  )