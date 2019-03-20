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

"""This module constructs the parser for the RIRPA module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from turbomoleparser.TurbomoleCommon import RE_FLOAT
import turbomoleparser.TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class RIRPAparser(object):

    def __init__(self, context, key="rirpa"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__no_hxx = False

    def purge_data(self):
        self.__no_hxx = False

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{5,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(r"\s*\*+\s*PROGRAM\s+RIRPA\s*\*+\s*$",
                    name="header",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        sub_matchers = [
            header,
            self.__context["geo"].build_qm_geometry_matcher(),
            self.__context["geo"].build_orbital_basis_matcher(
                "\s*\*\s*BASIS\s+SET\s+information:\s*$"),
            self.__context["geo"].build_auxiliary_basis_matcher(
                "\s*\*\s*AUXILIARY\s+BASIS\s+SET\s+information:\s*$"),
            # FIXME: figure out the difference between auxbasis and RI-J auxbasis
            self.__context["geo"].build_auxiliary_basis_matcher(
                "\s*\*\s*RIJ\s+AUXILIARY\s+BASIS\s+SET\s+information:\s*$"),
            self.__build_rpa_total_energy_matcher(),
            self.__context["gradient"].build_gradient_matcher(),
            Common.build_profiling_matcher(r"\s*Rirpa\s+profiling\s+cpu\s+wall\s+ratio\s*$"),
        ]

        return self.__context.build_module_matcher("rirpa", sub_matchers, "RIRPA")

    def __build_no_hxx_matcher(self):
        def set_flag(backend, groups):
            self.__no_hxx = True

        return SM(r"\s*The\s+HXX\s+energy\s+will\s+NOT\s+be\s+computed\s+due\s+to\s+nohxx"
                  r"\s+option.\s*$",
                  name="nohxx flag",
                  startReAction=set_flag
                  )

    def __build_rpa_total_energy_matcher(self):
        def get_total_energy(backend, groups):
            if not self.__no_hxx:
                backend.addRealValue("energy_total", float(groups[0]),
                                     self.__context.index_configuration(), unit="hartree")
            backend.addValue("electronic_structure_method", "RPA", self.__context.index_method())

        def get_correlation_energy(backend, groups):
            backend.addRealValue("energy_current", float(groups[0]),
                                 self.__context.index_configuration(), unit="hartree")

        def get_hf_energy(backend, groups):
            e_hf = float(groups[0])
            if abs(e_hf) < 1E-10:
                return
            index_config = backend.openSection("section_single_configuration_calculation")
            index_method = backend.openSection("section_method")
            references = {"section_method": index_method}
            backend.addValue("single_configuration_to_calculation_method_ref",
                             index_method, index_config)
            backend.addValue("single_configuration_calculation_to_system_ref",
                             self.__context.index_system(), index_config)
            backend.addValue("electronic_structure_method", "DFT", index_method)
            backend.addValue("calculation_method_kind", "absolute", index_method)
            index_xc = backend.openSection("section_XC_functionals")
            backend.setSectionInfo("section_XC_functionals", index_xc, references)
            backend.addValue('XC_functional_name', "HF_X", index_xc)
            backend.closeSection("section_XC_functionals", index_xc)
            backend.addRealValue("energy_total", e_hf, index_method, unit="hartree")

            index_link = backend.openSection("section_calculation_to_calculation_refs")
            references = {"section_single_configuration_calculation":
                              self.__context.index_configuration()}
            backend.setSectionInfo("section_calculation_to_calculation_refs", index_link,
                                   references)
            backend.addValue("calculation_to_calculation_kind", "source_calculation", index_link)
            backend.addValue("calculation_to_calculation_ref", index_config, index_link)
            backend.closeSection("section_calculation_to_calculation_refs", index_link)

            backend.closeSection("section_single_configuration_calculation", index_config)
            backend.closeSection("section_method", index_method)

        total_energy = SM(r"\s*\|\s*HXX\+RIRPA\s+total\s+energy\s*=\s*("+RE_FLOAT+")\s*\|\s*$",
                          name="RPA total energy",
                          startReAction=get_total_energy
                          )
        correlation = SM(r"\s*:\s*RIRPA\s+correlation\s+energy\s*=\s*("+RE_FLOAT+")\s*:\s*$",
                         name="RPA correlation ",
                         startReAction=get_correlation_energy
                         )
        hf_energy = SM(r"\s*:\s*HXX\s+total\s+energy\s*=\s*("+RE_FLOAT+")\s*:\s*$",
                       name="RPA total energy",
                       startReAction=get_hf_energy
                       )

        return SM(r"\s*Complex\s+frequency\s+integration\s*$",
                  subMatchers=[
                      total_energy,
                      correlation,
                      hf_energy
                  ]
                  )
