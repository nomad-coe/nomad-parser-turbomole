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

"""This module constructs the parser for the ESCF module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
import turbomoleparser.TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")

class ESCFparser(object):

    def __init__(self, context, key="escf"):
        context[key] = self
        self.__context = context
        self.__backend = None

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
        header = SM(r"\s*e s c f\s*$",
                    name="Credits",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        sub_matchers = [
            header,
            self.__context["geo"].build_qm_geometry_matcher(),
            self.__context["geo"].build_orbital_basis_matcher(),
            self.__build_gw_matcher(),
        ]

        return self.__context.build_module_matcher("escf", sub_matchers, "ESCF")

    def __build_gw_matcher(self):
        def get_gw_approximation(backend, groups):
            types = groups[0]
            regex = re.compile(r"\s*(?P<index>[0-9]+)\s*:\s*(?P<name>[^\s:]+)")
            approximations = dict()
            match = regex.match(types)
            while match:
                approximations[match.group("index")] = match.group("name")
                match = regex.match(types, pos=match.end(2))
            backend.addValue("electronic_structure_method", "G0W0")
            backend.addValue("calculation_method_kind", "perturbative")
            backend.addValue("x_turbomole_gw_approximation", approximations[groups[1]])

        params = SM(name="GW parameters",
                    startReStr="\s*par[ae]meters:",  # typo in Turbomole 6.6 output
                    sections=["section_method"],
                    fixedStartValues={"number_of_eigenvalues_kpoints": 1},  # TM has no periodic GW
                    subMatchers=[
                        SM(r"\s*number of levels to calculate\s+(?P<number_of_eigenvalues>[0-9]+)",
                           name="num states"),
                        SM(r"\s*number of spin channels\s+(?P<number_of_spin_channels>[0-9]+)",
                           name="num spin channels"),
                        SM(r"\s*type of gw((?:\s+[0-9]+\s*:\s*\S+)+)\s+([0-9]+)\s*$",
                           name="GW approximation",
                           startReAction=get_gw_approximation),
                        SM(r"\s*rpa response function\s+(?P<x_turbomole_gw_use_rpa_response>[TF])",
                           name="GW screened interaction"),
                        SM(r"\s*eta \(Hartree\)\s+(?P<x_turbomole_gw_eta_factor__hartree>"
                           r"[+-]?[0-9]+.?[0-9]*)",
                           name="GW eta factor")
                    ]
                    )

        return SM(name="GW",
                  startReStr=r"\s*GW version\s+[0-9]+",
                  # startReStr = r"\s*GW version\s+(?P<x_turbomole_version_GW>[0-9]+)",
                  sections=["section_eigenvalues"],
                  subMatchers=[
                      params,
                      self.__build_gw_qp_states_matcher_no_spin()
                  ]
                  )

    def __build_gw_qp_states_matcher_no_spin(self):
        state = SM(r"\s*[0-9]+\s+(?P<x_turbomole_eigenvalue_ks_GroundState__eV>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+"
                   r"(?P<x_turbomole_eigenvalue_quasiParticle_energy__eV>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+"
                   r"(?P<x_turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW__eV>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+"
                   r"(?P<x_turbomole_eigenvalue_ExactExchange_perturbativeGW__eV>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+"
                   r"(?P<x_turbomole_eigenvalue_correlation_perturbativeGW__eV>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+"
                   r"(?P<x_turbomole_eigenvalue_ks_ExchangeCorrelation__eV>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+"
                   r"(?P<x_turbomole_Z_factor>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+"
                   r"(?P<x_turbomole_ExchangeCorrelation_perturbativeGW_derivation>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)",
                   name="GW QP state",
                   repeats=True)
        return SM(r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
                  name="GW QP statelist",
                  sections=["x_turbomole_section_eigenvalues_GW"],
                  subMatchers=[
                      SM(r"\s*in\s*eV", required=True, name="GW output unit"),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      state,
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      state.copy(),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                  ])
