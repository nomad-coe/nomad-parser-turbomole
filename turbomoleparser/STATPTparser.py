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

"""This module constructs the parser for the STATPT module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from turbomoleparser.TurbomoleCommon import RE_FLOAT
import turbomoleparser.TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


# TODO: include support for transition state search once sample files are provided
class STATPTparser(object):

    def __init__(self, context, key="statpt"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__trustregion = False
        self.__geo_converged = None

    def purge_data(self):
        self.__trustregion = False
        self.__geo_converged = None

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{15,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(name="Credits",
                    startReStr=r"\s*this is S T A T P T\s*$",
                    coverageIgnore=True,
                    subMatchers=[references],
                    )

        # TODO: add geometry and force extractors for STATPT format
        sub_matchers = [
            self.__context.build_start_time_matcher(),
            header,
            self.build_hessian_update_matcher(),
            self.__context.build_end_time_matcher("statpt")
        ]

        return self.__context.build_module_matcher("statpt", sub_matchers, "STATPT")

    def build_hessian_update_matcher(self):
        def set_trust_region(backend, groups):
            self.__trustregion = True

        def set_hessian_update_method(backend, groups):
            if groups[0] == "BFGS":
                method = "trm_bfgs" if self.__trustregion else "bfgs"
                backend.addValue("geometry_optimization_method", method)

        trust_region_max = SM("\s*Maximum\s+allowed\s+trust\s+radius\s*:\s*"
                              "(?P<x_turbomole_geometry_optimization_trustregion_max__bohr>"
                              + RE_FLOAT+")",
                              name="trust region",
                              startReAction=set_trust_region)
        trust_region_min = SM("\s*Minimum\s+allowed\s+trust\s+radius\s*:\s*"
                              "(?P<x_turbomole_geometry_optimization_trustregion_min__bohr>"
                              + RE_FLOAT+")",
                              name="trust region",
                              startReAction=set_trust_region)
        trust_region_initial = SM("\s*Maximum\s+allowed\s+trust\s+radius\s*:\s*"
                                  "(?P<x_turbomole_geometry_optimization_trustregion_initial__bohr>"
                                  + RE_FLOAT+")",
                                  name="trust region",
                                  startReAction=set_trust_region)
        hessian_update = SM(r"\s*Hessian\s+update\s+method\s*:\s*([^\s].*[^\s]|[^\s])\s*$",
                            name="Hessian update",
                            startReAction=set_hessian_update_method
                            )

        def set_sampling_data(backend, gIndex, section):
            status = all(self.__geo_converged) if self.__geo_converged else None
            self.__context.set_sampling_mode_section(gIndex, status)

        cycle_index = SM(r"\s*ENERGY\s*=\s*"+RE_FLOAT+"\s*a.u.\s*;\s*#\s*of\s*cycle\s*=\s*"
                         r"(?P<x_turbomole_geometry_optimization_cycle_index>[0-9]+)\s*$",
                         name="cycle index")

        return SM(r"\s*\*+\s*Stationary\s*point\s*options\s*\*+\s*$",
                  name="Geo Opt Config",
                  sections=["section_sampling_method"],
                  onClose={"section_sampling_method": set_sampling_data},
                  fixedStartValues={"sampling_method": "geometry optimization"},
                  subMatchers=[
                      trust_region_max,
                      trust_region_min,
                      trust_region_initial,
                      hessian_update,
                      self.build_criteria_matcher(),
                      self.__context["geo"].build_qm_geometry_matcher_statpt(),
                      self.__context["gradient"].build_gradient_matcher_statpt(),
                      cycle_index,
                      self.build_convergence_matcher()
                  ]
                  )

    def build_criteria_matcher(self):
        conv_energy = SM(r"\s*Threshold\s+for\s+energy\s+change\s*:\s*"
                         r"(?P<geometry_optimization_energy_change__hartree>"+RE_FLOAT+")\s*$",
                         name="energy max")
        conv_displace = SM(r"\s*Threshold\s+for\s+max\s+displacement\s+element\s*:\s*"
                           r"(?P<geometry_optimization_geometry_change__bohr>"+RE_FLOAT+")\s*$",
                           name="geometry max")
        conv_force = SM(r"\s*Threshold\s+for\s+max\s+gradient\s+element\s*:\s*"
                        r"(?P<geometry_optimization_threshold_force__forceAu>"
                        + RE_FLOAT+")\s*$",
                        name="force max")
        conv_displace_rms = SM(r"\s*Threshold\s+for\s+RMS\s+of\s+displacement\s*:\s*"
                               r"(?P<x_turbomole_geometry_optimization_geometry_change_rms__bohr>"
                               + RE_FLOAT+")\s*$",
                               name="geometry RMS")
        conv_force_rms = SM(r"\s*Threshold\s+for\s+RMS\s+of\s+gradient\s*:\s*"
                            r"(?P<x_turbomole_geometry_optimization_threshold_force_rms__forceAu>"
                            + RE_FLOAT+")\s*$",
                            name="force RMS")
        return SM(r"\s*\*+\s*Convergence\s+criteria\s*\*+\s*$",
                  name="convergence criteria",
                  subMatchers=[
                      conv_energy,
                      conv_displace,
                      conv_force,
                      conv_displace_rms,
                      conv_force_rms
                  ]
                  )

    def build_convergence_matcher(self):
        def set_initial_convergence(backend, groups):
            self.__geo_converged = list()

        def update_convergence(backend, groups):
            self.__geo_converged.append(groups[0] == "yes")

        def build_check_matcher(criterion):
            return SM(r"\s*"+criterion+"\s+(yes|no)\s+"+RE_FLOAT+"\s+"+RE_FLOAT+"\s*$",
                      name="converged?",
                      startReAction=update_convergence
                      )

        return SM(r"\s*CONVERGENCE\s+INFORMATION\s*$",
                  name="Convergence Test",
                  startReAction=set_initial_convergence,
                  subMatchers=[
                      SM(r"\s*Converged\?\s+Value\s+Criterion\s*$", name="header"),
                      build_check_matcher(r"Energy\s+change"),
                      build_check_matcher(r"RMS\s+of\s+displacement"),
                      build_check_matcher(r"RMS\s+of\s+gradient"),
                      build_check_matcher(r"MAX\s+displacement"),
                      build_check_matcher(r"MAX\s+gradient"),
                  ]
                  )
