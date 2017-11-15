"""This module constructs the parser for the GRAD module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class GRADparser(object):

    def __init__(self, context, key="grad"):
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
        header = SM(r"\s*g r a d - program\s*$",
                    name="header",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        return SM(self.__context.get_module_invocation("grad"),
                  name="GRAD module",
                  sections=["section_single_configuration_calculation"],
                  startReAction=self.__context.process_module_invocation,
                  subMatchers=[
                      self.__context.build_start_time_matcher(),
                      header,
                      self.__context["geo"].build_qm_geometry_matcher(),
                      self.__context["geo"].build_orbital_basis_matcher(),
                      self.__context["method"].build_dft_functional_matcher(),
                      self.build_gradient_calculation_start_matcher(),
                      self.__context["gradient"].build_gradient_matcher(),
                      self.__context.build_end_time_matcher("grad")
                  ]
                  )

    def build_gradient_calculation_start_matcher(self):

        def finalize_system_data(backend, groups):
            """link the section_single_configuration to the method and system sections"""
            self.__context["geo"].finalize_sections()
            self.__context["geo"].write_basis_set_mapping()
            backend.addValue("single_configuration_to_calculation_method_ref",
                             self.__context["method"].index_method())
            backend.addValue("single_configuration_calculation_to_system_ref",
                             self.__context["geo"].index_qm_geo())
            self.__context["method"].close_method_section()

        return SM(r"\s*SCF\s+ENERGY\s+GRADIENT\s+with\s+respect\s+to\s+NUCLEAR\s+COORDINATES\s*$",
                  startReAction=finalize_system_data,
                  name="gradient start",
                  required=True
                  )
