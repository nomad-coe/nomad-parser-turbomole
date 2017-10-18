"""This module constructs the parser for the DSCF module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")

class DSCFparser(object):

    def __init__(self, context, key="dscf"):
        context[key] = self
        self.__context = context
        self.__index_map = dict()
        self.__backend = None

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{5,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(name="Credits",
                    startReStr=r"\s*idea & directorship : reinhart ahlrichs",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        return SM(name="DSCF module",
                  startReStr=r"\s*d s c f - program",
                  subMatchers=[
                      header,
                      self.__context["geo"].build_qm_geometry_matcher(),
                      self.__context["geo"].build_orbital_basis_matcher(),
                      # TODO: read optional DFT functional specification
                      self.__build_scf_matcher()
                  ]
                  )

    def __build_scf_matcher(self):

        def finalize_system_data(backend, groups):
            """close the system-related sections, add HF method if no DFT usage was specified and
            link the just opened single_configuration section to the method and system sections"""
            self.__index_map.update(self.__context["geo"].finalize_sections())
            self.__context["geo"].write_basis_set_mapping()

            if "method" not in self.__index_map:
                self.__index_map["method"] = self.__backend.openSection("section_method")
                self.__backend.addValue("electronic_structure_method", "HF")
                self.__backend.addValue("calculation_method_kind", "absolute")
                self.__backend.closeSection("section_method", self.__index_map["method"])
            self.__backend.addValue("single_configuration_to_calculation_method_ref",
                                    self.__index_map["method"])
            self.__backend.addValue("single_configuration_calculation_to_system_ref",
                                    self.__index_map["qm-geo"])

        scf_iteration = SM("\s*current damping\s*:\s*[+-]?[0-9]+\.?[0-9]*",
                           name="SCF iteration",
                           repeats=True,
                           )

        return SM(r"\s*STARTING INTEGRAL EVALUATION FOR 1st SCF ITERATION",
                  name="HF/DFT SCF",
                  sections=["section_single_configuration_calculation"],
                  subMatchers=[
                      SM("\s*time elapsed for pre-SCF steps : cpu\s+([0-9]+\.[0-9]+)\s+sec",
                         name="SCF preparation",
                         required=True),
                      SM("\s*wall\s+([0-9]+\.[0-9]+)\s+sec",
                         name="SCF preparation",
                         required=True),
                      scf_iteration
                  ],
                  startReAction=finalize_system_data
                  )

