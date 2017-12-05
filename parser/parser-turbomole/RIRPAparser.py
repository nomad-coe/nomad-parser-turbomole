"""This module constructs the parser for the RIRPA module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT
import TurbomoleCommon as Common

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
            self.__context.build_start_time_matcher(),
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
            self.__context.build_end_time_matcher("rirpa")
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

        total_energy = SM(r"\s*\|\s*HXX\+RIRPA\s+total\s+energy\s*=\s*("+RE_FLOAT+")\s*\|\s*$",
                          name="RPA total energy",
                          startReAction=get_total_energy
                          )
        correlation = SM(r"\s*:\s*RIRPA\s+correlation\s+energy\s*=\s*("+RE_FLOAT+")\s*:\s*$",
                         name="RPA correlation ",
                         startReAction=get_correlation_energy
                         )

        return SM(r"\s*Complex\s+frequency\s+integration\s*$",
                  subMatchers=[
                      total_energy,
                      correlation
                  ]
                  )
