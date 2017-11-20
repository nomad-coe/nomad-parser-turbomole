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

        sub_matchers = [
            self.__context.build_start_time_matcher(),
            header,
            self.__context["geo"].build_qm_geometry_matcher(),
            self.__context["geo"].build_orbital_basis_matcher(),
            self.__context["orbitals"].build_ir_rep_matcher(),
            self.__context["method"].build_dft_functional_matcher(),
            self.__context["method"].build_dftd3_vdw_matcher(),
            self.__context["gradient"].build_gradient_matcher(),
            Common.build_profiling_matcher(r"\s*grad(?:\.all)? profiling\s*$"),
            self.__context.build_end_time_matcher("grad")
        ]

        return self.__context.build_module_matcher("grad", sub_matchers, "GRAD",
                                                   self.__context["method"].add_default_functional)
