"""This module constructs the parser for the RICC2 module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class RICC2parser(object):

    def __init__(self, context, key="ricc2"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__ref_energy = None

    def purge_data(self):
        self.__ref_energy = None

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{5,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(r"\s*R I C C 2 - PROGRAM\s*$",
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
            self.__context["method"].build_wave_function_model_matcher(),
            self.__context["geo"].build_auxiliary_basis_matcher(),
            self.build_ref_energy_matcher(),
            self.build_ground_state_first_order_properties_matcher(),
            self.__context["gradient"].build_gradient_matcher(),
            self.__context.build_end_time_matcher("ricc2")
        ]

        return self.__context.build_module_matcher("ricc2", sub_matchers, "RICC2")

    def build_ref_energy_matcher(self):

        def get_reference_energy(backend, groups):
            self.__ref_energy = float(groups[0])

        return SM(r"\s*Energy\s+of\s+reference\s+wave\s+function\s+is\s+("+RE_FLOAT+")\s*$",
                  name="ref WF energy",
                  startReAction=get_reference_energy,
                  required=True
                  )

    def build_ground_state_first_order_properties_matcher(self):

        def compute_energy_correction(backend, groups):
            energy_diff = float(groups[0]) - self.__ref_energy
            backend.addRealValue("energy_total", float(groups[0]), unit="eV")
            backend.addRealValue("energy_current", energy_diff, unit="eV")

        method_matcher = SM(r"\s*Method\s*:\s*([^\s].*[^\s])\s*$",
                            name="method",
                            required=True
                            )
        energy_matcher = SM(r"\s*Total\s+Energy\s*:\s*("+RE_FLOAT+")\s*$",
                            name="total energy",
                            startReAction=compute_energy_correction,
                            required=True
                            )

        return SM(r"\s*\*<+\s+GROUND\s+STATE\s+FIRST-ORDER\s+PROPERTIES\s+>+\*\s*$",
                  name="GS 1st order props",
                  subMatchers=[
                      method_matcher,
                      energy_matcher
                  ]
                  )
