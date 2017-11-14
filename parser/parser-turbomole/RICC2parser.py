"""This module constructs the parser for the RICC2 module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class RICC2parser(object):

    # TODO: verify mapping of approximations to general Electronic Structure Method class
    __wavefunction_models_map = {
        "MP2": "MP2",
        "CCS": "CCS",
        "CIS": None,  # TODO: truncated CI is not listed as valid ESM on the wiki
        "CIS(D)": None,
        "CIS(Dinf)": None,
        "ADC(2)": "MP2",  # TODO: check paper to verify this mapping
        "CC2": "CCSD"
    }

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

        return SM(self.__context.get_module_invocation("ricc2"),
                  name="RICC2 module",
                  sections=["section_single_configuration_calculation"],
                  startReAction=self.__context.process_module_invocation,
                  subMatchers=[
                      self.__context.build_start_time_matcher(),
                      header,
                      self.__context["geo"].build_qm_geometry_matcher(),
                      self.__context["geo"].build_orbital_basis_matcher(),
                      self.build_wave_function_model_matcher(),
                      self.build_ref_energy_matcher(),
                      self.build_ground_state_first_order_properties_matcher(),
                      self.__context.build_end_time_matcher("ricc2")
                  ]
                  )

    def build_wave_function_model_matcher(self):

        def finalize_system_data(backend, gIndex, section):
            """link the section_single_configuration to the method and system sections"""
            self.__context["geo"].finalize_sections()
            self.__context["geo"].write_basis_set_mapping()
            backend.addValue("single_configuration_to_calculation_method_ref", gIndex)
            backend.addValue("single_configuration_calculation_to_system_ref",
                             self.__context["geo"].index_qm_geo())

        def determine_spin(backend, groups):
            if groups[0] == "restricted closed":
                pass
            elif groups[0] == "restricted open":
                self.__context["method"].spin_channels = 2
            elif groups[0] == "unrestricted open":
                self.__context["method"].spin_channels = 2
            else:
                logger.error("found unknown spin configuration in RICC2: " + groups[0])

        def extract_wf_method(backend, groups):
            method = self.__wavefunction_models_map.get(groups[0], None)
            if method:
                backend.addValue("electronic_structure_method", method)
            else:
                logger.error("unknown wave-function model encountered: %s - %s" % groups)

        method_matcher = SM(r"\s*([^\s].*[^\s])\s*-\s*([^\s].*[^\s])\s*$",
                            name="WF model",
                            required=True,
                            startReAction=extract_wf_method
                            )

        # TODO: extract further RICC2 parameters?
        return SM(r"\s*([^ ].+)\s+shell\s+calculation\s+for\s+the\s+wavefunction\s+models:\s*$",
                  startReAction=determine_spin,
                  sections=["section_method"],
                  name="spin treatment",
                  onClose={"section_method": finalize_system_data},
                  required=True,
                  subMatchers=[
                      method_matcher
                  ]
                  )

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
