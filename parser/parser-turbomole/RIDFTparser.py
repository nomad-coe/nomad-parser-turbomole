"""This module constructs the parser for the RIDFT module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT, build_total_energy_matcher

logger = logging.getLogger("nomad.turbomoleParser")


class RIDFTparser(object):

    def __init__(self, context, key="ridft"):
        context[key] = self
        self.__context = context
        self.__backend = None

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{5,}[^+ ].+$",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(name="Credits",
                    startReStr=r"\s*DFT program with RI approximation\s*$",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        return SM(name="RIDFT module",
                  startReStr=r"\s*r i d f t\s*$",
                  sections=["section_single_configuration_calculation"],
                  subMatchers=[
                      header,
                      self.__context["method"].build_uhf_matcher(),
                      self.__context["geo"].build_qm_geometry_matcher(),
                      self.__context["geo"].build_orbital_basis_matcher(),
                      self.__context["orbitals"].build_state_matcher(),
                      self.__context["method"].build_dft_functional_matcher(),
                      self.__build_scf_cycle_matcher(),
                      build_total_energy_matcher()
                  ]
                  )

    def __build_scf_cycle_matcher(self):

        def finalize_system_data(backend, groups):
            """close the system-related sections, add HF method if no DFT usage was specified and
            link the just opened single_configuration section to the method and system sections"""
            self.__context["geo"].finalize_sections()
            self.__context["geo"].write_basis_set_mapping()
            if self.__context["method"].index_method() == -1:
                self.__context["method"].add_default_functional()
            self.__backend.addValue("single_configuration_to_calculation_method_ref",
                                    self.__context["method"].index_method())
            self.__backend.addValue("single_configuration_calculation_to_system_ref",
                                    self.__context["geo"].index_qm_geo())
            self.__context["method"].close_method_section()

        class PreviousCycle(object):
            energy = None

        def compute_energy_difference(backend, groups):
            backend.addRealValue("energy_total_scf_iteration", float(groups[0]), unit="hartree")
            if PreviousCycle.energy:
                # TODO: doublecheck we're dealing with Hartrees here
                backend.addRealValue("energy_change_scf_iteration",
                                 float(groups[0]) - PreviousCycle.energy, unit="hartree")
            PreviousCycle.energy = float(groups[0])

        total_energy = SM(r"\s*[0-9]+\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")"
                          "\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")",
                          startReAction=compute_energy_difference,
                          required=True
                          )

        # TODO: verify if coulomb energy is really a derived meta-data
        xc_energy = SM(r"\s*Exc =\s*(?P<energy_XC_scf_iteration__hartree>"+RE_FLOAT+")"
                       r"\s+Coul =\s*("+RE_FLOAT+")",
                       )
        damping = SM(r"\s*current damping\s*=\s*"+RE_FLOAT)

        scf_iteration = SM("\s*ITERATION\s+ENERGY\s+1e-ENERGY\s+2e-ENERGY\s+"
                           "RMS\[dD\(SAO\)\]\s+TOL",
                           name="SCF iteration",
                           repeats=True,
                           sections=["section_scf_iteration"],
                           subMatchers=[
                               total_energy,
                               xc_energy,
                               damping
                           ]
                           )

        return SM(r"\s*Starting SCF iterations\s*$",
                  name="HF/DFT SCF",
                  required=True,
                  subMatchers=[
                      scf_iteration
                  ],
                  startReAction=finalize_system_data
                  )
