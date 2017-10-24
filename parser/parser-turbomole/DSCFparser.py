"""This module constructs the parser for the DSCF module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT

logger = logging.getLogger("nomad.turbomoleParser")

class DSCFparser(object):

    def __init__(self, context, key="dscf"):
        context[key] = self
        self.__context = context
        self.__index_map = dict()
        self.__backend = None
        self.__spin_polarized = False

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
                  sections=["section_single_configuration_calculation"],
                  subMatchers=[
                      header,
                      self.__build_uhf_matcher(),
                      self.__context["geo"].build_qm_geometry_matcher(),
                      self.__context["geo"].build_orbital_basis_matcher(),
                      # TODO: read optional DFT functional specification
                      self.__build_scf_cycle_matcher(),
                      self.__build_total_energy_matcher()
                  ]
                  )

    def __build_uhf_matcher(self):

        def set_spin_polarized(backend, groups):
            self.__spin_polarized = True

        return SM(r"\s*UHF mode switched on !",
                  name="UHF switch",
                  startReAction=set_spin_polarized
                  )

    def __build_scf_cycle_matcher(self):

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

        class PreviousCycle(object):
            energy = None

        def compute_energy_difference(backend, groups):
            backend.addRealValue("energy_total_scf_iteration", float(groups[0]), unit="hartree")
            if PreviousCycle.energy:
                # TODO: doublecheck we're dealing with Hartrees here
                backend.addRealValue("energy_change_scf_iteration",
                                 float(groups[0]) - PreviousCycle.energy, unit="hartree")
            PreviousCycle.energy = float(groups[0])

        total_energy_matcher = SM(r"\s*[0-9]+\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")"
                                  "\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")",
                                  startReAction=compute_energy_difference,
                                  required=True
                                  )
        xc_energy_matcher = SM(r"\s*Exc =\s*(?P<energy_XC_scf_iteration__hartree>"+RE_FLOAT+")"
                               r"\s+N =\s*("+RE_FLOAT+")",
                               )

        scf_iteration = SM("\s*current damping\s*:\s*[+-]?[0-9]+\.?[0-9]*",
                           name="SCF iteration",
                           repeats=True,
                           sections=["section_scf_iteration"],
                           subMatchers=[
                               SM("\s*ITERATION\s+ENERGY\s+1e-ENERGY\s+2e-ENERGY\s+"
                                  "NORM\[dD\(SAO\)\]\s+TOL",
                                  name="SCF iteration",
                                  required=True
                                  ),
                               total_energy_matcher,
                               xc_energy_matcher
                           ]
                           )

        return SM(r"\s*STARTING INTEGRAL EVALUATION FOR 1st SCF ITERATION",
                  name="HF/DFT SCF",
                  required=True,
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

    def __build_total_energy_matcher(self):
        def set_current_energy(backend, groups):
            backend.addRealValue("energy_current", float(groups[0]), unit="hartree")

        energy_total = SM(r"\s*\|\s*total energy\s*=\s*(?P<energy_total__hartree>"
                          +RE_FLOAT+")\s*\|",
                          name = "total energy",
                          required=True,
                          startReAction=set_current_energy
                          )
        energy_kinetic = SM(r"\s*:\s*kinetic energy\s*=\s*(?P<electronic_kinetic_energy__hartree>"
                            + RE_FLOAT+")\s*:",
                            name="kinetic energy",
                            required=True
                            )
        energy_potential = SM(r"\s*:\s*potential energy\s*=\s*"
                              r"(?P<x_turbomole_potential_energy_final__hartree>"+RE_FLOAT+")\s*:",
                              name="potential energy",
                              required=True
                              )

        return SM(r"\s*convergence criteria satisfied after\s+"
                  r"(?P<number_of_scf_iterations>[0-9]+)\s+iterations",
                  name="SCF end",
                  required=True,
                  subMatchers=[
                      energy_total,
                      energy_kinetic,
                      energy_potential
                  ]
                  )
