"""This module constructs the parser for the DSCF module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class DSCFparser(object):

    def __init__(self, context, key="dscf"):
        context[key] = self
        self.__context = context
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
                    startReStr=r"\s*d s c f - program\s*$",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        return SM(self.__context.get_module_invocation("dscf"),
                  name="DSCF module",
                  sections=["section_single_configuration_calculation"],
                  startReAction=self.__context.process_module_invocation,
                  subMatchers=[
                      self.__context.build_start_time_matcher(),
                      header,
                      self.__context["method"].build_uhf_matcher(),
                      self.__context["geo"].build_qm_geometry_matcher(),
                      self.__context["geo"].build_orbital_basis_matcher(),
                      self.__context["orbitals"].build_ir_rep_matcher(),
                      self.__context["method"].build_dft_functional_matcher(),
                      self.__context["geo"].build_embedding_matcher(),
                      self.__build_scf_cycle_matcher(),
                      Common.build_total_energy_matcher(),
                      self.__context["orbitals"].build_eigenstate_matcher(),
                      self.__build_profiling_matcher(),
                      self.__context.build_end_time_matcher("dscf")
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
            backend.addRealValue("x_turbomole_energy_1electron_scf_iteration",
                                 float(groups[1]), unit="hartree")
            backend.addRealValue("x_turbomole_energy_2electron_scf_iteration",
                                 float(groups[2]), unit="hartree")
            # TODO: clarify the meaning of the "NORM[dD(SAO)]" and "TOL" columns in the output
            if PreviousCycle.energy:
                backend.addRealValue("energy_change_scf_iteration",
                                 float(groups[0]) - PreviousCycle.energy, unit="hartree")
            PreviousCycle.energy = float(groups[0])

        def store_damping(backend, groups):
            backend.addRealValue("x_turbomole_damping_scf_iteration",
                                 1.0 / (1.0 + float(groups[0])))

        total_energy_matcher = SM(r"\s*[0-9]+\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")"
                                  "\s+("+RE_FLOAT+")\s+("+RE_FLOAT+")",
                                  startReAction=compute_energy_difference,
                                  name="SCF E total",
                                  required=True
                                  )
        xc_energy_matcher = SM(r"\s*Exc =\s*(?P<energy_XC_scf_iteration__hartree>"+RE_FLOAT+")"
                               r"\s+N =\s*("+RE_FLOAT+")",
                               name="SCF E xc"
                               )

        # TODO: identify the units of these norms
        norm_diis = SM(r"\s*Norm of current diis error:\s*(?P<x_turbomole_norm_diis_scf_iteration>"
                       + RE_FLOAT+r")\s*$",
                       name="SCF norm DIIS"
                       )
        norm_fia_block = SM(r"\s*max. resid. norm for Fia\-block\s*=\s*(?P"
                            r"<x_turbomole_norm_fia_scf_iteration>"+RE_FLOAT+")\s*for orbital\s+"
                            r"(?P<x_turbomole_norm_fia_orbital_scf_iteration>[0-9]+[a-z][0-9'\"]?"
                            r"\s*(?:alpha|beta)?)\s*$",
                            name="SCF norm Fia"
                            )
        norm_fock = SM(r"\s*max. resid. fock norm\s*=\s*(?P"
                       r"<x_turbomole_norm_fock_scf_iteration>"+RE_FLOAT+")\s*for orbital\s+"
                       r"(?P<x_turbomole_norm_fock_orbital_scf_iteration>[0-9]+[a-z][0-9'\"]?"
                       r"\s*(?:alpha|beta)?)\s*$",
                       name="SCF norm Fock"
                       )
        eigenval_change = SM(r"\s*Delta Eig\.\s*=\s*(?P<x_turbomole_delta_eigenvalues__eV>"
                             +RE_FLOAT+r")\s+eV\s*$",
                             name="SCF D eigenvals")

        scf_iteration = SM(r"\s*current damping\s*:\s*("+RE_FLOAT+")\s*$",
                           name="SCF damping",
                           repeats=True,
                           sections=["section_scf_iteration"],
                           startReAction=store_damping,
                           subMatchers=[
                               SM("\s*ITERATION\s+ENERGY\s+1e-ENERGY\s+2e-ENERGY\s+"
                                  "NORM\[dD\(SAO\)\]\s+TOL",
                                  name="SCF iteration",
                                  required=True
                                  ),
                               total_energy_matcher,
                               xc_energy_matcher,
                               norm_diis,
                               norm_fia_block,
                               norm_fock,
                               eigenval_change
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

    def __build_profiling_matcher(self):

        return SM(r"\s*dscf profiling\s*$",
                  name="profiling",
                  coverageIgnore=True,
                  subMatchers=[
                      SM(r"\s*-{20,}\s*$", name="<format>", coverageIgnore=True),
                      SM(r"\s*module\s+cpu\s+total\s+\(s\)\s+%\s+wall\s+total\s+\(s\)\s+%\s*$",
                         name="profiling",
                         coverageIgnore=True
                         ),
                      SM(r"\s*[A-z0-9\._+-]+(?: [A-z0-9\._+-]+)?(?:\s+"+RE_FLOAT+"){4}\s*$",
                         name="profiling",
                         repeats=True,
                         coverageIgnore=True
                         ),
                  ]
                  )