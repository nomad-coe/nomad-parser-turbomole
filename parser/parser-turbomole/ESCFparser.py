"""This module constructs the parser for the ESCF module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")

class ESCFparser(object):

    def __init__(self, context, key="escf"):
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
        header = SM(r"\s*e s c f\s*$",
                    name="Credits",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        return SM(self.__context.get_module_invocation("escf"),
                  name="ESCF module",
                  startReAction=self.__context.process_module_invocation,
                  sections=["section_single_configuration_calculation"],
                  subMatchers=[
                      self.__context.build_start_time_matcher(),
                      header,
                      self.__context["geo"].build_qm_geometry_matcher(),
                      self.__context["geo"].build_orbital_basis_matcher(),
                      self.__build_gw_matcher(),
                      self.__context.build_end_time_matcher("escf")
                  ]
                  )

    def __build_gw_matcher(self):
        def get_gw_approximation(backend, groups):
            types = groups[0]
            regex = re.compile(r"\s*(?P<index>[0-9]+)\s*:\s*(?P<name>[^\s:]+)")
            approximations = dict()
            match = regex.match(types)
            while match:
                approximations[match.group("index")] = match.group("name")
                match = regex.match(types, pos=match.end(2))
            backend.addValue("electronic_structure_method", "G0W0")
            backend.addValue("calculation_method_kind", "perturbative")
            backend.addValue("x_turbomole_gw_approximation", approximations[groups[1]])

        def finalize_system_data(backend, groups):
            self.__context["geo"].finalize_sections()
            self.__context["geo"].write_basis_set_mapping()
            self.__backend.addValue("single_configuration_to_calculation_method_ref",
                                    backend.get_latest_section("section_method").gIndex)
            self.__backend.addValue("single_configuration_calculation_to_system_ref",
                                    self.__context["geo"].index_qm_geo())

        params = SM(name="GW parameters",
                    startReStr="\s*par[ae]meters:",  # typo in Turbomole 6.6 output
                    sections=["section_method"],
                    fixedStartValues={"number_of_eigenvalues_kpoints": 1},  # TM has no periodic GW
                    startReAction=finalize_system_data,
                    subMatchers=[
                        SM(r"\s*number of levels to calculate\s+(?P<number_of_eigenvalues>[0-9]+)",
                           name="num states"),
                        SM(r"\s*number of spin channels\s+(?P<number_of_spin_channels>[0-9]+)",
                           name="num spin channels"),
                        SM(r"\s*type of gw((?:\s+[0-9]+\s*:\s*\S+)+)\s+([0-9]+)\s*$",
                           name="GW approximation",
                           startReAction=get_gw_approximation),
                        SM(r"\s*rpa response function\s+(?P<x_turbomole_gw_use_rpa_response>[TF])",
                           name="GW screened interaction"),
                        SM(r"\s*eta \(Hartree\)\s+(?P<x_turbomole_gw_eta_factor__hartree>"
                           r"[+-]?[0-9]+.?[0-9]*)",
                           name="GW eta factor")
                    ]
                    )

        return SM(name="GW",
                  startReStr=r"\s*GW version\s+[0-9]+",
                  # startReStr = r"\s*GW version\s+(?P<x_turbomole_version_GW>[0-9]+)",
                  sections=["section_eigenvalues"],
                  subMatchers=[
                      params,
                      self.__build_gw_qp_states_matcher_no_spin()
                  ]
                  )

    def __build_gw_qp_states_matcher_no_spin(self):
        state = SM(r"\s*[0-9]+\s+(?P<x_turbomole_eigenvalue_ks_GroundState__eV>[-+0-9.eEdD]+)\s+"
                   r"(?P<x_turbomole_eigenvalue_quasiParticle_energy__eV>[-+0-9.eEdD]+)\s+"
                   r"(?P<x_turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                   r"(?P<x_turbomole_eigenvalue_ExactExchange_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                   r"(?P<x_turbomole_eigenvalue_correlation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                   r"(?P<x_turbomole_eigenvalue_ks_ExchangeCorrelation__eV>[-+0-9.eEdD]+)\s+"
                   r"(?P<x_turbomole_Z_factor>[-+0-9.eEdD]+)\s+"
                   r"(?P<x_turbomole_ExchangeCorrelation_perturbativeGW_derivation>[-+0-9.eEdD]+)",
                   name="GW QP state",
                   repeats=True)
        return SM(r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
                  name="GW QP statelist",
                  sections=["x_turbomole_section_eigenvalues_GW"],
                  subMatchers=[
                      SM(r"\s*in\s*eV", required=True, name="GW output unit"),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      state,
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      state.copy(),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                  ])
