"""This module constructs the parser for the ESCF module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
import TurbomoleCommon as common
from SystemParser import SystemParser

logger = logging.getLogger("nomad.turbomoleParser")


def build_escf_parser(context):
    return SM(name="ESCF module",
              startReStr=r"\s*escf\s*\([a-zA-Z0-9.]+\)\s+\: TURBOMOLE [a-zA-Z0-9.]+",
              subMatchers=[
                  common.build_credits_matcher("e s c f"),
                  context["geo"].build_qm_geometry_matcher(),
                  context["geo"].build_orbital_basis_matcher(),
                  # common.build_controlinout_matcher(),
                  build_gw_matcher(context)
              ]
              )


def build_gw_matcher(context):
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
        context["geo"].finalize_sections()
        context["geo"].write_basis_set_mapping()

    params = SM(name="GW parameters",
                startReStr="\s*par[ae]meters:",  # typo in Turbomole 6.6 output
                sections=["section_method"],
                fixedStartValues={"number_of_eigenvalues_kpoints": 1},  # no periodic GW available
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
              sections=["section_single_configuration_calculation", "section_eigenvalues"],
              subMatchers=[
                  params,
                  build_gw_qp_states_matcher_no_spin()
              ]
              )


def build_gw_qp_states_matcher_no_spin():
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
    return SM(name="GW QP statelist",
              startReStr=r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
              sections=["x_turbomole_section_eigenvalues_GW"],
              subMatchers=[
                  SM(r"\s*in\s*eV", required=True, name="GW output unit"),
                  SM(r"\s*----*", name="<format>", coverageIgnore=True),
                  state,
                  SM(r"\s*----*", name="<format>", coverageIgnore=True),
                  state.copy(),
                  SM(r"\s*----*", name="<format>", coverageIgnore=True),
              ])
