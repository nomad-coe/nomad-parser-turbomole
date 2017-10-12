"""This module constructs the parser for the ESCF module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
import TurbomoleCommon as common

logger = logging.getLogger("nomad.turbomoleParser")


def build_escf_parser():
    return SM(name="ESCF module",
              startReStr=r"\s*escf\s*\([a-zA-Z0-9.]+\)\s+\: TURBOMOLE [a-zA-Z0-9.]+",
              subMatchers=[
                  common.build_credits_matcher("e s c f"),
                  common.build_geometry_matcher(),
                  common.build_controlinout_matcher(),
                  build_gw_matcher()
              ]
              )


def build_gw_matcher():
    def get_gw_approximation(parser):
        regex = re.compile(r"\s*type of gw(?P<types>(\s+[0-9]+\s*:\s*\S+)+)\s+(?P<id>[0-9]+)\s*$")
        match = regex.match(parser.fIn.readline())
        # id = match.group("id")
        regex = re.compile(r"\s*(?P<index>[0-9]+)\s*:\s*(?P<name>\S+)")
        types = match.group("types")
        approximations = dict()
        match = regex.match(types)
        while match:
            approximations[match.group("index")] = match.group("name")
            match = regex.match(types, pos=match.end(2))
        parser.backend.addValue("electronic_structure_method", "G0W0")
        parser.backend.addValue("calculation_method_kind", "perturbative")
        # TODO: capture implementation-specific parameters in a subsection
        # parser.backend.addValue("x_turbomole_GW_approximation", approximations[id])

    params = SM(name="GW parameters",
                startReStr="\s*par[ae]meters:",  # typo in Turbomole 6.6 output
                sections=["section_method"],
                fixedStartValues={"number_of_eigenvalues_kpoints": 1},  # no periodic GW available
                subMatchers=[
                    SM(r"\s*number of levels to calculate\s+(?P<number_of_eigenvalues>[0-9]+)",
                       name="num states"),
                    SM(r"\s*number of spin channels\s+(?P<number_of_spin_channels>[0-9]+)",
                       name="num spin channels"),
                    SM(r"\s*type of gw((\s+[0-9]+\s*:\s*\S+)+)\s+([0-9]+)\s*$",
                       name="GW approximation",
                       forwardMatch=True,
                       adHoc=get_gw_approximation)
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
