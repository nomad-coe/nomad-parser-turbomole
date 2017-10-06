"""This module constructs the parser for the ESCF module from TurboMole"""

from nomadcore.simple_parser import SimpleMatcher as SM
import TurbomoleCommon as common

def build_escf_parser():
    return SM(name = "ESCF module",
              startReStr = r"\s*escf\s*\([a-zA-Z0-9.]+\)\s+\: TURBOMOLE [a-zA-Z0-9.]+",
              subMatchers = [
                  common.build_geometry_matcher(),
                  common.build_controlinout_matcher(),
                  build_gw_matcher()
              ]
    )

def build_gw_matcher():
    params = SM (name = 'GW parameters',
                 startReStr = "\s*par[ae]meters:", #typo in Turbomole 6.6 output
                 sections = ["section_method"],
                 fixedStartValues = {"number_of_eigenvalues_kpoints": 1}, #no periodic GW available
                 subMatchers = [
                     SM(r"\s*number of levels to calculate\s+(?P<number_of_eigenvalues>[0-9]+)",
                        name = "num states"),
                     SM(r"\s*number of spin channels\s+(?P<number_of_spin_channels>[0-9]+)",
                        name = "num spin channels")
                 ]
                 )

    return SM(name = 'GW',
              startReStr = r"\s*GW version\s+[0-9]+",
              # startReStr = r"\s*GW version\s+(?P<x_turbomole_version_GW>[0-9]+)",
              sections = ["section_single_configuration_calculation", "section_eigenvalues"],
              subMatchers = [
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
               name = "GW QP state",
               repeats = True)
    return SM(name = 'GW QP statelist',
              startReStr = r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
              sections = ["x_turbomole_section_eigenvalues_GW"],
              subMatchers = [
                  SM (r"\s*in\s*eV", required=True, name="GW output unit"),
                  SM (r"\s*----*", name = "<format>"),
                  state,
                  SM (r"\s*----*", name = "<format>"),
                  state.copy(),
              ])
