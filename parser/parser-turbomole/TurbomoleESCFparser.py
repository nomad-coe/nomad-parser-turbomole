"""This module constructs the parser for the ESCF module from TurboMole"""

from nomadcore.simple_parser import SimpleMatcher as SM
import TurbomoleCommon as common

def build_ESCF_parser():
    return SM(name = 'ESCF module',
              startReStr = r"\s*escf\s*\([a-zA-Z0-9.]+\)\s+\: TURBOMOLE [a-zA-Z0-9.]+",
              subMatchers = [
                  common.build_geometry_matcher(),
                  common.build_controlinout_matcher(),
                  build_gw_matcher("_perturbativeGW")
              ]
    )

def build_gw_matcher(section_suffix):
    GWEigenvaluesListSubMatcher = SM (name = 'x_turbomole_perturbativeGW_EigenvaluesLists',
                                      #	   startReStr = r"\s*in\s*eV",
                                      startReStr = r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
                                      sections = ['x_turbomole_section_eigenvalues_list%s' % section_suffix],
                                      subMatchers = [
                                          SM (r"\s*in\s*eV"),
                                          SM (r"\s*------------------------------------------------------------------------------------"),
                                          SM (r"\s*(?P<x_turbomole_eigenstate_number>[0-9]+)\s+(?P<x_turbomole_eigenvalue_ks_GroundState__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_quasiParticle_energy__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_ExactExchange_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_correlation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_ks_ExchangeCorrelation__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_Z_factor>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_ExchangeCorrelation_perturbativeGW_derivation>[-+0-9.eEdD]+)",
                                              adHoc = lambda parser: parser.superContext.setStartingPointCalculation(parser),
                                              repeats = True),
                                          SM (r"\s*------------------------------------------------------------------------------------"),
                                          SM (r"\s*(?P<x_turbomole_eigenstate_number>[0-9]+)\s+(?P<x_turbomole_eigenvalue_ks_GroundState__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_quasiParticle_energy__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_ExactExchange_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_correlation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_eigenvalue_ks_ExchangeCorrelation__eV>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_Z_factor>[-+0-9.eEdD]+)\s+"
                                              "(?P<x_turbomole_ExchangeCorrelation_perturbativeGW_derivation>[-+0-9.eEdD]+)",
                                              adHoc = lambda parser: parser.superContext.setStartingPointCalculation(parser),
                                              repeats = True)
                                      ])
    return SM (name = 'x_turbomole_perturbativeGW_EigenvaluesGroup',
               startReStr = r"\s*GW\s*version:",
               sections = ['x_turbomole_section_eigenvalues_group%s' % section_suffix],
               subMatchers = [
                   # non-spin-polarized
                   SM (name = 'x_turbomole_GW_EigenvaluesNoSpinNonPeriodic',
                       startReStr = r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
                       sections = ['x_turbomole_section_eigenvalues_spin%s' % section_suffix],
                       forwardMatch = True,
                       subMatchers = [
                           #               SM (r"\s*-+"),
                           #               SM (r"\s*-+"),
                           GWEigenvaluesListSubMatcher.copy()
                       ]), # END EigenvaluesNoSpinNonPeriodic
               ])
