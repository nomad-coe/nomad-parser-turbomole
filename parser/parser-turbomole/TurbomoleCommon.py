import logging

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.simple_parser import SimpleMatcher as SM

logger = logging.getLogger("nomad.turbomoleParser")

############################################################
# This file contains functions that are needed
# by more than one parser.
############################################################

RE_FLOAT = r"(?:[+-]?(?:[0-9]+.?[0-9]*|[0-9]*.[0-9]+)(?:[DEde][+-]?[0-9]+)?)"
RE_DATE = r"(?:[0-9]{4}-[0-9]{2}-[0-9]{2})"
RE_TIME = r"(?:(?:[01][0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9]\.[0-9]{3})"


def build_total_energy_matcher():
    def set_current_energy(backend, groups):
        backend.addRealValue("energy_current", float(groups[0]), unit="hartree")

    energy_total = SM(r"\s*\|\s*total energy\s*=\s*(?P<energy_total__hartree>"
                      +RE_FLOAT+")\s*\|",
                      name="total energy",
                      required=True,
                      startReAction=set_current_energy
                      )
    energy_kinetic = SM(r"\s*:\s*kinetic energy\s*=\s*(?P<electronic_kinetic_energy__hartree>"
                        + RE_FLOAT+")\s*:\s*$",
                        name="kinetic energy",
                        required=True
                        )
    energy_potential = SM(r"\s*:\s*potential energy\s*=\s*"
                          r"(?P<x_turbomole_potential_energy_final__hartree>"+RE_FLOAT+")\s*:\s*$",
                          name="potential energy",
                          required=True
                          )
    virial_theorem = SM(r"\s*\:\s*virial theorem\s*\=\s*"
                        r"(?P<x_turbomole_virial_theorem>"+RE_FLOAT+")\s*:\s*$",
                        name="virial theorem",
                        required=True
                        )
    wavefunction_norm = SM(r"\s*\:\s*wavefunction norm\s*\=\s*"
                           r"(?P<x_turbomole_wave_func_norm>"+RE_FLOAT+")\s*:\s*$",
                           name="wavefunction norm",
                           required=True
                           )

    return SM(r"\s*convergence criteria satisfied after\s+"
              r"(?P<number_of_scf_iterations>[0-9]+)\s+iterations",
              name="SCF end",
              required=True,
              subMatchers=[
                  energy_total,
                  energy_kinetic,
                  energy_potential,
                  virial_theorem,
                  wavefunction_norm
              ]
              )


def build_profiling_matcher(title_regex):

    return SM(title_regex,
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



def build_controlinout_matcher():
    return SM (name = 'ControlInOut',
               startReStr = r"\s*\|\s*basis set information\s",
               #startReStr = r"\s*SCF run will be profiled",

               subMatchers = [
                   SM (name = 'ControlInOutLines',
                       startReStr = r"\s*we will work with the",
                       sections = ['section_topology','x_turbomole_section_functionals'],
                       weak = True,
                       subFlags = SM.SubFlags.Unordered,
                       subMatchers = [
                           # Now follows the list to match the aims output from the parsed control.in.
                           # The search is done unordered since the output is not in a specific order.
                           # Repating occurrences of the same keywords are captured.
                           # List the matchers in alphabetical order according to metadata name.
                           #
                           SM (name = 'Basis set informations',
                               startReStr = r"\s*type   atoms  prim   cont   basis",
                               #repeats = True,
                               sections = ['section_basis_set'],
                               subMatchers = [
                                   # SM (r"\s*-{20}-*", weak = True),
                                   SM (r"\s*(?P<x_turbomole_controlInOut_atom_labels>[a-zA-Z]+)\s*[0-9]+\s*(?P<x_turbomole_controlInOut_basis_prim_number>[0-9]+)\s*(?P<x_turbomole_controlInOut_basis_cont_number>[0-9]+)\s*(?P<x_turbomole_controlInOut_basis_type>[a-zA-Z-a-zA-Z]+)"
                                       ,repeats = True)
                               ]),
                           # only the first character is important for aims
                           SM (r"\s*total number of primitive shells\s*:\s*(?P<x_turbomole_controlInOut_tot_primitive_shells>[0-9]+)",sections = ['section_basis_set'], repeats = True),
                           SM (r"\s*total number of contracted shells\s*:\s*(?P<x_turbomole_controlInOut_tot_contracted_shells>[0-9]+)",sections = ['section_basis_set'], repeats = True),
                           SM (r"\s*total number of cartesian basis functions\s*:\s*(?P<x_turbomole_controlInOut_tot_cartesian_func>[0-9]+)",sections = ['section_basis_set'], repeats = True),
                           SM (r"\s*total number of SCF-basis functions\s*:\s*(?P<x_turbomole_controlInOut_tot_scf_basis_func>[0-9]+)",sections = ['section_basis_set'], repeats = True),
                           SM (r"\s*density functional"), # XC functional matching follows for turbomole_section_functionals
                           SM (r"\s*\+------------------\+\s*"),
                           SM (r"\s*(?P<x_turbomole_XC_functional_type>[a-zA-Z-a-zA-Z0-9]+)\s*(?: functional)"),
                           SM (r"\s*(?P<x_turbomole_XC_functional_type>[a-zA-Z-a-zA-Z0-9]+)\s*(?: meta-GGA functional\s)"),
                           SM (r"(?:[a-zA-Z-a-zA-Z0-9\s]+)\s*functional\:\s*(?P<x_turbomole_XC_functional_type>[a-zA-Z-a-zA-Z0-9]+)"),
                           SM (r"\s*exchange:\s*(?P<x_turbomole_controlInOut_functional_type_exchange>[a-zA-Z-+a-zA-Z0-9\(\)\s.\*]+)"),
                           SM (r"\s*correlation:\s*(?P<x_turbomole_controlInOut_functional_type_correlation>[a-zA-Z-+a-zA-Z0-9\(\)\s.\*]+)"),
                           SM (r"\s*spherical integration\s*:\s*(?P<x_turbomole_controlInOut_grid_integration>[a-zA-Z\'\s]+)"),
                           SM (r"\s*spherical gridsize\s*:\s*(?P<x_turbomole_controlInOut_grid_size>[0-9]+)"),
                           SM (r"\s*i\.e\. gridpoints\s*:\s*(?P<x_turbomole_controlInOut_grid_points_number>[0-9]+)"),
                           SM (r"\s*radial integration\s*:\s*(?P<x_turbomole_controlInOut_grid_radial_integration>[a-zA-Z0-9\(\)\s]+)"),
                           SM (r"\s*radial gridsize\s*:\s*(?P<x_turbomole_controlInOut_grid_radial_grid_size>[0-9]+)"),
                           SM (r"\s*integration cells\s*:\s*(?P<x_turbomole_controlInOut_grid_integration_cells>[0-9]+)"),
                           SM (r"\s*partition function\s*:\s*(?P<x_turbomole_controlInOut_grid_partition_func>[a-zA-Z]+)"),
                           SM (r"\s*partition sharpness\s*:\s*(?P<x_turbomole_controlInOut_grid_partition_sharpness>[0-9]+)"),
                       ]), # END ControlInOutLines
                   SM (name = 'post-HF',
                       startReStr = r"\s*(?:[a-zA-Z-a-zA-Z0-9\s]+)\s*shell calculation for the wavefunction models",
                       #sections = ['section_method'],
                       subMatchers = [
                           SM (r"\s*(?P<electronic_structure_method>[a-zA-Z-a-zA-Z0-9\(\)]+)\s*\-")
                       ]),
                   #        SM (name = "smearing",
                   #            startReStr = r"\s*and increment of one",
                   #            #sections = ["section_method"],
                   #            subMatchers = [
                   #                SM (r"\s*(?P<smearing_kind>[a-zA-Z]+)\s*smearing switched on"),
                   #                SM (r"\s*Final electron temperature\:\s*(?P<smearing_width>[0-9.eEdD]+)")
                   #            ])
               ])


def get_metaInfo(filePath):
    """Loads metadata.

    Args:
        filePath: Location of metadata.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    meta_info_env, warnings = loadJsonFile(filePath=filePath,
                                           dependencyLoader=None,
                                           extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS,
                                           uri=None)
    return meta_info_env
