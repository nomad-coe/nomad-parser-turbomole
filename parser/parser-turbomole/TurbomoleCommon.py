from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.simple_parser import SimpleMatcher as SM

############################################################
# This file contains functions that are needed
# by more than one parser.
############################################################

def build_geometry_matcher():
    atom = SM (r"\s*(?P<x_turbomole_geometry_atom_positions_x__angstrom>[-+0-9.]+)\s+"
               "(?P<x_turbomole_geometry_atom_positions_y__angstrom>[-+0-9.]+)\s+"
               "(?P<x_turbomole_geometry_atom_positions_z__angstrom>[-+0-9.]+)\s+"
               "(?P<x_turbomole_geometry_atom_labels>[a-zA-Z]+)\s+"
               "(?P<x_turbomole_geometry_atom_charge>[0-9.]+)", repeats = True,
               name="atom data")
    return SM (name = 'Geometry',
               startReStr = r"\s*\|\s*Atomic coordinate",
               sections = ['section_system'],
               subMatchers = [
                   SM (r"\s*-{20}-*", weak = True),
                   SM (startReStr = r"\s*atomic coordinates",subMatchers = [atom])]
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
    metaInfoEnv, warnings = loadJsonFile(filePath = filePath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)
    return metaInfoEnv