

from builtins import object
import setup_paths
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import get_metaInfo
import logging, os, re, sys
from nomadcore.unit_conversion.unit_conversion import convert_unit_function
import TurbomoleCommon as common
from SystemParser import SystemParser
from OrbitalParser import OrbitalParser
from MethodParser import MethodParser
from ESCFparser import ESCFparser
from DSCFparser import DSCFparser

eV2J = convert_unit_function("eV","J")

############################################################
# This is the parser for the main file of turbomole.
############################################################

logger = logging.getLogger("nomad.turbomoleParser")


class TurbomoleParserContext(object):

    def __init__(self):
        self.__data = dict()
        self.functionals = []
        self.generic = False

    def __getitem__(self, item):
        return self.__data[item]

    def __setitem__(self, key, value):
        if key in self.__data:
            raise Exception("sub data storage '%s' already declared!" % key)
        self.__data[key] = value

    def __iter__(self):
        for key, sub_parser in self.__data.items():
            yield key, sub_parser

    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        self.secMethodIndex = None
        self.secSystemDescriptionIndex = None
        self.lastCalculationGIndex = None
        self.singleConfCalcs = []
        self.geoConvergence = None
        self.eigenvalues = []
        self.occupation = []
        self.evSymm = []
        self.alphaEv = None

    def switchEVSpins(self):
        """stores alpha spin and prepares for beta spin"""
        self.alphaEv = {
            'eigenvalues': self.eigenvalues,
            'occupation': self.occupation,
            'evSymm': self.evSymm
        }
        self.eigenvalues = []
        self.occupation = []
        self.evSymm = []

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Get compiled parser, filename and metadata.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        self.fName = fInName
        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()

    def onClose_section_eigenvalues(self, backend, gIndex, section):
        """write eigenvalues to the backend and then cleans local storage"""
        pass

    def onClose_section_run(self, backend, gIndex, section):

        if self.geoConvergence is not None:
            backend.addValue('x_turbomole_geometry_optimization_converged', self.geoConvergence)
            if self.geoConvergence is True:
                sampling_method = "geometry_optimization"
            elif len(self.singleConfCalcs) > 1:
                pass
            else:
                return
            samplingGIndex = backend.openSection("section_sampling_method")
            backend.addValue("sampling_method", sampling_method)
            backend.closeSection("section_sampling_method", samplingGIndex)
            frameSequenceGIndex = backend.openSection("section_frame_sequence")
            backend.addValue("frame_sequence_to_sampling_ref", samplingGIndex)
            backend.addArrayValues("frame_sequence_local_frames_ref", np.asarray(self.singleConfCalcs))
            backend.closeSection("section_frame_sequence", frameSequenceGIndex)

    def onOpen_section_method(self, backend, gIndex, section):
        # keep track of the latest method section
        self.secMethodIndex = gIndex

    ###################################################################
    # (3.4) onClose for geometry and force (section_system)
    # todo: maybe we can move the force to onClose_section_single_configuration_calculation in the future. 
    ###################################################################

    def onOpen_section_system(self, backend, gIndex, section):
        # keep track of the latest system description section
        self.secSystemDescriptionIndex = gIndex

    def onClose_x_turbomole_section_eigenvalues_list(self, backend, gIndex, section):
        irrep_name = section['x_turbomole_irreducible_representation_state_str']
        for item in range(len(irrep_name)):
            Irrepresent = irrep_name[item].split()
            self.evSymm += Irrepresent

        eigenvalues_name = section['x_turbomole_eigenvalue_eigenvalue_str']
        for mem in range(len(eigenvalues_name)):
            Eigenval = eigenvalues_name[mem].split()
            self.eigenvalues += map(lambda x: eV2J(float(x)), Eigenval)

        occupation_name = section['x_turbomole_eigenvalue_occupation_str']
        if not occupation_name == None:
            for ele in range(len(occupation_name)):
                Occupat = occupation_name[ele].split()
                self.occupation += map(float, Occupat)
        self.occupation += [0.0 for i in range(len(self.evSymm)-len(self.occupation))]

    def onOpen_section_single_configuration_calculation(self, backend, gIndex, section):
        self.singleConfCalcs.append(gIndex)

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        if section['x_turbomole_geometry_optimization_converged'] is not None:
            if section['x_turbomole_geometry_optimization_converged'][-1] == 'FULFILLED':
                self.geoConvergence = True
            else:
                self.geoConvergence = False

        if self.generic:
            if self.secMethodIndex:
                backend.addValue('single_configuration_to_calculation_method_ref',
                                 self.secMethodIndex)
            if self.secSystemDescriptionIndex:
                backend.addValue('single_configuration_calculation_to_system_ref',
                                 self.secSystemDescriptionIndex)

    def setStartingPointCalculation(self, parser):
        backend = parser.backend
        backend.openSection('section_calculation_to_calculation_refs')
        if self.lastCalculationGIndex:
            backend.addValue('calculation_to_calculation_ref', self.lastCalculationGIndex)
        backend.addValue('calculation_to_calculation_kind', 'pertubative GW')
        #        backend.closeSection('section_calculation_to_calculation_refs')
        return None

#############################################################
#################[2] MAIN PARSER STARTS HERE  ###############
#############################################################

def build_root_parser(context):
    """Builds the SimpleMatcher to parse the main file of turbomole.
    Matches for subsections of the output file are generated in dedicated
    functions below.

    Returns:
       SimpleMatcher that parses main file of Turbomole.
    """

    def set_backends(backend, groups):
        for key, sub_parser in context:
            sub_parser.set_backend(backend)

    def finalize_system_data(backend, groups):
        context["geo"].finalize_sections()
        context["geo"].write_basis_set_mapping()
        context["method"].close_method_section()

    # shared subparsers created here are automatically stored in the context
    SystemParser(context)
    OrbitalParser(context)
    MethodParser(context)

    def set_generic(backend, groups):
        context.generic = True

    # matches only those subprograms without dedicated parser
    generic = SM(r"\s*(?:aoforce|cosmoprep|egrad|evib|frog|gradsammel|"
                 r"hesssammel|moloch|odft|relax|ridft|rirpa|sdg|thirdsammel|"
                 r"vibration|atbandbta|define|eigerf|fdetools|grad|haga|"
                 r"intense|mpgrad|proper|ricc2|rimp2|ruecker|statpt|tm2molden|"
                 r"woelfling|bsseenergy|freeh|gradruecker|riper|"
                 r"hessruecker|mdprep|mpshift|rdgrad|ricctools|rimp2prep|"
                 r"sammler|thirdruecker|uff)\s*"
                 r"\([a-zA-Z0-9.]+\) \: TURBOMOLE [a-zA-Z0-9.]+",
                 name="NewRun",
                 startReAction=set_generic,
                 subMatchers=[
                     SM(name="general info",
                         startReStr=r"\s*Copyright \(C\) ",
                         subMatchers=[
                             context["geo"].build_qm_geometry_matcher(),
                             context["geo"].build_orbital_basis_matcher(),
                             context["method"].build_dft_functional_matcher()
                         ]),
                     # the actual section for a single configuration calculation starts here
                     SM(r"\s*1e\-*integrals will be neglected if expon",
                        name="Single Config",
                        startReAction=finalize_system_data,
                        sections=["section_single_configuration_calculation"],
                        # repeats = True,
                        subMatchers=[
                            SM(r"\s*\|\s*EMBEDDING IN PERIODIC POINT CHARGES\s*\|",
                               name = "Embedding",
                               subMatchers=[
                                   #SmearingOccupation,
                                   build_embedding_matcher()
                               ]
                               ),
                            SM(name='TotalEnergyForEachScfCycle',
                                startReStr = r"\s*scf convergence criterion",
                                subMatchers=[
                                    #SmearingOccupation,
                                    build_total_energy_scf_matcher(),
                                    build_total_energy_final_matcher()
                                ]),
                            build_eigenvalues_matcher(),
                            build_forces_matcher(),
                        ]),
                     SM(r"\s*Energy of reference wave function is",
                        name="PostHFTotalEnergies",
                        sections = ["section_single_configuration_calculation"],
                        subMatchers=[
                            build_total_energy_coupled_cluster_matcher()
                        ]
                        ),
                     SM(r"\s*\|\s*MP2 relaxed",
                        name="PTTotalEnergies",
                        sections=["section_single_configuration_calculation"],
                        subMatchers=[
                            build_total_energy_perturbation_theory_matcher()
                        ]
                        )
                 ]
                 )
    modules = [
        ESCFparser(context).build_parser(),
        DSCFparser(context).build_parser(),
        generic,
        build_relaxation_matcher()
    ]

    return SM (name = 'Root',
               startReStr = "",
               forwardMatch = True,
               weak = True,
               sections = ['section_run'],
               subMatchers = [
                   SM(name = 'ProgramHeader',
                      startReStr = r"\s*(?:aoforce|cosmoprep|egrad|evib|frog|gradsammel|"
                                   r"hesssammel|moloch|odft|relax|ridft|rirpa|sdg|thirdsammel|"
                                   r"vibration|atbandbta|define|eigerf|fdetools|grad|haga|"
                                   r"intense|mpgrad|proper|ricc2|rimp2|ruecker|statpt|tm2molden|"
                                   r"woelfling|bsseenergy|dscf|escf|freeh|gradruecker|"
                                   r"hessruecker|mdprep|mpshift|rdgrad|ricctools|rimp2prep|"
                                   r"sammler|thirdruecker|uff)\s*"
                                   r"\((?P<x_turbomole_nodename>[a-zA-Z0-9.]+)\) \: "
                                   r"TURBOMOLE (?P<program_version>[a-zA-Z0-9.]+)",
                      forwardMatch=True,  # necessary to match runs without dedicated subparser
                      startReAction=set_backends,
                      fixedStartValues={'program_name': 'turbomole', 'program_basis_set_type': 'GTOs'},
                      subMatchers=modules
                      )
               ]
               )

def build_eigenvalues_matcher():
    return SM(name = 'Eigenvalues',
              repeats = False,
              sections = ["section_eigenvalues"],
              startReStr = r"\s*orbitals .*  will be written to file.*",
              subMatchers = [
                  SM(r"\s*(?:alpha:)\s*"),
                  SM(r"\s*(?: irrep)\s*(?P<x_turbomole_irreducible_representation_state_str>[0-9a-z\s]+)",
                     repeats = True,
                     sections = ['x_turbomole_section_eigenvalues_list'],
                     subMatchers = [
                         SM(r"\s*eigenvalues H"),
                         SM (r"\s*(?: eV)\s*(?P<x_turbomole_eigenvalue_eigenvalue_str>[-+0-9a-z.eEdD\s]+)", repeats = True),
                         SM (r"\s*(?: occupation)\s*(?P<x_turbomole_eigenvalue_occupation_str>[0-9.\s]+)", repeats = True)
                     ]),
                  SM(r"\s*(?:beta:)\s*",
                     adHoc = lambda parser: parser.superContext.switchEVSpins(),
                     subMatchers = [
                         SM(r"\s*(?: irrep)\s*(?P<x_turbomole_irreducible_representation_state_str>[0-9a-z\s]+)",
                            repeats = True,
                            sections = ['x_turbomole_section_eigenvalues_list'],
                            subMatchers = [
                                SM(r"\s*eigenvalues H"),
                                SM (r"\s*(?: eV)\s*(?P<x_turbomole_eigenvalue_eigenvalue_str>[-+0-9a-z.eEdD\s]+)", repeats = True),
                                SM (r"\s*(?: occupation)\s*(?P<x_turbomole_eigenvalue_occupation_str>[0-9.\s]+)", repeats = True)
                            ])
                     ])
              ])

def build_forces_matcher():
    return SM (name = 'AtomicForces',
               repeats =True,
               startReStr = r"\s*ATOM\s*CARTESIAN GRADIENTS",
               #forwardMatch = True,
               subMatchers = [
                   SM (r"\s*(?:[0-9]+)\s*(?:[a-z]+)\s*(?P<x_turbomole_atom_forces_raw_x__hartree_bohr_1>[-+0-9.eEdD]+)\s*(?P<x_turbomole_atom_forces_raw_y__hartree_bohr_1>[-+0-9.eEdD]+)\s*(?P<x_turbomole_atom_forces_raw_z__hartree_bohr_1>[-+0-9.eEdD]+)", repeats = True)
               ])

def build_total_energy_scf_matcher():
    return SM (name = 'TotalEnergyScf',
               repeats =True,
               #startReStr = r"\s*scf convergence criterion",
               #startReStr = r"\s*current damping\s*:\s*",
               sections = ['section_scf_iteration'],
               startReStr = r"\s*ITERATION  ENERGY\s*",
               #forwardMatch = True,
               subMatchers = [
                   SM (r"\s*current damping\s*:\s*(?P<x_turbomole_energy_scf_damping>[0-9.eEdD]+)"),
                   SM (r"\s*(?P<x_turbomole_iteration_number>[0-9]+)\s*(?P<x_turbomole_energy_total_scf_iteration__eV>[-+0-9.eEdD]+)\s*(?P<x_turbomole_energy_one_scf_iteration__eV>[-+0-9.eEdD]+)"
                       "\s*(?P<x_turbomole_energy_two_scf_iteration__eV>[-+0-9.eEdD]+)\s*(?P<x_turbomole_energy_norm_scf_iteration__eV>[-+0-9.eEdD]+)\s*(?P<x_turbomole_energy_tolerance_scf_iteration__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*max. resid. norm for Fia\-block\=\s*(?P<x_turbomole_max_res_norm_fia_block>[-+0-9.eEdD]+)\s*for orbital\s*(?P<x_turbomole_orbital_name_fia_block>[a-z0-9\s]+)"),
                   SM (r"\s*max. resid. fock norm\s*\=\s*(?P<x_turbomole_max_res_norm_fock_norm>[-+0-9.eEdD]+)\s*for orbital\s*(?P<x_turbomole_orbital_name_fock_norm>[a-z0-9\s]+)"),
                   SM (r"\s*irrep a   \: virtual orbitals shifted by\s*(?P<x_turbomole_virtual_orbital_shift>[0-9.]+)"),
                   SM (r"\s*Delta Eig\.\s*\=\s*(?P<x_turbomole_delta_eigenvalues__eV>[-+0-9.eEdD]+)\s*eV")
               ])

def build_total_energy_final_matcher():
    return SM (name = 'TotalEnergyFinal',
               startReStr = r"\s*\|\s*total energy\s*\=",
               forwardMatch = True,
               subMatchers = [
                   SM (r"\s*\|\s*total energy\s*\=\s*(?P<energy_total__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\:\s*kinetic energy\s*\=\s*(?P<x_turbomole_kinetic_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\:\s*potential energy\s*\=\s*(?P<x_turbomole_potential_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\:\s*virial theorem\s*\=\s*(?P<x_turbomole_virial_theorem_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\:\s*wavefunction norm\s*\=\s*(?P<x_turbomole_wave_func_norm__eV>[-+0-9.eEdD]+)")

               ])

def build_occupation_smearing_matcher():
    return SM (name = "smearing",
               startReStr = r"\s*and increment of one",
               sections = ["section_method"],
               subMatchers = [
                   SM (r"\s*(?P<smearing_kind>[a-zA-Z]+)\s*smearing switched on"),
                   SM (r"\s*Final electron temperature\:\s*(?P<smearing_width>[0-9.eEdD]+)")
               ])

def build_total_energy_coupled_cluster_matcher():
    return SM (name = 'TotalEnergyCC',
               #startReStr = r"\s*Calculate\s*integrals\s*\(*ia\|*jb\)*\s*for MP2 start guess",
               startReStr = r"\s*=========================================================================",
               forwardMatch = True,
               subMatchers = [
                   #SM (r"\s*\*\s*RHF  energy\s*\:\s*(?P<x_turbomole_HF_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*RHF  energy\s*\:\s*(?P<energy_total__eV>[-+0-9.eEdD]+)"),
                   #SM (r"\s*\*\s*UHF  energy\s*\:\s*(?P<x_turbomole_HF_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*UHF  energy\s*\:\s*(?P<energy_total__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final MP2 energy\s*\:\s*(?P<x_turbomole_MP2_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final CCSD energy\s*\:\s*(?P<x_turbomole_CCSD_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final CC2 energy\s*\:\s*(?P<x_turbomole_CC2_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final CCSD\(T\) energy\s*\:\s*(?P<x_turbomole_CCSDparT_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*D1 diagnostic \(CCSD\)\s*\:\s*(?P<x_turbomole_D1_diagnostic>[-+0-9.eEdD]+)")

               ])

def build_total_energy_perturbation_theory_matcher():
    return SM (name = 'TotalEnergyPT',
               startReStr = r"\s*\|*\s*| natural orb",
               forwardMatch = True,
               subMatchers = [
                   SM (r"\s*Total Energy\s*\:\s*(?P<x_turbomole_PT_total_energy_final__eV>[-+0-9.eEdD]+)")

               ])

def build_embedding_matcher():
    return SM (name = 'PeriodicEmbedding',
               startReStr = r"\s*\+------------------------ Parameters ------------------------\+",
               forwardMatch = True,
               subMatchers = [
                   SM (r"\s*Maximum multipole moment used               :\s*(?P<x_turbomole_max_multipole_moment>[0-9]+)"),
                   SM (r"\s*Multipole precision parameter               :\s*(?P<x_turbomole_multipole_precision_parameter>[-+0-9.eEdD]+)"),
                   SM (r"\s*Minimum separation between cells            :\s*(?P<x_turbomole_min_separation_cells>[-+0-9.eEdD]+)"),
                   SM (r"\s*\+-----------------------------------------------------------\+\s*"),
                   SM (r"\s*Charge Neutrality tolerance :\s*(?P<x_turbomole_charge_neutrality_tol>[-+0-9.eEdD]+)"),
                   SM (r"\s*Total charge                :\s*(?P<x_turbomole_total_charge>[-+0-9.eEdD]+)"),
                   SM (startReStr = r"\s*\|\s*Coordinates of all systems centered about cell 0\s*\|",
                       subMatchers = [
                           SM (r"\s*Redefined unit cell content (au):"),
                           SM (r"\s*Label               Cartesian Coordinates            Charge"),
                           SM (r"\s*(?P<x_turbomole_embed_geometry_atom_label>[a-zA-Z]+)\s+"
                               "(?P<x_turbomole_embed_geometry_atom_positions_x__angstrom>[-+0-9.]+)\s+"
                               "(?P<x_turbomole_embed_geometry_atom_positions_y__angstrom>[-+0-9.]+)\s+"
                               "(?P<x_turbomole_embed_geometry_atom_positions_z__angstrom>[-+0-9.]+)\s+"
                               "(?P<x_turbomole_embed_geometry_atom_charge>[-+0-9.]+)", repeats = True),
                           SM (r"\s*QM cluster transformed to the center of cell 0 \(au\)\:"),
                           SM (r"\s*Atom               Cartesian Coordinates"),
                           SM (r"\s*(?P<x_turbomole_embed_qm_cluster_geometry_atom_label>[a-zA-Z]+)\s+"
                               "(?P<x_turbomole_embed_qm_cluster_geometry_atom_positions_x__angstrom>[-+0-9.]+)\s+"
                               "(?P<x_turbomole_embed_qm_cluster_geometry_atom_positions_y__angstrom>[-+0-9.]+)\s+"
                               "(?P<x_turbomole_embed_qm_cluster_geometry_atom_positions_z__angstrom>[-+0-9.]+)\s+",repeats = True)
                       ])

               ])

def build_relaxation_matcher():
    return SM (name = "relaxation",
               #sections = ["section_single_configuration_calculation"],
               sections = ["section_single_configuration_calculation"],
               startReStr = r"\s*CONVERGENCY CRITERIA (?P<x_turbomole_geometry_optimization_converged>FULFILLED) IN CYCLE",
               subMatchers = [])

def get_cachingLevelForMetaName(metaInfoEnv):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.

    Returns:
        Dictionary with metaname as key and caching level as value. 
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
        'eigenvalues_eigenvalues': CachingLevel.Cache,
        'eigenvalues_kpoints':CachingLevel.Cache,
        'x_turbomole_geometry_optimization_converged': CachingLevel.Cache
    }

    # Set caching for temparary storage variables
    for name in metaInfoEnv.infoKinds:
        if (   name.startswith('x_turbomole_controlInOut')
               or name.startswith('x_turbomole_geometry')
               or name.startswith('x_turbomole_embed')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main():
    """Main function.

    Set up everything for the parsing of the turbomole main file and run the parsing.
    """
    # get main file description
    context = TurbomoleParserContext()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json
    path = "../../../../nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json"
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), path))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name': 'turbomole-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv)
    # start parsing
    mainFunction(mainFileDescription=build_root_parser(context),
                 metaInfoEnv=metaInfoEnv,
                 parserInfo=parserInfo,
                 cachingLevelForMetaName=cachingLevelForMetaName,
                 superContext=context)


if __name__ == "__main__":
    main()
