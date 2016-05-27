import setup_paths
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import get_metaInfo
import logging, os, re, sys

############################################################
# This is the parser for the main file of qbox.
############################################################

logger = logging.getLogger("nomad.turbomoleParser") 

class TurbomoleParserContext(object):

    def __init__(self):
        self.functionals                       = []

    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """

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

    ###################################################################
    # (3.4) onClose for geometry and force (section_system)
    # todo: maybe we can move the force to onClose_section_single_configuration_calculation in the future. 
    ###################################################################
    def onClose_section_system(self, backend, gIndex, section):
        """Trigger called when section_system is closed.
        Writes atomic positions, atom labels and lattice vectors.
        """
        # keep track of the latest system description section
        self.secSystemDescriptionIndex = gIndex

       #------1.atom_position
        atom_pos = []
        for i in ['x', 'y', 'z']:
            api = section['turbomole_geometry_atom_positions_' + i]
            if api is not None:
               atom_pos.append(api)
        if atom_pos:
            # need to transpose array since its shape is [number_of_atoms,3] in the metadata
           backend.addArrayValues('atom_positions', np.transpose(np.asarray(atom_pos)))

        #------2.atom labels
        atom_labels = section['turbomole_geometry_atom_labels']
        if atom_labels is not None:
           backend.addArrayValues('atom_labels', np.asarray(atom_labels))

                

#############################################################
#################[2] MAIN PARSER STARTS HERE  ###############
#############################################################

def build_TurbomoleMainFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the main file of qbox.

    First, several subMatchers are defined, which are then used to piece together
    the final SimpleMatcher.
    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses main file of Turbomole. 
    """

    ########################################                                    
    # submatcher for aims output from the parsed control.in                     
    controlInOutSubMatcher = SM (name = 'ControlInOut',                         
        startReStr = r"\s*\|\s*basis set information\s*\|\s",                          
        subMatchers = [                                                         
        SM (name = 'ControlInOutLines',                                         
            startReStr = r"\s*we will work with the",                                         
            sections = ['section_topology'],                                    
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
                sections = ['section_basis_set','section_system'],                                    
                subMatchers = [
                # SM (r"\s*-{20}-*", weak = True),                                                 
                SM (r"\s*(?P<turbomole_geometry_atom_labels>[a-zA-Z]+)\s*[0-9]\s*(?P<turbomole_basis_prim_number>[0-9]+)\s*(?P<turbomole_basis_cont_number>[0-9]+)\s*(?P<turbomole_basis_type>[a-zA-Z-a-zA-Z]+)"
                   ,repeats = True)
                ]),
            # only the first character is important for aims                    
            SM (r"\s*total number of primitive shells\s*:\s*(?P<turbomole_tot_primitive_shells>[0-9]+)",sections = ['section_basis_set'], repeats = True),
            SM (r"\s*total number of contracted shells\s*:\s*(?P<turbomole_tot_contracted_shells>[0-9]+)",sections = ['section_basis_set'], repeats = True),
            SM (r"\s*total number of cartesian basis functions\s*:\s*(?P<turbomole_tot_cartesian_func>[0-9]+)",sections = ['section_basis_set'], repeats = True),   
            SM (r"\s*total number of SCF-basis functions\s*:\s*(?P<turbomole_tot_scf_basis_func>[0-9]+)",sections = ['section_basis_set'], repeats = True),
            SM (name = 'Density functional informations',                                
                startReStr = r"\s*density functional",           
                sections = ['section_system'],  
                subMatchers = [                                                 
                SM (r"\s*(?P<turbomole_functional_type>[a-zA-Z-a-zA-Z0-9]+)\s*(?: functional)"),
                SM (r"\s*exchange:\s*(?P<turbomole_functional_type_exchange>[a-zA-Z-+a-zA-Z0-9\(\)\s]+)"),
                SM (r"\s*correlation:\s*(?P<turbomole_functional_type_correlation>[a-zA-Z-+a-zA-Z0-9\(\)\s]+)"),
                SM (r"\s*spherical integration\s*:\s*(?P<turbomole_grid_integration>[a-zA-Z\'\s]+)"),
                SM (r"\s*spherical gridsize\s*:\s*(?P<turbomole_grid_size>[0-9]+)"),
                SM (r"\s*i\.e\. gridpoints\s*:\s*(?P<turbomole_grid_points_number>[0-9]+)"),
                SM (r"\s*radial integration\s*:\s*(?P<turbomole_grid_radial_integration>[a-zA-Z0-9\(\)\s]+)"),
                SM (r"\s*radial gridsize\s*:\s*(?P<turbomole_grid_radial_gridsize>[0-9]+)"),
                SM (r"\s*integration cells\s*:\s*(?P<turbomole_grid_integration_cells>[0-9]+)"),
                SM (r"\s*partition function\s*:\s*(?P<turbomole_grid_partition_func>[a-zA-Z]+)"),
                SM (r"\s*partition sharpness\s*:\s*(?P<turbomole_grid_partition_sharpness>[0-9]+)")
                  ])  
                ]) # END ControlInOutLines                                     
        #SM (r"\s*-{20}-*", weak = True)                                         
        ])    

    #####################################################################
    # subMatcher for geometry                                                   
    # the verbatim writeout of the geometry.in is not considered for getting the structure data
    # using the geometry output of aims has the advantage that it has a clearer structure
    geometrySubMatcher = SM (name = 'Geometry',                                 
        startReStr = r"\s*\|\s*Atomic coordinate\, charge and isotop information\s\|\s",
        sections = ['section_system'],                              
        subMatchers = [                                                         
        SM (r"\s*-{20}-*", weak = True),                                        
        SM (startReStr = r"\s*atomic coordinates\s*atom    charge  isotop\s",   
            subMatchers = [                                                     
            SM (r"\s*(?P<turbomole_geometry_atom_positions_x__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_geometry_atom_positions_y__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_geometry_atom_positions_z__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_geometry_atom_labels>[a-zA-Z]+)", repeats = True)
            ])                                                                  
        ])                                                                      

    ########################################                                    
    # submatcher for total energy components during SCF interation              
    TotalEnergyScfSubMatcher = SM (name = 'TotalEnergyScf',                    
        repeats =True, 
        startReStr = r"\s*current damping\s*:\s*",                          
        forwardMatch = True,
        subMatchers = [                                                         
        SM (r"\s*current damping\s*:\s*(?P<turbomole_energy_scf_damping>[0-9.eEdD]+)"),
        SM (r"\s*(?P<turbomole_iteration_number>[0-9]+)\s*(?P<turbomole_energy_total_scf_iteration__eV>[-+0-9.eEdD]+)\s*(?P<turbomole_energy_one_scf_iteration__eV>[-+0-9.eEdD]+)"
             "\s*(?P<turbomole_energy_two_scf_iteration__eV>[-+0-9.eEdD]+)\s*(?P<turbomole_energy_norm_scf_iteration__eV>[-+0-9.eEdD]+)\s*(?P<turbomole_energy_tolerance_scf_iteration__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*max. resid. norm for Fia\-block\=\s*(?P<turbomole_max_res_norm_fia_block>[-+0-9.eEdD]+)\s*for orbital\s*(?P<turbomole_orbital_name_fia_block>[a-z0-9\s]+)"),
        SM (r"\s*max. resid. fock norm\s*\=\s*(?P<turbomole_max_res_norm_fock_norm>[-+0-9.eEdD]+)\s*for orbital\s*(?P<turbomole_orbital_name_fock_norm>[a-z0-9\s]+)"),
        SM (r"\s*irrep a   \: virtual orbitals shifted by\s*(?P<turbomole_virtual_orbital_shift>[0-9.]+)"),
        SM (r"\s*Delta Eig\.\s*\=\s*(?P<turbomole_delta_eigenvalues__eV>[-+0-9.eEdD]+)\s*eV")
        ])   
    ########################################                                    
    # submatcher for final total energy components               
    TotalEnergySubMatcher = SM (name = 'TotalEnergyFinal',                     
        startReStr = r"\s*\|\s*total energy\s*\=",                            
        forwardMatch = True,
        subMatchers = [                                                         
        SM (r"\s*\|\s*total energy\s*\=\s*(?P<turbomole_total_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*kinetic energy\s*\=\s*(?P<turbomole_kinetic_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*potential energy\s*\=\s*(?P<turbomole_potential_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*virial theorem\s*\=\s*(?P<turbomole_virial_theorem_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*wavefunction norm\s*\=\s*(?P<turbomole_wave_func_norm__eV>[-+0-9.eEdD]+)")
        
        ]) 
    ########################################                                    
    # return main Parser                                                        
    ########################################                                    
    return SM (name = 'Root',                                                   
                                                                                
        startReStr = "",                                                        
        forwardMatch = True,                                                    
        weak = True,
        sections = ['section_run'],                                                            
        subMatchers = [                                                         
            SM (name = 'ProgramHeader',                                         
                startReStr = r"\s*RUNNING PROGRAM",                    
                subMatchers = [                                                 
                SM (r"\s*dscf \(thnec244\) \: TURBOMOLE V(?P<turbomole_program_version>[0-9.()\s]+)")       
                ]), # END ProgramHeader
        #=============================================================================
        #  read OUPUT file *.r, the method part comes from INPUT file *.i,  so we 
        #  do not need to parser INPUT file, the OUTPUT file contains all information
        #=============================================================================
        SM (name = 'NewRun',                                                    
            startReStr = r"\s*SCF run will be profiled \!",                        
            endReStr = r"\s*\*\*\*\*  dscf \: all done  \*\*\*\*",                                 
            repeats = True,                                                     
            required = True,                                                    
            forwardMatch = True,                                                
            fixedStartValues={'program_name': 'Turbomole', 'program_basis_set_type': 'numeric AOs' },
            subMatchers = [                                                     
	    #controlInOutSubMatcher,
            SM (name = 'SectionMethod',                                         
                startReStr = r"\s*SCF run will be profiled \!",
                sections = ['section_method'],                                  
                subMatchers = [                                                 
                # parse geometry writeout of aims                               
                geometrySubMatcher,
                controlInOutSubMatcher                                              
                ]),

              # the actual section for a single configuration calculation starts here
            SM (name = 'SingleConfigurationCalculation',                    
                  startReStr = r"\s*start vectors will be provided from a core hamilton",
                  repeats = True,                                             
                  subMatchers = [
                  SM (name = 'TotalEnergyForEachScfCycle',                            
                      startReStr = r"\s*STARTING INTEGRAL EVALUATION FOR 1st SCF ITERATION",
                      sections = ['section_scf_iteration'],                   
                      subMatchers = [                                         
                      TotalEnergyScfSubMatcher,
                      TotalEnergySubMatcher                                          
                      ])#, # END ScfInitialization  
                   ])#, # END SingleConfigurationCalculation
           ]), # CLOSING SM NewRun                                               
        ]) # END Root  

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
                                'eigenvalues_kpoints':CachingLevel.Cache
                                }

    # Set caching for temparary storage variables
    for name in metaInfoEnv.infoKinds:
        if (   name.startswith('turbomole_')
            or name.startswith('turbomole_')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName




def main():
    """Main function.

    Set up everything for the parsing of the qbox main file and run the parsing.
    """
    # get main file description
    TurbomoleMainFileSimpleMatcher = build_TurbomoleMainFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name':'turbomole-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv)
    # start parsing
    mainFunction(mainFileDescription = TurbomoleMainFileSimpleMatcher,
                 metaInfoEnv = metaInfoEnv,
                 parserInfo = parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = TurbomoleParserContext())

if __name__ == "__main__":
    main()

