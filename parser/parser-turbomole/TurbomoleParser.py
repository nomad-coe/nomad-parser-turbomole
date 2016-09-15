

from builtins import object
import setup_paths
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import get_metaInfo
import logging, os, re, sys

############################################################
# This is the parser for the main file of turbomole.
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
        self.lastCalculationGIndex = None

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

    def onClose_turbomole_section_functionals(self, backend, gIndex, section):
        functional_names = section["XC_functional_type"]

        s = backend.openSection("section_XC_functionals")
        backend.addValue('XC_functional_name', functional_names[-1])
        backend.closeSection("section_XC_functionals", s)

    def onClose_section_system(self, backend, gIndex, section):
        """Trigger called when section_system is closed.
        Writes atomic positions, atom labels and lattice vectors.
        """
        # keep track of the latest system description section
        self.secSystemDescriptionIndex = gIndex

        method_name = section['electronic_structure_method']
        if method_name is None:
                match = 'DFT'
                backend.addValue('electronic_structure_method', match)

        smear_type = section['smearing_kind']
        if smear_type is None:
                value = ''
                backend.addValue('smearing_kind', value)

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

def build_TurbomoleMainFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the main file of turbomole.

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
        startReStr = r"\s*\|\s*basis set information\s",                      
	#startReStr = r"\s*SCF run will be profiled",

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
                sections = ['section_basis_set'],
                subMatchers = [
                # SM (r"\s*-{20}-*", weak = True),                                                 
                SM (r"\s*(?P<turbomole_controlInOut_atom_labels>[a-zA-Z]+)\s*[0-9]+\s*(?P<turbomole_controlInOut_basis_prim_number>[0-9]+)\s*(?P<turbomole_controlInOut_basis_cont_number>[0-9]+)\s*(?P<turbomole_controlInOut_basis_type>[a-zA-Z-a-zA-Z]+)"
                   ,repeats = True)
                ]),
            # only the first character is important for aims                    
            SM (r"\s*total number of primitive shells\s*:\s*(?P<turbomole_controlInOut_tot_primitive_shells>[0-9]+)",sections = ['section_basis_set'], repeats = True),
            SM (r"\s*total number of contracted shells\s*:\s*(?P<turbomole_controlInOut_tot_contracted_shells>[0-9]+)",sections = ['section_basis_set'], repeats = True),
            SM (r"\s*total number of cartesian basis functions\s*:\s*(?P<turbomole_controlInOut_tot_cartesian_func>[0-9]+)",sections = ['section_basis_set'], repeats = True),   
            SM (r"\s*total number of SCF-basis functions\s*:\s*(?P<turbomole_controlInOut_tot_scf_basis_func>[0-9]+)",sections = ['section_basis_set'], repeats = True),
	]), # END ControlInOutLines
        SM (name = 'post-HF',
            startReStr = r"\s*(?:[a-zA-Z-a-zA-Z0-9\s]+)\s*shell calculation for the wavefunction models",
            sections = ['section_system'],
            subMatchers = [
                SM (r"\s*(?P<electronic_structure_method>[a-zA-Z-a-zA-Z0-9\(\)]+)\s*\-")
            ]),
#	SM (name = "smearing",
#	    startReStr = r"\s*Fermi\s*smearing\s*switched\s*on",
#	    sections = ["section_system"],
#	    subMatchers = [
#		SM (r"\s*(?P<smearing_kind>[a-zA-Z]+)\s*smearing switched on"),
#		SM (r"\s*Final electron temperature\:\s*(?P<smearing_width>[0-9.eEdD]+)")
#	    ]),
	SM (name = 'XC functional',
	    startReStr = r"\s*density functional",
	    sections = ['turbomole_section_functionals'],
	    subMatchers = [ 
	    SM (r"\s*(?P<XC_functional_type>[a-zA-Z-a-zA-Z0-9]+)\s*(?: functional)"),
	    SM (r"\s*(?P<XC_functional_type>[a-zA-Z-a-zA-Z0-9]+)\s*(?: meta-GGA functional)"),
	    SM (r"(?:[a-zA-Z-a-zA-Z0-9\s]+)\s*functional\:\s*(?P<XC_functional_type>[a-zA-Z-a-zA-Z0-9]+)"),
            SM (r"\s*exchange:\s*(?P<turbomole_controlInOut_functional_type_exchange>[a-zA-Z-+a-zA-Z0-9\(\)\s.\*]+)"),
            SM (r"\s*correlation:\s*(?P<turbomole_controlInOut_functional_type_correlation>[a-zA-Z-+a-zA-Z0-9\(\)\s.\*]+)"),
            SM (r"\s*spherical integration\s*:\s*(?P<turbomole_controlInOut_grid_integration>[a-zA-Z\'\s]+)"),
            SM (r"\s*spherical gridsize\s*:\s*(?P<turbomole_controlInOut_grid_size>[0-9]+)"),
            SM (r"\s*i\.e\. gridpoints\s*:\s*(?P<turbomole_controlInOut_grid_points_number>[0-9]+)"),
            SM (r"\s*radial integration\s*:\s*(?P<turbomole_controlInOut_grid_radial_integration>[a-zA-Z0-9\(\)\s]+)"),
            SM (r"\s*radial gridsize\s*:\s*(?P<turbomole_controlInOut_grid_radial_grid_size>[0-9]+)"),
            SM (r"\s*integration cells\s*:\s*(?P<turbomole_controlInOut_grid_integration_cells>[0-9]+)"),
            SM (r"\s*partition function\s*:\s*(?P<turbomole_controlInOut_grid_partition_func>[a-zA-Z]+)"),
            SM (r"\s*partition sharpness\s*:\s*(?P<turbomole_controlInOut_grid_partition_sharpness>[0-9]+)")
              ])
       #SM (r"\s*-{20}-*", weak = True)                                         
        ])    

    ########################################                                    
    # submatcher for eigenvalues                                                
    # first define function to build subMatcher for normal case and scalar ZORA 

    # Here the separation between the parsed string of eigenvalues needs to be done! 
    def build_eigenvaluesGroupSubMatcher(addStr):                               
        """Builds the SimpleMatcher to parse the normal and the scalar ZORA eigenvalues in aims.
                                                                                
        Args:                                                                   
            addStr: String that is appended to the metadata names.              
                                                                                
        Returns:                                                                
            SimpleMatcher that parses eigenvalues with metadata according to addStr. 
        """                                                                     
        # submatcher for eigenvalue list                                        
        EigenvaluesListSubMatcher =  SM (name = 'EigenvaluesLists',            
	    repeats =True, 
#	    startReStr = r"\s*orbitals $uhfmo_alpha  will be written to file alpha\s*",
            startReStr = r"\s*(?: alpha|beta)\:\s*",
            sections = ['turbomole_section_eigenvalues_list%s' % addStr],        
            subMatchers = [                                                     
		 SM (r"\s*(?: irrep)\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)" % (1 * (addStr,)), repeats = True),
		 #SM (r"\s*(?: irrep)\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)" % (addStr, addStr, addStr, addStr, addStr), repeats = True), 	+	#SM (r"\s*(?: irrep)\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)" % addStr, r"\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)" % addStr, r"\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)" % addStr, r"\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)" % addStr, r"\s*(?P<turbomole_irreducible_representation_state%s>[0-9a\s]+)" % addStr, repeats = True),
		 SM (r"\s*(?: eigenvalues H)\s*[-+0-9.eEdD\s]+", repeats = True),
                 SM (r"\s*(?: eV)\s*(?P<turbomole_eigenvalue_eigenvalue%s>[-+0-9.eEdD]+)" % (1 * (addStr,)), repeats = True)#,
#                 SM (r"\s*(?: occupation)\s*(?P<turbomole_eigenvalue_occupation%s>[0-9.\s]+)" % (1 * (addStr,)), repeats = True)
            ]) 
        return SM (name = 'EigenvaluesGroup',                                   
	    startReStr = "\s*orbitals\s*\$",
            #startReStr = "\s*(?: alpha|beta)\:\s*",
            sections = ['turbomole_section_eigenvalues_group%s' % addStr],       
            subMatchers = [                                                     
            SM (name = 'EigenvaluesNoSpin',                          
                startReStr = r"\s*(?: alpha|beta)\:\s*",
                sections = ['turbomole_section_eigenvalues_spin%s' % addStr],    
                forwardMatch = True,                                            
                subMatchers = [                                                 
                EigenvaluesListSubMatcher.copy()                                
                ]), # END EigenvaluesNoSpinNonPeriodic                          
            ])                                                                  
   # now construct the two subMatchers                                         
    EigenvaluesGroupSubMatcher = build_eigenvaluesGroupSubMatcher('')   
    #####################################################################
    # subMatcher for geometry                                                   
    # the verbatim writeout of the geometry.in is not considered for getting the structure data
    # using the geometry output of aims has the advantage that it has a clearer structure
    geometrySubMatcher = SM (name = 'Geometry',                                 
	startReStr = r"\s*\|\s*Atomic coordinate",
        sections = ['section_system'],                              
        subMatchers = [                                                         
        SM (r"\s*-{20}-*", weak = True),                                        
	SM (startReStr = r"\s*atomic coordinates",
            subMatchers = [                                                     
            SM (r"\s*(?P<turbomole_geometry_atom_positions_x__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_geometry_atom_positions_y__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_geometry_atom_positions_z__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_geometry_atom_labels>[a-zA-Z]+)\s+(?P<turbomole_geometry_atom_charge>[0-9.]+)", repeats = True)
            ])                                                                  
        ])                                                                      
    ########################################
    # submatcher for atomic forces
    ForcesMatcher = SM (name = 'AtomicForces',
	repeats =True, 
	startReStr = r"\s*ATOM\s*CARTESIAN GRADIENTS",
	#forwardMatch = True,
	subMatchers = [
	    SM (r"\s*(?:[0-9]+)\s*(?:[a-z]+)\s*(?P<atom_forces_raw_x__hartree_bohr_1>[-+0-9.eEdD]+)\s*(?P<atom_forces_raw_y__hartree_bohr_1>[-+0-9.eEdD]+)\s*(?P<atom_forces_raw_z__hartree_bohr_1>[-+0-9.eEdD]+)", repeats = True)
	])
    ########################################                                    
    # submatcher for total energy components during SCF interation              
    TotalEnergyScfSubMatcher = SM (name = 'TotalEnergyScf',                    
        repeats =True, 
	#startReStr = r"\s*scf convergence criterion",
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
        SM (r"\s*\|\s*total energy\s*\=\s*(?P<energy_total__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*kinetic energy\s*\=\s*(?P<turbomole_kinetic_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*potential energy\s*\=\s*(?P<turbomole_potential_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*virial theorem\s*\=\s*(?P<turbomole_virial_theorem_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\:\s*wavefunction norm\s*\=\s*(?P<turbomole_wave_func_norm__eV>[-+0-9.eEdD]+)")
        
        ]) 
    SmearingOccupation = SM (name = "smearing",
	startReStr = r"\s*and increment of one",
        sections = ["section_system"],
        subMatchers = [
            SM (r"\s*(?P<smearing_kind>[a-zA-Z]+)\s*smearing switched on"),
            SM (r"\s*Final electron temperature\:\s*(?P<smearing_width>[0-9.eEdD]+)")
        ])
    ########################################
    # submatcher for coupled-cluster and MP2 energy
    CCEnergySubMatcher = SM (name = 'TotalEnergyCC',
        #startReStr = r"\s*Calculate\s*integrals\s*\(*ia\|*jb\)*\s*for MP2 start guess",
	startReStr = r"\s*=========================================================================",
        forwardMatch = True,
        subMatchers = [
        SM (r"\s*\*\s*RHF  energy\s*\:\s*(?P<turbomole_HF_total_energy_final__eV>[-+0-9.eEdD]+)"),
	SM (r"\s*\*\s*UHF  energy\s*\:\s*(?P<turbomole_HF_total_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\*\s*Final MP2 energy\s*\:\s*(?P<turbomole_MP2_total_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\*\s*Final CCSD energy\s*\:\s*(?P<turbomole_CCSD_total_energy_final__eV>[-+0-9.eEdD]+)"),
	SM (r"\s*\*\s*Final CC2 energy\s*\:\s*(?P<turbomole_CC2_total_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\*\s*Final CCSD\(T\) energy\s*\:\s*(?P<turbomole_CCSDparT_total_energy_final__eV>[-+0-9.eEdD]+)"),
        SM (r"\s*\*\s*D1 diagnostic \(CCSD\)\s*\:\s*(?P<turbomole_D1_diagnostic>[-+0-9.eEdD]+)")
        
        ])
    ########################################
    # submatcher for perturbation theory total energy
    PTEnergySubMatcher = SM (name = 'TotalEnergyPT',
	startReStr = r"\s*\|*\s*| natural orb",
	forwardMatch = True,
	subMatchers = [
	SM (r"\s*Total Energy\s*\:\s*(?P<turbomole_PT_total_energy_final__eV>[-+0-9.eEdD]+)")

	])
    ########################################                                    
    # submatcher for final total energy components                              
    EmbeddingSubMatcher = SM (name = 'PeriodicEmbedding',                      
        startReStr = r"\s*\+------------------------ Parameters ------------------------\+",                              
        forwardMatch = True,                                                    
        subMatchers = [                                                         
        SM (r"\s*Maximum multipole moment used               :\s*(?P<turbomole_max_multipole_moment>[0-9]+)"),
        SM (r"\s*Multipole precision parameter               :\s*(?P<turbomole_multipole_precision_parameter>[-+0-9.eEdD]+)"),
        SM (r"\s*Minimum separation between cells            :\s*(?P<turbomole_min_separation_cells>[-+0-9.eEdD]+)"),
        SM (r"\s*\+-----------------------------------------------------------\+\s*"),
        SM (r"\s*Charge Neutrality tolerance :\s*(?P<turbomole_charge_neutrality_tol>[-+0-9.eEdD]+)"),
        SM (r"\s*Total charge                :\s*(?P<turbomole_total_charge>[-+0-9.eEdD]+)"),
        SM (startReStr = r"\s*\|\s*Coordinates of all systems centered about cell 0\s*\|",   
            subMatchers = [
            SM (r"\s*Redefined unit cell content (au):"),                                                     
            SM (r"\s*Label               Cartesian Coordinates            Charge"),                                                     
            SM (r"\s*(?P<turbomole_embed_geometry_atom_label>[a-zA-Z]+)\s+"
                 "(?P<turbomole_embed_geometry_atom_positions_x__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_embed_geometry_atom_positions_y__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_embed_geometry_atom_positions_z__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_embed_geometry_atom_charge>[-+0-9.]+)", repeats = True),
            SM (r"\s*QM cluster transformed to the center of cell 0 \(au\)\:"),                       
            SM (r"\s*Atom               Cartesian Coordinates"),
            SM (r"\s*(?P<turbomole_embed_qm_cluster_geometry_atom_label>[a-zA-Z]+)\s+"     
                 "(?P<turbomole_embed_qm_cluster_geometry_atom_positions_x__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_embed_qm_cluster_geometry_atom_positions_y__angstrom>[-+0-9.]+)\s+"
                 "(?P<turbomole_embed_qm_cluster_geometry_atom_positions_z__angstrom>[-+0-9.]+)\s+",repeats = True)
            ]) 
                                                                                
        ])  
    ########################################
    # submatcher for pertubative GW eigenvalues
    # first define function to build subMatcher
    def build_GWeigenvaluesGroupSubMatcher(addStr):
       """Builds the SimpleMatcher to parse the perturbative GW eigenvalues in turbomole.

       Args:
           addStr: String that is appended to the metadata names.

       Returns:
           SimpleMatcher that parses eigenvalues with metadata according to addStr.
       """
       # submatcher for eigenvalue list
       GWEigenvaluesListSubMatcher = SM (name = 'perturbativeGW_EigenvaluesLists',
#	   startReStr = r"\s*in\s*eV",
	   startReStr = r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
           sections = ['turbomole_section_eigenvalues_list%s' % addStr],
           subMatchers = [
	   SM (r"\s*in\s*eV"),
	   SM (r"\s*------------------------------------------------------------------------------------"),
	   SM (r"\s*(?P<eigenstate_number>[0-9]+)\s+(?P<turbomole_eigenvalue_ks_GroundState__eV>[-+0-9.eEdD]+)\s+"
			     "(?P<turbomole_eigenvalue_quasiParticle_energy__eV>[-+0-9.eEdD]+)\s+"
			     "(?P<turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
			     "(?P<turbomole_eigenvalue_ExactExchange_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
			     "(?P<turbomole_eigenvalue_correlation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
			     "(?P<turbomole_eigenvalue_ks_ExchangeCorrelation__eV>[-+0-9.eEdD]+)\s+"
			     "(?P<turbomole_Z_factor>[-+0-9.eEdD]+)\s+"
			     "(?P<turbomole_ExchangeCorrelation_perturbativeGW_derivation>[-+0-9.eEdD]+)", 
           adHoc = lambda parser: parser.superContext.setStartingPointCalculation(parser),
           repeats = True),
	   SM (r"\s*------------------------------------------------------------------------------------"),
           SM (r"\s*(?P<eigenstate_number>[0-9]+)\s+(?P<turbomole_eigenvalue_ks_GroundState__eV>[-+0-9.eEdD]+)\s+"
                             "(?P<turbomole_eigenvalue_quasiParticle_energy__eV>[-+0-9.eEdD]+)\s+"
                             "(?P<turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                             "(?P<turbomole_eigenvalue_ExactExchange_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                             "(?P<turbomole_eigenvalue_correlation_perturbativeGW__eV>[-+0-9.eEdD]+)\s+"
                             "(?P<turbomole_eigenvalue_ks_ExchangeCorrelation__eV>[-+0-9.eEdD]+)\s+"
                             "(?P<turbomole_Z_factor>[-+0-9.eEdD]+)\s+"
                             "(?P<turbomole_ExchangeCorrelation_perturbativeGW_derivation>[-+0-9.eEdD]+)", 
	   adHoc = lambda parser: parser.superContext.setStartingPointCalculation(parser),
	   repeats = True)
           ])
       return SM (name = 'perturbativeGW_EigenvaluesGroup',
           startReStr = r"\s*GW\s*version:",
           sections = ['turbomole_section_eigenvalues_group%s' % addStr],
           subMatchers = [
           # non-spin-polarized
           SM (name = 'GW_EigenvaluesNoSpinNonPeriodic',
               startReStr = r"\s*orb\s+eps\s+QP-eps\s+Sigma\s+Sigma_x\s+Sigma_c\s+Vxc\s+Z\s+dS\/de",
               sections = ['turbomole_section_eigenvalues_spin%s' % addStr],
               forwardMatch = True,
               subMatchers = [
#               SM (r"\s*-+"),
#               SM (r"\s*-+"),
               GWEigenvaluesListSubMatcher.copy()
               ]), # END EigenvaluesNoSpinNonPeriodic
           ])
    # now construct the two subMatchers
    GWEigenvaluesGroupSubMatcher = build_GWeigenvaluesGroupSubMatcher('_perturbativeGW')

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
		startReStr = r"",
                subMatchers = [                                                 
		SM (r"\s*(?:aoforce|cosmoprep|egrad|evib|frog|gradsammel|hesssammel|moloch|odft|relax|ridft|rirpa|sdg|thirdsammel|vibration|atbandbta|define|eigerf|fdetools|grad|haga|intense|mpgrad|proper|ricc2|rimp2|ruecker|statpt|tm2molden|woelfling|bsseenergy|dscf|escf|freeh|gradruecker|hessruecker|mdprep|mpshift|rdgrad|ricctools|rimp2prep|sammler|thirdruecker|uff)\s*\((?P<turbomole_nodename>[a-zA-Z0-9.]+)\) \: TURBOMOLE (?P<program_version>[a-zA-Z0-9.]+)")
                ]), # END ProgramHeader
        #=============================================================================
        #  read OUPUT file *.r, the method part comes from INPUT file *.i,  so we 
        #  do not need to parser INPUT file, the OUTPUT file contains all information
        #=============================================================================
        SM (name = 'NewRun',  
            startReStr = r"\s*Copyright \(C\) ",
	    endReStr = r"\s*\*\*\*\*\s",
            repeats = True,                                                     
            required = True,                                                    
            forwardMatch = True,                                                
            fixedStartValues={'program_name': 'Turbomole', 'program_basis_set_type': 'GTOs' },
            subMatchers = [                                                     
	    #controlInOutSubMatcher,
            SM (name = 'SectionMethod',                                         
                startReStr = r"\s*Copyright \(C\) ",
                sections = ['section_method'],                                  
                subMatchers = [                                                 
                # parse geometry writeout of aims                               
		geometrySubMatcher,
                controlInOutSubMatcher
		]),

              # the actual section for a single configuration calculation starts here
            SM (name = 'SingleConfigurationCalculation',                    
                  #startReStr = r"\s*start vectors will be provided from a core hamilton",
		  startReStr = r"\s*1e\-*integrals will be neglected if expon",
		  #startReStr = r"\s*\|",
                  repeats = True,                                             
                  subMatchers = [
                  SM (name = 'PeriodicEmbeddingSettings',                      
                      startReStr = r"\s*\|\s*EMBEDDING IN PERIODIC POINT CHARGES\s*\|",
                      sections = ['section_method'],                     
                      subMatchers = [                                           
                      EmbeddingSubMatcher                                     
                      ]),
                  SM (name = 'TotalEnergyForEachScfCycle',                            
		      startReStr = r"\s*scf convergence criterion",
                      #startReStr = r"\s*STARTING INTEGRAL EVALUATION FOR 1st SCF ITERATION",
                      sections = ['section_scf_iteration'],                   
                      subMatchers = [                                         
                      #EigenvaluesGroupSubMatcher,  
                      SmearingOccupation,
		      TotalEnergyScfSubMatcher,
                      TotalEnergySubMatcher#,
		      #SmearingOccupation
                      ]), # END ScfInitialization 
                  SM (name = 'EigenvaluesGroupSubMatcher',                      
                      #startReStr = r"\s+orbitals [a-zA-Z\$\_]+ (?: will be written to file) [a-zA-Z]+",
		      startReStr = r"\s*orbitals\s*\$",
                      #sections = ['section_scf_iteration'],                     
                      subMatchers = [                                           
                      EigenvaluesGroupSubMatcher                               
                      ])
                  ]), # END SingleConfigurationCalculation
	    ForcesMatcher,
            SM (name = 'PostHFTotalEnergies',
                startReStr = r"\s*Energy of reference wave function is",
                sections = ['section_scf_iteration'],
                subMatchers = [
                CCEnergySubMatcher
                ]),
	    SM (name = 'PTTotalEnergies',
		startReStr = r"\s*\|\s*MP2 relaxed",
		sections = ['section_scf_iteration'],
		subMatchers = [
		PTEnergySubMatcher
		]),
	    GWEigenvaluesGroupSubMatcher
           ]) # CLOSING SM NewRun                                               
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
        if (   name.startswith('turbomole_controlInOut')
            or name.startswith('turbomole_geometry')
            or name.startswith('turbomole_embed')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName




def main():
    """Main function.

    Set up everything for the parsing of the turbomole main file and run the parsing.
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

