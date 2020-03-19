import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import public

m_package = Package(
    name='turbomole_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='turbomole.nomadmetainfo.json'))


class x_turbomole_section_eigenvalues_GW(MSection):
    '''
    section for the eigenvalues of a GW calculation (at present only pertubative G0W0)
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_turbomole_section_eigenvalues_GW'))

    x_turbomole_eigenvalue_correlation_perturbativeGW = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Correlation energy at a given eigenstate from perturbative GW
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_eigenvalue_correlation_perturbativeGW'))

    x_turbomole_eigenvalue_ExactExchange_perturbativeGW = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Exact exchange energy at given eigenstate from perturbative GW
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_eigenvalue_ExactExchange_perturbativeGW'))

    x_turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Self-energy at a given eigenstate from perturbative GW
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_eigenvalue_ExchangeCorrelation_perturbativeGW'))

    x_turbomole_eigenvalue_ks_ExchangeCorrelation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        KS exchange correlation energy at a given eigenstate needed to calculate the
        quasi-particle energy in perturbative GW
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_eigenvalue_ks_ExchangeCorrelation'))

    x_turbomole_eigenvalue_ks_GroundState = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        KS ground state energy at a given eigenstate needed in perturbative GW
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_eigenvalue_ks_GroundState'))

    x_turbomole_eigenvalue_quasiParticle_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Quasiparticle energy at a given eigenstate from perturbative GW
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_eigenvalue_quasiParticle_energy'))

    x_turbomole_ExchangeCorrelation_perturbativeGW_derivation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        TODO:
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_ExchangeCorrelation_perturbativeGW_derivation'))

    x_turbomole_Z_factor = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        TODO:
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_Z_factor'))


class x_turbomole_section_functionals(MSection):
    '''
    section for one list of XC functionals
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_turbomole_section_functionals'))

    x_turbomole_controlInOut_grid_integration_cells = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The integration cells
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_integration_cells'))

    x_turbomole_controlInOut_grid_integration = Quantity(
        type=str,
        shape=[],
        description='''
        type of the used grid integration
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_integration'))

    x_turbomole_controlInOut_grid_partition_func = Quantity(
        type=str,
        shape=[],
        description='''
        Type of the partition function used
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_partition_func'))

    x_turbomole_controlInOut_grid_partition_sharpness = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Sharpness of the partition function
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_partition_sharpness'))

    x_turbomole_controlInOut_grid_points_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Grid points number
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_points_number'))

    x_turbomole_controlInOut_grid_radial_grid_size = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The size of the radial grid
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_radial_grid_size'))

    x_turbomole_controlInOut_grid_radial_integration = Quantity(
        type=str,
        shape=[],
        description='''
        The radial integration type
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_radial_integration'))

    x_turbomole_controlInOut_grid_size = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        The size of the used grid
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlInOut_grid_size'))

    x_turbomole_XC_functional_type = Quantity(
        type=str,
        shape=[],
        description='''
        XC functional type
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_XC_functional_type'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_turbomole_controlIn_atom_label = Quantity(
        type=str,
        shape=[],
        description='''
        The label of the atoms in the system
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_atom_label'))

    x_turbomole_controlIn_atom_number = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        The number of atoms in the system
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_atom_number'))

    x_turbomole_controlIn_basis_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status mean here ON or OFF
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_basis_status'))

    x_turbomole_controlIn_cartesian_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status mean here ON or OFF
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_cartesian_status'))

    x_turbomole_controlIn_damping_parameter_min = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_damping_parameter_min'))

    x_turbomole_controlIn_damping_parameter_start = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_damping_parameter_start'))

    x_turbomole_controlIn_damping_parameter_step = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_damping_parameter_step'))

    x_turbomole_controlIn_dipole_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status mean here ON or OFF
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_dipole_status'))

    x_turbomole_controlIn_global_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status mean here ON or OFF
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_global_status'))

    x_turbomole_controlIn_hessian_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status mean here ON or OFF
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_hessian_status'))

    x_turbomole_controlIn_interconversion_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status mean here ON or OFF
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_interconversion_status'))

    x_turbomole_controlIn_number_of_integral_stored = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_number_of_integral_stored'))

    x_turbomole_controlIn_operating_system = Quantity(
        type=str,
        shape=[],
        description='''
        The kind of operating system
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_operating_system'))

    x_turbomole_controlIn_pople_kind = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_pople_kind'))

    x_turbomole_controlIn_scf_conv = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_scf_conv'))

    x_turbomole_controlIn_scf_iter_limit = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_scf_iter_limit'))

    x_turbomole_controlIn_scfintunit_file = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_scfintunit_file'))

    x_turbomole_controlIn_scfintunit_size = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_scfintunit_size'))

    x_turbomole_controlIn_scfintunit_unit = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_scfintunit_unit'))

    x_turbomole_controlIn_symmetry = Quantity(
        type=str,
        shape=[],
        description='''
        The given symmetry of the system
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_symmetry'))

    x_turbomole_controlIn_time_for_integral_calc = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_controlIn_time_for_integral_calc'))

    x_turbomole_dft_d3_version = Quantity(
        type=str,
        shape=[],
        description='''
        version of the DFT-D3 van-der-Waals correction that is used
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_dft_d3_version'))

    x_turbomole_functional_type_correlation = Quantity(
        type=str,
        shape=[],
        description='''
        type of the used correlation functional
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_functional_type_correlation'))

    x_turbomole_functional_type_exchange = Quantity(
        type=str,
        shape=[],
        description='''
        type of the used exchange functional
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_functional_type_exchange'))

    x_turbomole_geometry_optimization_cycle_index = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        By default Turbomole only keeps the output of the final iteration once the
        geometry has been converged, thus the entire optimization trajectory cannot be
        rebuild in most cases. Instead, this value contains the optimization cycle index
        to indicate how many iterations have preceded this one.
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_geometry_optimization_cycle_index'))

    x_turbomole_gw_approximation = Quantity(
        type=str,
        shape=[],
        description='''
        The employed GW approximation.
        ''',
        categories=[public.settings_XC, public.settings_post_hartree_fock, public.settings_GW, public.settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='x_turbomole_gw_approximation'))

    x_turbomole_gw_eta_factor = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        [TO BE VERIFIED]Infinitesimal complex energy shift. Negative value switches to
        calculating at that value but extrapolating to 0 in linear approximation.
        ''',
        categories=[public.settings_XC, public.settings_post_hartree_fock, public.settings_GW, public.settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='x_turbomole_gw_eta_factor'))

    x_turbomole_gw_use_rpa_response = Quantity(
        type=bool,
        shape=[],
        description='''
        If true, the pure RPA response function is calculated. Otherwise, the TDDFT
        response function is calculated and used to screen the coulomb interaction.
        ''',
        categories=[public.settings_XC, public.settings_post_hartree_fock, public.settings_GW, public.settings_potential_energy_surface],
        a_legacy=LegacyDefinition(name='x_turbomole_gw_use_rpa_response'))

    x_turbomole_uhfmo_type = Quantity(
        type=str,
        shape=[],
        description='''
        Type of UHF molecular orbital
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_uhfmo_type'))

    x_turbomole_section_functionals = SubSection(
        sub_section=SectionProxy('x_turbomole_section_functionals'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_turbomole_section_functionals'))


class section_scf_iteration(public.section_scf_iteration):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_scf_iteration'))

    x_turbomole_damping_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Damping of the two-electron contributions to Fock matrix in the present SCF
        iteration
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_damping_scf_iteration'))

    x_turbomole_delta_eigenvalues = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        change of the eigenvalues in the current SCF iteration
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_turbomole_delta_eigenvalues'))

    x_turbomole_energy_1electron_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Total energy contribution from one-electron integrals
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_turbomole_energy_1electron_scf_iteration'))

    x_turbomole_energy_2electron_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Total energy contribution from two-electron integrals
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_turbomole_energy_2electron_scf_iteration'))

    x_turbomole_norm_diis_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Norm of the DIIS error in an SCF-iteration
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_norm_diis_scf_iteration'))

    x_turbomole_norm_fia_orbital_scf_iteration = Quantity(
        type=str,
        shape=[],
        description='''
        orbital with the largest residual norm for the Fia block in this iteration
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_norm_fia_orbital_scf_iteration'))

    x_turbomole_norm_fia_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Maximal resid. norm for Fia-block in an SCF-iteration
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_norm_fia_scf_iteration'))

    x_turbomole_norm_fock_orbital_scf_iteration = Quantity(
        type=str,
        shape=[],
        description='''
        orbital with the largest residual Fock norm in this iteration
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_norm_fock_orbital_scf_iteration'))

    x_turbomole_norm_fock_scf_iteration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Maximal resid. fock norm in an SCF-iteration
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_norm_fock_scf_iteration'))


class section_eigenvalues(public.section_eigenvalues):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_eigenvalues'))

    x_turbomole_eigenvalues_irreducible_representation = Quantity(
        type=np.dtype('U'),
        shape=['number_of_spin_channels', 'number_of_eigenvalues_kpoints', 'number_of_eigenvalues'],
        description='''
        Irreducible representation the eigenstates belong to.
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_eigenvalues_irreducible_representation'))

    x_turbomole_section_eigenvalues_GW = SubSection(
        sub_section=SectionProxy('x_turbomole_section_eigenvalues_GW'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_turbomole_section_eigenvalues_GW'))


class section_sampling_method(public.section_sampling_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_sampling_method'))

    x_turbomole_geometry_optimization_geometry_change_rms = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        geometry optimization convergence criterion - Root Mean Square of displacements
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_geometry_optimization_geometry_change_rms'))

    x_turbomole_geometry_optimization_threshold_force_rms = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='newton',
        description='''
        geometry optimization convergence criterion - Root Mean Square of forces
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_geometry_optimization_threshold_force_rms'))

    x_turbomole_geometry_optimization_trustregion_initial = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        geometry optimization trust region - initial radius
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_geometry_optimization_trustregion_initial'))

    x_turbomole_geometry_optimization_trustregion_max = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        geometry optimization trust region - maximum radius
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_geometry_optimization_trustregion_max'))

    x_turbomole_geometry_optimization_trustregion_min = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        geometry optimization trust region - minimum radius
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_geometry_optimization_trustregion_min'))


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    x_turbomole_module = Quantity(
        type=str,
        shape=[],
        description='''
        The name of the Turbomole module used for this single configuration calculation.
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_module'))

    x_turbomole_potential_energy_final = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        Final potential energy
        ''',
        categories=[public.energy_value, public.energy_component],
        a_legacy=LegacyDefinition(name='x_turbomole_potential_energy_final'))

    x_turbomole_vibrations_infrared_activity = Quantity(
        type=bool,
        shape=['x_turbomole_vibrations_num_modes'],
        description='''
        IR activicity for vibration modes
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_vibrations_infrared_activity'))

    x_turbomole_vibrations_raman_activity = Quantity(
        type=bool,
        shape=['x_turbomole_vibrations_num_modes'],
        description='''
        Raman activicity for vibration modes
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_vibrations_raman_activity'))

    x_turbomole_vibrations_intensities = Quantity(
        type=np.dtype(np.float64),
        shape=['x_turbomole_vibrations_num_modes'],
        description='''
        IR Intensity for each vibrational normal mode
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_vibrations_intensities'))

    x_turbomole_vibrations_mode_energies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_turbomole_vibrations_num_modes'],
        description='''
        Excitation energy associated with the vibrational normal modes.
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_vibrations_mode_energies'))

    x_turbomole_vibrations_normal_modes = Quantity(
        type=np.dtype(np.float64),
        shape=['x_turbomole_vibrations_num_modes', 'number_of_atoms', 3],
        description='''
        Nuclear displacements for each vibrational normal mode
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_vibrations_normal_modes'))

    x_turbomole_vibrations_num_modes = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of vibrational normal modes
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_vibrations_num_modes'))

    x_turbomole_virial_theorem = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Final value from the virial theorem
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_virial_theorem'))

    x_turbomole_wave_func_norm = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Final Wave Function Norm
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_wave_func_norm'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_turbomole_nodename = Quantity(
        type=str,
        shape=[],
        description='''
        compute node
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_nodename'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_turbomole_pceem_charges = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Charges of the point charges in the unit cell used by the PCEEM embedding model
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_pceem_charges'))

    x_turbomole_pceem_max_multipole = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Maximum multipole moment used in the PCEEM embedding
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_pceem_max_multipole'))

    x_turbomole_pceem_min_separation_cells = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Minimum separation between cells in PCEEM embedding for periodic fast multipole
        treatment
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_pceem_min_separation_cells'))

    x_turbomole_pceem_multipole_precision = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Multipole precision parameter for PCEEM embedding
        ''',
        a_legacy=LegacyDefinition(name='x_turbomole_pceem_multipole_precision'))


m_package.__init_metainfo__()
