#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from turbomoleparser.turbomole_parser import TurbomoleParser


@pytest.fixture(scope='module')
def parser():
    return TurbomoleParser()


def test_aoforce(parser):
    archive = EntryArchive()
    parser.parse('tests/data/aoforce/vib.out', archive, None)

    assert archive.section_run[0].program_version == '7.2 ( 21285 )'
    assert archive.section_run[0].time_run_date_start.magnitude == 1532973127.689

    sec_method = archive.section_run[0].section_method[0]
    assert sec_method.electronic_structure_method == 'DFT'
    assert sec_method.section_XC_functionals[1].XC_functional_name == 'GGA_C_P86'
    assert len(sec_method.section_method_atom_kind) == 4
    sec_basis_set = sec_method.section_method_basis_set[0]
    assert sec_basis_set.number_of_basis_sets_atom_centered == 4
    assert len(sec_basis_set.mapping_section_method_basis_set_atom_centered) == 4
    assert sec_method.van_der_Waals_method == 'DFT-D3'
    assert sec_method.x_turbomole_controlIn_scf_conv == 8

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert sec_scc.energy_total.magnitude == pytest.approx(-3.58404386e-15)
    assert sec_scc.energy_zero_point == pytest.approx(1.02171533e-18)
    assert sec_scc.energy_current == pytest.approx(-3.58302215e-15)
    assert np.shape(sec_scc.hessian_matrix) == (31, 31, 3, 3)
    assert sec_scc.hessian_matrix[3][2][2][0] == pytest.approx(-38.1728237)
    assert np.shape(sec_scc.x_turbomole_vibrations_normal_modes) == (93, 31, 3)
    assert sec_scc.x_turbomole_vibrations_normal_modes[0][4][1] == pytest.approx(-0.02383)
    assert sec_scc.x_turbomole_vibrations_mode_energies[46] == pytest.approx(1005.73)
    assert sec_scc.x_turbomole_vibrations_intensities[72] == pytest.approx(41.61)
    assert sec_scc.x_turbomole_vibrations_infrared_activity[49]
    assert not sec_scc.x_turbomole_vibrations_raman_activity[5]

    sec_system = archive.section_run[0].section_system[0]
    assert len(sec_system.atom_positions) == 31
    assert sec_system.atom_positions[7][1].magnitude == pytest.approx(-9.34235013e-11)
    assert sec_system.atom_labels[21] == 'O'


def test_ccsdf12(parser):
    archive = EntryArchive()
    parser.parse('tests/data/ccsdf12.out', archive, None)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 3
    assert sec_sccs[0].energy_total.magnitude == pytest.approx(-2.99865659e-15)
    assert sec_sccs[1].energy_total.magnitude == pytest.approx(-2.99841134e-15)
    assert sec_sccs[2].energy_total.magnitude == pytest.approx(-2.99844594e-15)
    assert sec_sccs[1].energy_current.magnitude == pytest.approx(-5.78479974e-18)
    sec_scfs = sec_sccs[0].section_scf_iteration
    assert len(sec_scfs) == 13
    assert sec_scfs[8].energy_total_scf_iteration.magnitude == pytest.approx(-2.99844594e-15)
    assert sec_scfs[2].time_scf_iteration.magnitude == 2.09


def test_grad_statpt_dscf(parser):
    archive = EntryArchive()
    parser.parse('tests/data/acrolein_grad_statpt_dscf.out', archive, None)

    assert archive.section_run[0].section_basis_set_atom_centered[1].basis_set_atom_centered_short_name == 'def2-SVP'

    sec_methods = archive.section_run[0].section_method
    assert len(sec_methods) == 3
    assert sec_methods[0].section_XC_functionals[0].XC_functional_name == 'HYB_GGA_XC_B3LYP'

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 3
    assert sec_systems[1].atom_positions[5][1].magnitude == pytest.approx(1.22377337e-10,)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert sec_sccs[0].atom_forces_raw[6][0].magnitude == pytest.approx(-4.2984543e-12)
    sec_scfs = sec_sccs[2].section_scf_iteration
    assert len(sec_scfs) == 3
    assert sec_scfs[1].energy_total_scf_iteration.magnitude == pytest.approx(-8.35592725e-16)
    assert sec_scfs[0].x_turbomole_delta_eigenvalues == pytest.approx(2.92683961e-22)
    assert sec_sccs[2].electronic_kinetic_energy.magnitude == pytest.approx(8.27834082e-16)

    sec_sampling = archive.section_run[0].section_sampling_method[0]
    assert sec_sampling.x_turbomole_geometry_optimization_trustregion_min.magnitude == pytest.approx(5.29177211e-14)
    assert sec_sampling.geometry_optimization_method == 'BFGS'
    assert sec_sampling.geometry_optimization_threshold_force.magnitude == pytest.approx(8.2387235e-11)


def test_escf(parser):
    archive = EntryArchive()
    parser.parse('tests/data/benzene_escf.out', archive, None)

    sec_method = archive.section_run[0].section_method[0]
    assert sec_method.electronic_structure_method == 'G0W0'
    assert sec_method.x_turbomole_gw_eta_factor == pytest.approx(4.35974472e-21)
    assert sec_method.x_turbomole_gw_approximation == 'G0W0'

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    sec_gw_eigs = sec_scc.section_eigenvalues[0].x_turbomole_section_eigenvalues_GW[0]
    assert sec_gw_eigs.x_turbomole_eigenvalue_ks_GroundState[9].magnitude == pytest.approx(-3.59608546e-18)
    assert sec_gw_eigs.x_turbomole_eigenvalue_ExactExchange_perturbativeGW[1].magnitude == pytest.approx(-1.55874163e-17)
    assert sec_gw_eigs.x_turbomole_Z_factor[19] == 0.786


def test_freeh(parser):
    archive = EntryArchive()
    parser.parse('tests/data/freeh.out', archive, None)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 2
    assert sec_sccs[0].energy_zero_point.magnitude == pytest.approx(4.89692971e-19)
    assert sec_sccs[1].energy_correction_entropy.magnitude == pytest.approx(2.00144971e-19)
    assert sec_sccs[1].heat_capacity_C_v.magnitude == pytest.approx(2.27860167e-22)
    assert sec_sccs[1].pressure.magnitude == 100000.0


def test_pnoccsd(parser):
    archive = EntryArchive()
    parser.parse('tests/data/pnoccsd.out', archive, None)

    assert np.shape(archive.section_run[0].section_system[0].atom_positions) == (51, 3)

    sec_methods = archive.section_run[0].section_method
    assert len(sec_methods) == 4
    assert sec_methods[0].electronic_structure_method == 'CCSD(T)'
    assert sec_methods[1].electronic_structure_method == 'MP2'
    assert sec_methods[2].electronic_structure_method == 'CCSD'
    assert sec_methods[3].electronic_structure_method == 'CCSD(T0)'

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 4
    assert sec_sccs[0].energy_total.magnitude == pytest.approx(-5.63810959e-15)
    assert sec_sccs[1].energy_total.magnitude == pytest.approx(-5.63669838e-15)
    assert sec_sccs[2].energy_current.magnitude == pytest.approx(-2.19140251e-17)
    assert sec_sccs[3].energy_total.magnitude == pytest.approx(-5.6380984e-15)

    sec_scfs = sec_sccs[0].section_scf_iteration
    assert len(sec_scfs) == 13
    assert sec_scfs[6].energy_total_scf_iteration.magnitude == pytest.approx(-5.63708622e-15)


def test_ricc2(parser):
    archive = EntryArchive()
    parser.parse('tests/data/MgO_embedding_ricc2.out', archive, None)

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 4
    assert sec_systems[0].atom_positions[4][2].magnitude == pytest.approx(6.38760003e-10)
    assert sec_systems[1].atom_positions[18][0].magnitude == pytest.approx(4.25840002e-10)
    assert sec_systems[2].atom_positions[25][2].magnitude == pytest.approx(2.12920001e-10)
    assert sec_systems[3].atom_positions[-2][1].magnitude == pytest.approx(8.51680003e-10)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 3
    assert sec_sccs[1].energy_total.magnitude == pytest.approx(-8.69567802e-15)


def test_ridft(parser):
    archive = EntryArchive()
    parser.parse('tests/data/ridft.out', archive, None)

    sec_method = archive.section_run[0].section_method[0]
    assert sec_method.x_turbomole_dft_d3_version == '3.1 Rev 0'

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert sec_scc.section_energy_van_der_Waals[0].energy_van_der_Waals.magnitude == pytest.approx(-1.32811671e-18)
    assert sec_scc.energy_total.magnitude == pytest.approx(-2.25881721e-14)
    assert sec_scc.x_turbomole_virial_theorem == pytest.approx(1.94918952771)

    sec_scf = sec_scc.section_scf_iteration
    assert len(sec_scf) == 28
    assert sec_scf[3].x_turbomole_energy_2electron_scf_iteration.magnitude == pytest.approx(1.02566632e-13)
    assert sec_scf[23].energy_XC_scf_iteration.magnitude == pytest.approx(-2.28814098e-15)
