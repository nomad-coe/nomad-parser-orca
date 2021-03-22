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
from orcaparser import OrcaParser


@pytest.fixture(scope='module')
def parser():
    return OrcaParser()


def test_scf(parser):
    archive = EntryArchive()
    parser.parse('tests/data/CO_scf/orca3.2985087.out', archive, None)

    assert archive.section_run[0].program_version == '3.0.3 - RELEASE   -'

    sec_method = archive.section_run[0].section_method[0]
    assert sec_method.electronic_structure_method == 'DFT'
    assert sec_method.x_orca_nelectrons == 14.0
    assert pytest.approx(sec_method.x_orca_energy_change_tolerance, 4.35974472220717e-26)
    assert sec_method.x_orca_radial_grid_type == 'Gauss-Chebyshev'
    assert len(sec_method.section_XC_functionals) == 4
    assert sec_method.section_XC_functionals[2].XC_functional_name == 'GGA_C_LYP'

    sec_system = archive.section_run[0].section_system[0]
    assert sec_system.atom_labels == ['C', 'O']
    assert pytest.approx(sec_system.atom_positions[1][0].magnitude, 1.25e-10)

    assert len(archive.section_run[0].section_single_configuration_calculation) == 1
    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert len(sec_scc.section_basis_set) == 3
    assert sec_scc.section_basis_set[1].x_orca_basis_set == '11s6p2d1f'
    assert sec_scc.section_basis_set[2].x_orca_nb_of_primitive_gaussian_functions == 92
    assert pytest.approx(sec_scc.energy_total.magnitude, -4.94114851e-16)
    assert pytest.approx(sec_scc.x_orca_potential_energy, -9.84253575e-16)
    assert pytest.approx(sec_scc.x_orca_nb_elect_total, 14.000005402207)
    assert len(sec_scc.section_scf_iteration) == 8
    assert pytest.approx(sec_scc.section_scf_iteration[3].energy_total_scf_iteration.magnitude, -4.94114809e-16)
    assert pytest.approx(sec_scc.section_scf_iteration[-1].x_orca_last_max_density_change, 9.441463170411847e-22)
    assert np.shape(sec_scc.section_eigenvalues[0].eigenvalues_values) == (1, 1, 62)
    assert pytest.approx(sec_scc.section_eigenvalues[0].eigenvalues_values[0][0][28].magnitude, 6.53237991e-18)
    assert sec_scc.section_eigenvalues[0].eigenvalues_occupation[0][0][6] == 2.0
    assert len(sec_scc.section_dos) == 2
    assert sec_scc.section_dos[0].x_orca_mulliken_atom_charge == 0.131793
    assert sec_scc.section_dos[1].x_orca_mulliken_partial_orbital_charge_mroc[7] == 0.027488
    assert sec_scc.x_orca_diis_solution == 0.003


def test_geomopt(parser):
    archive = EntryArchive()
    parser.parse('tests/data/CHO_geomopt/orca3.2985006.out', archive, None)

    sec_run = archive.section_run[0]
    assert len(sec_run.section_method) == 6
    assert len(sec_run.section_system) == 6
    assert len(sec_run.section_single_configuration_calculation) == 6

    assert sec_run.section_method[2].x_orca_nb_grid_pts_after_weights_screening == 34298
    assert sec_run.section_method[4].x_orca_integr_weight_cutoff == 1e-14
    assert pytest.approx(sec_run.section_system[1].atom_positions[2][1].magnitude, 9.54068e-11)
    assert len(sec_run.section_single_configuration_calculation[0].section_scf_iteration) == 13
    assert pytest.approx(sec_run.section_single_configuration_calculation[-1].x_orca_elec_energy, -6.34048432e-16)


def test_spinpol(parser):
    archive = EntryArchive()
    parser.parse('tests/data/BO_spinpol/orca3.2984863.out', archive, None)

    assert archive.section_run[0].section_method[0].x_orca_multiplicity == 2
    sec_eig = archive.section_run[0].section_single_configuration_calculation[0].section_eigenvalues[0]
    assert np.shape(sec_eig.eigenvalues_values) == (2, 1, 28)
    assert pytest.approx(sec_eig.eigenvalues_values[1][0][22].magnitude, 7.57745431e-18)
    assert sec_eig.eigenvalues_occupation[0][0][2] == 1.0
    sec_dos = archive.section_run[0].section_single_configuration_calculation[0].section_dos
    assert sec_dos[0].x_orca_mulliken_atom_charge == -0.01143
    assert sec_dos[1].x_orca_mulliken_partial_orbital_charge_mroc[2] == 1.450488


def test_ci(parser):
    archive = EntryArchive()
    parser.parse('tests/data/FeMgO_ci/orca3.2713636.out', archive, None)

    sec_method = archive.section_run[0].section_method[1]
    assert sec_method.electronic_structure_method == 'CCSD'
    assert sec_method.x_orca_single_excitations_on_off == 'ON'
    assert pytest.approx(sec_method.x_orca_t1_diagnostic, 6.77481921e-20)

    sec_system = archive.section_run[0].section_system[0]
    assert sec_system.atom_labels == ['Fe'] + ['O'] * 6 + ['Mg'] * 18 + ['Q'] * 704
    assert pytest.approx(sec_system.atom_positions[39][1].magnitude, -4.211228e-10)

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert pytest.approx(sec_scc.x_orca_ccsd_total_energy, -1.70953359e-14)


def test_tddft(parser):
    archive = EntryArchive()
    parser.parse('tests/data/ClTi_tddft/orca3.2706823.out', archive, None)

    assert archive.section_run[0].section_method[1].electronic_structure_method == 'TDDFT'
    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert sec_scc.section_excited_states[0].x_orca_excitation_energy[8] == 3956596900.0
    assert sec_scc.section_excited_states[0].x_orca_transition_dipole_moment_y[21] == -0.00035
