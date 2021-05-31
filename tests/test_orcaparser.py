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


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


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
    assert sec_method.x_orca_energy_change_tolerance == approx(4.35974472e-24)
    assert sec_method.x_orca_radial_grid_type == 'Gauss-Chebyshev'
    assert len(sec_method.section_XC_functionals) == 4
    assert sec_method.section_XC_functionals[2].XC_functional_name == 'GGA_C_LYP'

    sec_system = archive.section_run[0].section_system[0]
    assert sec_system.atom_labels == ['C', 'O']
    assert sec_system.atom_positions[1][0].magnitude == approx(1.25e-10)

    assert len(archive.section_run[0].section_single_configuration_calculation) == 1
    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert len(sec_scc.section_basis_set) == 3
    assert sec_scc.section_basis_set[1].x_orca_basis_set == '11s6p2d1f'
    assert sec_scc.section_basis_set[2].x_orca_nb_of_primitive_gaussian_functions == 92
    assert sec_scc.energy_total.magnitude == approx(-4.94114851e-16)
    assert sec_scc.x_orca_potential_energy == approx(-9.84253575e-16)
    assert sec_scc.x_orca_nb_elect_total == approx(14.000005402207)
    assert len(sec_scc.section_scf_iteration) == 8
    assert sec_scc.section_scf_iteration[3].energy_total_scf_iteration.magnitude == approx(-4.94113458e-16)
    assert sec_scc.section_scf_iteration[-1].x_orca_last_max_density_change == approx(9.441463170411847e-22)
    assert np.shape(sec_scc.eigenvalues[0].band_energies[0].band_energies_values) == (62,)
    assert sec_scc.eigenvalues[0].band_energies[0].band_energies_values[28].magnitude == approx(6.53237991e-18)
    assert sec_scc.eigenvalues[0].band_energies[0].band_energies_occupations[6] == 2.0
    assert len(sec_scc.atom_charges[0].charges_total) == 2
    assert sec_scc.atom_charges[0].charges_total[0].charges_value.magnitude == 0.131793
    assert sec_scc.atom_charges[0].charges_partial[27].charges_value.magnitude == 0.027488
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
    assert sec_run.section_system[1].atom_positions[2][1].magnitude == approx(9.54068e-11)
    assert len(sec_run.section_single_configuration_calculation[0].section_scf_iteration) == 13
    assert sec_run.section_single_configuration_calculation[-1].x_orca_elec_energy == approx(-6.34048432e-16)


def test_spinpol(parser):
    archive = EntryArchive()
    parser.parse('tests/data/BO_spinpol/orca3.2984863.out', archive, None)

    assert archive.section_run[0].section_method[0].x_orca_multiplicity == 2
    sec_eig = archive.section_run[0].section_single_configuration_calculation[0].eigenvalues[0]
    assert np.shape(sec_eig.band_energies[1].band_energies_values) == (28,)
    assert sec_eig.band_energies[1].band_energies_values[22].magnitude == approx(7.57745431e-18)
    assert sec_eig.band_energies[0].band_energies_occupations[2] == 1.0
    sec_charges = archive.section_run[0].section_single_configuration_calculation[0].atom_charges[0]
    assert sec_charges.charges_total[0].charges_value.magnitude == -0.01143
    assert sec_charges.charges_partial[14].charges_value.magnitude == 1.450488


def test_ci(parser):
    archive = EntryArchive()
    parser.parse('tests/data/FeMgO_ci/orca3.2713636.out', archive, None)

    sec_method = archive.section_run[0].section_method[1]
    assert sec_method.electronic_structure_method == 'CCSD'
    assert sec_method.x_orca_single_excitations_on_off == 'ON'
    assert sec_method.x_orca_t1_diagnostic == approx(6.77481921e-20)

    sec_system = archive.section_run[0].section_system[0]
    assert sec_system.atom_labels == ['Fe'] + ['O'] * 6 + ['Mg'] * 18 + ['Q'] * 704
    assert sec_system.atom_positions[39][1].magnitude == approx(-4.211228e-10)

    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert sec_scc.x_orca_ccsd_total_energy == approx(-1.70953359e-14)


def test_tddft(parser):
    archive = EntryArchive()
    parser.parse('tests/data/ClTi_tddft/orca3.2706823.out', archive, None)

    assert archive.section_run[0].section_method[1].electronic_structure_method == 'TDDFT'
    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert sec_scc.section_excited_states[0].x_orca_excitation_energy[8] == 3956596900.0
    assert sec_scc.section_excited_states[0].x_orca_transition_dipole_moment_y[21] == -0.00035
