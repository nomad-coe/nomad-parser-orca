import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import common
from nomad.datamodel.metainfo import public

m_package = Package(
    name='orca_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='orca.nomadmetainfo.json'))


class x_orca_atom_positions(MSection):
    '''
    positions of an atom
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_orca_atom_positions'))

    x_orca_atom_labels = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_labels'))

    x_orca_atom_positions_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_positions_x'))

    x_orca_atom_positions_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_positions_y'))

    x_orca_atom_positions_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_positions_z'))


class x_orca_final_geometry(MSection):
    '''
    final optimized positions of an atom
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_orca_final_geometry'))

    x_orca_atom_labels_geo_opt = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_labels_geo_opt'))

    x_orca_atom_positions_x_geo_opt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_positions_x_geo_opt'))

    x_orca_atom_positions_y_geo_opt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_positions_y_geo_opt'))

    x_orca_atom_positions_z_geo_opt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_positions_z_geo_opt'))


class x_orca_section_functionals(MSection):
    '''
    XC functional
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_orca_section_functionals'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_orca_1_elect_energy_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_1_elect_energy_change'))

    x_orca_angular_grid = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_angular_grid'))

    x_orca_atomic_orbital_integral_source = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atomic_orbital_integral_source'))

    x_orca_avg_nb_grid_pts_per_atom = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_avg_nb_grid_pts_per_atom'))

    x_orca_avg_nb_points_per_batch = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_avg_nb_points_per_batch'))

    x_orca_basis_fn_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_basis_fn_cutoff'))

    x_orca_beckes_beta_param = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_beckes_beta_param'))

    x_orca_Brueckner_orbitals_calc_on_off = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_Brueckner_orbitals_calc_on_off'))

    x_orca_ccsdt_aaa_triples_contribution = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsdt_aaa_triples_contribution'))

    x_orca_ccsdt_aab_triples_contribution = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsdt_aab_triples_contribution'))

    x_orca_ci_half_s_and_s_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ci_half_s_and_s_energy'))

    x_orca_ci_iteration_nb = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ci_iteration_nb'))

    x_orca_convergence_check_mode = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_convergence_check_mode'))

    x_orca_convergence_tol_max_residuum = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_convergence_tol_max_residuum'))

    x_orca_correl_functional = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_correl_functional'))

    x_orca_coulomb_transformation_dimension_basis = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_coulomb_transformation_dimension_basis'))

    x_orca_coulomb_transformation_type = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_coulomb_transformation_type'))

    x_orca_energy_change_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_energy_change_tolerance'))

    x_orca_exchange_functional = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_exchange_functional'))

    x_orca_f12_correction_on_off = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_f12_correction_on_off'))

    x_orca_frozen_core_treatment = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_frozen_core_treatment'))

    x_orca_gral_integ_accuracy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_gral_integ_accuracy'))

    x_orca_grid_pruning_method = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_grid_pruning_method'))

    x_orca_hf_method = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_hf_method'))

    x_orca_hf_type = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_hf_type'))

    x_orca_integr_weight_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_integr_weight_cutoff'))

    x_orca_integral_package_used = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_integral_package_used'))

    x_orca_integral_transformation = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_integral_transformation'))

    x_orca_K_C_formation = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_K_C_formation'))

    x_orca_lda_part_of_gga_corr = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_lda_part_of_gga_corr'))

    x_orca_level_shift_amplitude_update = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_level_shift_amplitude_update'))

    x_orca_mp2_aux_basis_dimension = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_aux_basis_dimension'))

    x_orca_mp2_basis_dimension = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_basis_dimension'))

    x_orca_mp2_initial_guess = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_initial_guess'))

    x_orca_multiplicity = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_multiplicity'))

    x_orca_nb_alpha_pairs_included = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_alpha_pairs_included'))

    x_orca_nb_beta_pairs_included = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_beta_pairs_included'))

    x_orca_nb_grid_pts_after_initial_pruning = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_grid_pts_after_initial_pruning'))

    x_orca_nb_grid_pts_after_weights_screening = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_grid_pts_after_weights_screening'))

    x_orca_nb_internal_alpha_mol_orbitals = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_internal_alpha_mol_orbitals'))

    x_orca_nb_internal_beta_mol_orbitals = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_internal_beta_mol_orbitals'))

    x_orca_nb_of_atomic_orbitals = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_atomic_orbitals'))

    x_orca_nb_of_correlated_electrons = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_correlated_electrons'))

    x_orca_nb_of_electrons = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_electrons'))

    x_orca_nelectrons = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nelectrons'))

    x_orca_nuclear_repulsion = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nuclear_repulsion'))

    x_orca_orbital_opt_on_off = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_orbital_opt_on_off'))

    x_orca_pair_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_pair_cutoff'))

    x_orca_perturbative_triple_excitations_on_off = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_perturbative_triple_excitations_on_off'))

    x_orca_radial_grid_type = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_radial_grid_type'))

    x_orca_reference_wave_function = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_reference_wave_function'))

    x_orca_scalar_relativistic_method = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_scalar_relativistic_method'))

    x_orca_scaling_mp2_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_scaling_mp2_energy'))

    x_orca_single_excitations_on_off = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_single_excitations_on_off'))

    x_orca_single_norm_half_ss = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_single_norm_half_ss'))

    x_orca_speed_of_light_used = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_speed_of_light_used'))

    x_orca_total_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_charge'))

    x_orca_total_nb_batches = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_nb_batches'))

    x_orca_total_nb_grid_pts = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_nb_grid_pts'))

    x_orca_total_nb_pairs_included = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_nb_pairs_included'))

    x_orca_wave_function_correlation_treatment = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_wave_function_correlation_treatment'))

    x_orca_weight_gener_scheme = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_weight_gener_scheme'))

    x_orca_xalpha_param = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_xalpha_param'))

    x_orca_XC_functional_type = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_XC_functional_type'))

    x_orca_z_vector_calc_on_off = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_z_vector_calc_on_off'))

    x_orca_section_functionals = SubSection(
        sub_section=SectionProxy('x_orca_section_functionals'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_orca_section_functionals'))


class section_scf_iteration(public.section_scf_iteration):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_scf_iteration'))

    x_orca_angular_grid_final = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_angular_grid_final'))

    x_orca_avg_nb_grid_pts_per_atom_final = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_avg_nb_grid_pts_per_atom_final'))

    x_orca_avg_nb_points_per_batch_final = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_avg_nb_points_per_batch_final'))

    x_orca_basis_fn_cutoff_final = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_basis_fn_cutoff_final'))

    x_orca_gral_integ_accuracy_final = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_gral_integ_accuracy_final'))

    x_orca_grid_pruning_method_final = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_grid_pruning_method_final'))

    x_orca_integr_weight_cutoff_final = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_integr_weight_cutoff_final'))

    x_orca_iteration_nb = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_iteration_nb'))

    x_orca_last_energy_change_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_last_energy_change_tolerance'))

    x_orca_last_energy_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_last_energy_change'))

    x_orca_last_max_density_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_last_max_density_change'))

    x_orca_last_max_density_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_last_max_density_tolerance'))

    x_orca_last_rms_density_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_last_rms_density_change'))

    x_orca_last_rms_density_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_last_rms_density_tolerance'))

    x_orca_nb_grid_pts_after_initial_pruning_final = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_grid_pts_after_initial_pruning_final'))

    x_orca_nb_grid_pts_after_weights_screening_final = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_grid_pts_after_weights_screening_final'))

    x_orca_radial_grid_type_final = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_radial_grid_type_final'))

    x_orca_total_nb_batches_final = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_nb_batches_final'))

    x_orca_total_nb_grid_pts_final = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_nb_grid_pts_final'))

    x_orca_weight_gener_scheme_final = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_weight_gener_scheme_final'))


class section_dos(public.section_dos):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_dos'))

    x_orca_atom_nb_mroc = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_nb_mroc'))

    x_orca_atom_nb = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_nb'))

    x_orca_atom_orbital_mroc = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_orbital_mroc'))

    x_orca_atom_species_mroc = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_species_mroc'))

    x_orca_atom_species = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_atom_species'))

    x_orca_mulliken_atom_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mulliken_atom_charge'))

    x_orca_mulliken_partial_orbital_charge_mroc = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mulliken_partial_orbital_charge_mroc'))

    x_orca_mulliken_total_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mulliken_total_charge'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_orca_atom_positions = SubSection(
        sub_section=SectionProxy('x_orca_atom_positions'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_orca_atom_positions'))


class section_basis_set(public.section_basis_set):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_basis_set'))

    x_orca_auxiliary_basis_set_contracted = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_auxiliary_basis_set_contracted'))

    x_orca_auxiliary_basis_set = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_auxiliary_basis_set'))

    x_orca_basis_set_atom_labels = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_basis_set_atom_labels'))

    x_orca_basis_set_contracted = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_basis_set_contracted'))

    x_orca_basis_set = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_basis_set'))

    x_orca_highest_angular_moment_aux = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_highest_angular_moment_aux'))

    x_orca_highest_angular_moment = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_highest_angular_moment'))

    x_orca_maximum_contraction_depth_aux = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_maximum_contraction_depth_aux'))

    x_orca_maximum_contraction_depth = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_maximum_contraction_depth'))

    x_orca_nb_of_contracted_basis_functions_aux = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_contracted_basis_functions_aux'))

    x_orca_nb_of_contracted_basis_functions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_contracted_basis_functions'))

    x_orca_nb_of_contracted_shells_aux = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_contracted_shells_aux'))

    x_orca_nb_of_contracted_shells = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_contracted_shells'))

    x_orca_nb_of_primitive_gaussian_functions_aux = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_primitive_gaussian_functions_aux'))

    x_orca_nb_of_primitive_gaussian_functions = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_primitive_gaussian_functions'))

    x_orca_nb_of_primitive_gaussian_shells_aux = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_primitive_gaussian_shells_aux'))

    x_orca_nb_of_primitive_gaussian_shells = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_of_primitive_gaussian_shells'))

    x_orca_nb_primitive_gaussian_shells_aux = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_primitive_gaussian_shells_aux'))

    x_orca_nb_primitive_gaussian_shells = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_primitive_gaussian_shells'))


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    x_orca_basis_fn_evaluation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_basis_fn_evaluation'))

    x_orca_ccsd_correlation_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsd_correlation_energy'))

    x_orca_ccsd_final_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsd_final_energy'))

    x_orca_ccsd_t_final_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsd_t_final_energy'))

    x_orca_ccsd_total_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsd_total_energy'))

    x_orca_ccsdt_aba_triples_contribution = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsdt_aba_triples_contribution'))

    x_orca_ccsdt_bbb_triples_contribution = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsdt_bbb_triples_contribution'))

    x_orca_ccsdt_final_corr_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsdt_final_corr_energy'))

    x_orca_ccsdt_total_triples_correction = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ccsdt_total_triples_correction'))

    x_orca_ci_correl_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ci_correl_energy'))

    x_orca_ci_deltaE_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ci_deltaE_energy'))

    x_orca_ci_iteration_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ci_iteration_time'))

    x_orca_ci_residual_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ci_residual_energy'))

    x_orca_ci_total_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_ci_total_energy'))

    x_orca_correlation_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_correlation_energy'))

    x_orca_coulomb_formation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_coulomb_formation'))

    x_orca_density_evaluation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_density_evaluation'))

    x_orca_density_matrix_formation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_density_matrix_formation'))

    x_orca_diagonalization = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_diagonalization'))

    x_orca_diis_solution = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_diis_solution'))

    x_orca_elec_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_elec_energy'))

    x_orca_exchange_correlation_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_exchange_correlation_energy'))

    x_orca_exchange_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_exchange_energy'))

    x_orca_final_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_final_time'))

    x_orca_fock_matrix_formation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_fock_matrix_formation'))

    x_orca_geo_opt_cycle = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_geo_opt_cycle'))

    x_orca_grid_generation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_grid_generation'))

    x_orca_initial_guess = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_initial_guess'))

    x_orca_kinetc_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_kinetc_energy'))

    x_orca_mp2_corr_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_corr_energy'))

    x_orca_mp2_energy_spin_aa = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_energy_spin_aa'))

    x_orca_mp2_energy_spin_ab = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_energy_spin_ab'))

    x_orca_mp2_energy_spin_bb = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_energy_spin_bb'))

    x_orca_mp2_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_energy'))

    x_orca_mp2_total_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_mp2_total_energy'))

    x_orca_nb_elect_alpha_channel = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_elect_alpha_channel'))

    x_orca_nb_elect_beta_channel = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_elect_beta_channel'))

    x_orca_nb_elect_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nb_elect_total'))

    x_orca_nuc_repulsion = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_nuc_repulsion'))

    x_orca_one_elec_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_one_elec_energy'))

    x_orca_orbital_orthonormalization = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_orbital_orthonormalization'))

    x_orca_orbital_transformation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_orbital_transformation'))

    x_orca_population_analysis = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_population_analysis'))

    x_orca_potential_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_potential_energy'))

    x_orca_potential_evaluation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_potential_evaluation'))

    x_orca_split_rj = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_split_rj'))

    x_orca_sum_individual_times = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_sum_individual_times'))

    x_orca_t1_diagnostic = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_t1_diagnostic'))

    x_orca_T_and_T_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_T_and_T_energy'))

    x_orca_total_days_time = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_days_time'))

    x_orca_total_hours_time = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_hours_time'))

    x_orca_total_mins_time = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_mins_time'))

    x_orca_total_secs_time = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_total_secs_time'))

    x_orca_two_elec_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_two_elec_energy'))

    x_orca_virial_ratio = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_virial_ratio'))

    x_orca_xc_functional_evaluation = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_xc_functional_evaluation'))

    x_orca_xc_integration = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_xc_integration'))


class section_sampling_method(public.section_sampling_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_sampling_method'))

    x_orca_coords_choice_name = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_coords_choice_name'))

    x_orca_coords_choice = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_coords_choice'))

    x_orca_energy_change_tol_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_energy_change_tol_value'))

    x_orca_energy_change_tol = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_energy_change_tol'))

    x_orca_initial_hessian_name = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_initial_hessian_name'))

    x_orca_initial_hessian = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_initial_hessian'))

    x_orca_max_displacement_tol_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_max_displacement_tol_value'))

    x_orca_max_displacement_tol = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_max_displacement_tol'))

    x_orca_max_gradient_tol_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_max_gradient_tol_value'))

    x_orca_max_gradient_tol = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_max_gradient_tol'))

    x_orca_rms_displacement_tol_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_rms_displacement_tol_value'))

    x_orca_rms_displacement_tol = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_rms_displacement_tol'))

    x_orca_rms_gradient_tol_value = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_rms_gradient_tol_value'))

    x_orca_rms_gradient_tol = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_rms_gradient_tol'))

    x_orca_update_method_name = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_update_method_name'))

    x_orca_update_method = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_update_method'))

    x_orca_final_geometry = SubSection(
        sub_section=SectionProxy('x_orca_final_geometry'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_orca_final_geometry'))


class section_excited_states(common.section_excited_states):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_excited_states'))

    x_orca_excitation_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_excitation_energy'))

    x_orca_oscillator_strength = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_oscillator_strength'))

    x_orca_transition_dipole_moment_x = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_transition_dipole_moment_x'))

    x_orca_transition_dipole_moment_y = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_transition_dipole_moment_y'))

    x_orca_transition_dipole_moment_z = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_transition_dipole_moment_z'))


class section_eigenvalues(public.section_eigenvalues):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_eigenvalues'))

    x_orca_orbital_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_orbital_energy'))

    x_orca_orbital_nb = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_orbital_nb'))

    x_orca_orbital_occupation_nb = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_orbital_occupation_nb'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_orca_program_compilation_date = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_program_compilation_date'))

    x_orca_program_compilation_time = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_program_compilation_time'))

    x_orca_program_svn = Quantity(
        type=str,
        shape=[],
        description='''
        svn revision of the program
        ''',
        a_legacy=LegacyDefinition(name='x_orca_program_svn'))

    x_orca_program_version = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_orca_program_version'))


m_package.__init_metainfo__()
