from builtins import object
import setup_paths
from nomadcore.simple_parser import SimpleMatcher, mainFunction
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import os, sys, json, logging
SM=SimpleMatcher

##########################################################
#							 #
#		Parser for ORCA output file		 #
#							 #
##########################################################

##########################################################
###############[1] transfer PARSER CONTEXT ###############
##########################################################

logger = logging.getLogger("nomad.orcaParser")

class OrcaContext(object):
    """context for the sample parser"""

    def __init__(self):
        self.parser = None

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse different files"""
        pass

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser

        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv

        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()

##########################################################
############    [2] MAIN PARSER STARTS HERE   ############
##########################################################

def build_OrcaMainFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the main output file of ORCA.

    First, several subMatchers are defined, which are then used to piece together
    the final SimpleMatcher.
    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses main file of ORCA. 
    """
#
# a) SimpleMatcher for header and ORCA version:
# *********************************************
#
    return SM(
        name = 'root',
        weak = True,
        startReStr = r"\s*\* O  R  C  A \*\s",
        forwardMatch = True,
        sections = ["section_run"],
        subMatchers = [
            SM(name = 'ProgramHeader',
               startReStr = r"\s*\* O  R  C  A \*\s",
               subMatchers = [
                    SM(r"\s*Program Version\s*(?P<program_version>[0-9a-zA-Z_.]*)"),
                    SM(r" *\(\$Date\: *(?P<orca_program_compilation_date>[0-9/]+) at (?P<orca_program_compilation_time>[0-9:]+)")
                    ]),
            buildSinglePointMatcher()
            ])

def buildSinglePointMatcher():
#
# b) SimpleMatcher for Single Point Calculation:
# **********************************************
#
  return SM(name = 'SinglePointCalculation',
     startReStr = r"\s*\* Single Point Calculation \*\s",
     sections = ["section_system", "section_single_configuration_calculation"],
     subMatchers = [
       # Get atomic positions:
       SM(name = 'Atomic Coordinates',
          startReStr = r"CARTESIAN COORDINATES (ANGSTROEM)\s*",
          sections = ["x_orca_atom_positions"],
          subMatchers = [
          SM(r"\s+(?P<x_orca_atom_labels>[a-zA-Z]+)\s+(?P<x_orca_atom_positions_x>[-+0-9.]+)\s+(?P<orca_atom_positions_y>[-+0-9.]+)\s+(?P<x_orca_atom_positions_z>[-+0-9.]+)", repeats = True)
          ] 
       ),  
       # Get basis set information:
       SM(name = 'Basis set information',
          startReStr = r"BASIS SET INFORMATION\s*",
          sections = ["x_orca_basis_set_info"],
          subMatchers = [
          # Atom labels and basis set:
          SM(r"\s*Group\s+[0-9]+\s+Type\s+(?P<x_orca_atom_labels>[a-zA-Z]+)\s+:\s+(?P<x_orca_basis_set>[0-9a-z]+)\s+contracted\s+to\s+(?P<x_orca_basis_set_contracted>[0-9a-z]+)\s+pattern\s+\{[0-9/]+\}", repeats = True),
          # Auxiliary basis set information:
          SM(name = 'Auxiliary basis set information',
             startReStr = r"AUXILIARY BASIS SET INFORMATION\s*",
             sections = ["x_orca_auxiliary_basis_set_info"],
             subMatchers = [
             SM(r"\s*Group\s+[0-9]+\s+Type\s+(?P<x_orca_atom_labels>[a-zA-Z]+)\s+:\s+(?P<x_orca_auxiliary_basis_set>[0-9a-z]+)\s+contracted\s+to\s+(?P<x_orca_auxiliary_basis_set_contracted>[0-9a-z]+)\s+pattern\s+\{[0-9/]+\}", repeats = True)
             ]
          )
          ]   
       ),
       # Basis set statistics and startup info:
       SM(name = 'Basis set statistics and startup info',
          startReStr = r"\s*BASIS SET STATISTICS AND STARTUP INFO\s*",
          sections = ["x_orca_basis_set_statistics_and_startup_info"],
          subMatchers = [
          # Gaussian basis set:
          SM(name = 'Gaussian basis set',
             startReStr = r"\s*Gaussian basis set:\s*",
             subMatchers = [
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<orca_x_nb_of_primitive_gaussian_shells>[-+0-9.eEdD]+)"),
             SM(r"\s+# of primitive gaussian functions\s+\.\.\.\s+(?P<orca_x_nb_of_primitive_gaussian_functions>[-+0-9.eEdD]+)"),
             SM(r"\s+# of contracted shells\s+\.\.\.\s+(?P<orca_x_nb_of_contracted_shells>[-+0-9.eEdD]+)"),
             SM(r"\s+# of contracted basis functions\s+\.\.\.\s+(?P<orca_x_nb_of_contracted_basis_functions>[-+0-9.eEdD]+)"),
             SM(r"\s+Highest angular momentum\s+\.\.\.\s+(?P<orca_x_highest_angular_moment>[-+0-9.eEdD]+)"),
             SM(r"\s+Maximum contraction depth\s+\.\.\.\s+(?P<orca_x_maximum_contraction_depth>[-+0-9.eEdD]+)"),
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<orca_x_primitive_gaussian_shells>[-+0-9.eEdD]+)"),
             ]
          ),
          # Gaussian auxiliary basis set:
          SM(name = 'Gaussian auxiliary basis set',
             startReStr = r"\s*Auxiliary gaussian basis set:\s*",
             subMatchers = [
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<orca_x_nb_of_primitive_gaussian_shells_aux>[-+0-9.eEdD]+)"),
             SM(r"\s+# of primitive gaussian functions\s+\.\.\.\s+(?P<orca_x_nb_of_primitive_gaussian_functions_aux>[-+0-9.eEdD]+)"),
             SM(r"\s+# of contracted shells\s+\.\.\.\s+(?P<orca_x_nb_of_contracted_shells_aux>[-+0-9.eEdD]+)"),
             SM(r"\s+# of contracted aux-basis functions\s+\.\.\.\s+(?P<orca_x_nb_of_contracted_basis_functions_aux>[-+0-9.eEdD]+)"),
             SM(r"\s+Highest angular momentum\s+\.\.\.\s+(?P<orca_x_highest_angular_moment_aux>[-+0-9.eEdD]+)"),
             SM(r"\s+Maximum contraction depth\s+\.\.\.\s+(?P<orca_x_maximum_contraction_depth_aux>[-+0-9.eEdD]+)"),
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<orca_x_primitive_gaussian_shells_aux>[-+0-9.eEdD]+)"),
             ]
          )
          ]
       ),
       # SCF Settings:
       SM(name = 'Orca SCF settings',
          startReStr = r"\s*ORCA SCF\s*",
          sections = ["x_orca_scf_settings"],
          subMatchers = [
          SM(r"\s+Density Functional\s+Method\s+\.\.\.\s+(?P<orca_x_dft_method>[a-zA-Z()]+)"),
          SM(r"\s+Exchange Functional\s+Exchange\s+\.\.\.\s+(?P<orca_x_exchange_functional>[a-zA-Z0-9]+)"),
          SM(r"\s+X-Alpha parameter\s+XAlpha\s+\.\.\.\s+(?P<orca_x_xalpha_param>[-+0-9.eEdD]+)"),
          SM(r"\s+Becke's b parameter\s+XBeta\s+\.\.\.\s+(?P<orca_x_beckes_beta_param>[-+0-9.eEdD]+)"),
          SM(r"\s+Correlation Functional\s+Correlation\s+\.\.\.\s+(?P<orca_x_correl_functional>[a-zA-Z0-9]+)"),
          SM(r"\s+LDA part of GGA corr\.\s+LDAOpt\s+\.\.\.\s+(?P<orca_x_lda_part_of_gga_corr>[a-zA-Z-+0-9]+)"),
          SM(r"\s+Scalar relativistic method\s+\.\.\.\s+(?P<orca_x_scalar_relativistic_method>[a-zA-Z-+0-9]+)"),
          SM(r"\s+Speed of light used\s+Velit\s+\.\.\.\s+(?P<orca_x_speed_of_light_used>[0-9.]+)"),
          SM(r"\s+Hartree-Fock type\s+HFTyps+\.\.\.\s+(?P<orca_x_hf_type>[a-zA-Z]+)"),
          SM(r"\s+Total Charge\s+Charge\s+\.\.\.\s+(?P<orca_x_total_charge>[-+0-9.eEdD]+)"),
          SM(r"\s+Multiplicity\s+Mult\s+\.\.\.\s+(?P<orca_x_multiplicity>[-+0-9.eEdD]+)"),
          SM(r"\s+Number of Electrons\s+NEL\s+\.\.\.\s+(?P<orca_x_nelectrons>[-+0-9.eEdD]+)"),
          SM(r"\s+Nuclear Repulsion\s+ENuc\s+\.\.\.\s+(?P<orca_x_nuclear_repulsion__hartree>[-+0-9.eEdD]+)"),
          # Convergence Tolerance:
          SM(r"\s+Convergence Check Mode ConvCheckMode\s*\.\.\.\s+(?P<orca_x_convergence_check_mode>[a-zA-Z-+0-9.eEdD_-]+)"),
          SM(r"\s+Energy Change\s+TolE\s*\.\.\.\s+(?P<orca_x_energy_change_tolerance__hartree>[-+0-9.eEdD]+)"),
          SM(r"\s+1-El\. energy change\s*\.\.\.\s+(?P<orca_x_1_elect_energy_change__hartree>[-+0-9.eEdD]+)"),
          # DFT Grid generation:
          SM(r"\s+General Integration Accuracy\s+IntAcc\s*\.\.\.\s+(?P<orca_x_gral_integ_accuracy>[-+0-9.eEdD]+)"),
          SM(r"\s+Radial Grid Type\s+RadialGrid\s*\.\.\.\s+(?P<orca_x_radial_grid_type>[a-zA-Z-_]+)"),
          SM(r"\s+Angular Grid \(max\. acc\.\)\s+AngularGrid\s*\.\.\.\s+(?P<orca_x_angular_grid>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Angular grid pruning method\s+GridPruning\s*\.\.\.\s+(?P<orca_x_grid_pruning_method>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Weight generation scheme\s*WeightScheme\s*\.\.\.\s+(?P<orca_x_weight_gener_scheme>[a-zA-Z0-9]+)"),
          SM(r"\s+Basis function cutoff\s+BFCut\s*\.\.\.\s+(?P<orca_x_basis_fn_cutoff>[-+0-9.eEdD]+)"),
          SM(r"\s+Integration weight cutoff\s+WCut\s*\.\.\.\s+(?P<orca_x_integr_weight_cutoff>[-+0-9.eEdD]+)"),
          SM(r"\s+# of grid points \(after initial pruning\)\s+WCut\s*\.\.\.\s+(?P<orca_x_nb_grid_pts_after_initial_pruning>[-+0-9.eEdD]+)"),
          SM(r"\s+# of grid points \(after weights\+screening\)\s*\.\.\.\s+(?P<orca_x_nb_grid_pts_after_weights_screening>[-+0-9.eEdD]+)"),
          SM(r"\s+Total number of grid points\s*\.\.\.\s+(?P<orca_x_total_nb_grid_pts>[-+0-9.eEdD]+)"),
          SM(r"\s+Total number of batches\s*\.\.\.\s+(?P<orca_x_total_nb_batches>[-+0-9.eEdD]+)"),
          SM(r"\s+Average number of points per batch\s*\.\.\.\s+(?P<orca_x_avg_nb_points_per_batch>[-+0-9.eEdD]+)"),
          SM(r"\s+Average number of grid points per atom\s*\.\.\.\s+(?P<orca_x_avg_nb_grid_pts_per_atom>[-+0-9.eEdD]+)")
          ]
       ),
       # SCF iterations:
       SM(name = 'Orca SCF iterations',
          startReStr = r"\s*SCF ITERATIONS\s*",
          sections = ["section_scf_iteration"],
          subMatchers = [
          SM(r"\s*(?P<orca_x_iteration_nb>[0-9]+)\s+(?P<orca_x_total_energy__hartree>[-+0-9.eEdD]+)", repeats = True),
          ]
       ),
#      *****************************************************
#      *                     SUCCESS                       *
#      *           SCF CONVERGED AFTER  XY CYCLES          *
#      *****************************************************
       SM(name = 'Final step after convergence',
          startReStr = r"Setting up the final grid:",
          sections = ["x_orca_final_run_after_convergence"],
          subMatchers = [
          # Final DFT Grid generation:
          SM(r"\s+General Integration Accuracy\s+IntAcc\s*\.\.\.\s+(?P<orca_x_gral_integ_accuracy_final>[-+0-9.eEdD]+)"),
          SM(r"\s+Radial Grid Type\s+RadialGrid\s*\.\.\.\s+(?P<orca_x_radial_grid_type_final>[a-zA-Z-_]+)"),
          SM(r"\s+Angular Grid \(max\. acc\.\)\s+AngularGrid\s*\.\.\.\s+(?P<orca_x_angular_grid_final>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Angular grid pruning method\s+GridPruning\s*\.\.\.\s+(?P<orca_x_grid_pruning_method_final>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Weight generation scheme\s*WeightScheme\s*\.\.\.\s+(?P<orca_x_weight_gener_scheme_final>[a-zA-Z0-9]+)"),
          SM(r"\s+Basis function cutoff\s+BFCut\s*\.\.\.\s+(?P<orca_x_basis_fn_cutoff_final>[-+0-9.eEdD]+)"),
          SM(r"\s+Integration weight cutoff\s+WCut\s*\.\.\.\s+(?P<orca_x_integr_weight_cutoff_final>[-+0-9.eEdD]+)"),
          SM(r"\s+# of grid points \(after initial pruning\)\s+WCut\s*\.\.\.\s+(?P<orca_x_nb_grid_pts_after_initial_pruning_final>[-+0-9.eEdD]+)"),
          SM(r"\s+# of grid points \(after weights\+screening\)\s*\.\.\.\s+(?P<orca_x_nb_grid_pts_after_weights_screening_final>[-+0-9.eEdD]+)"),
          SM(r"\s+Total number of grid points\s*\.\.\.\s+(?P<orca_x_total_nb_grid_pts_final>[-+0-9.eEdD]+)"),
          SM(r"\s+Total number of batches\s*\.\.\.\s+(?P<orca_x_total_nb_batches_final>[-+0-9.eEdD]+)"),
          SM(r"\s+Average number of points per batch\s*\.\.\.\s+(?P<orca_x_avg_nb_points_per_batch_final>[-+0-9.eEdD]+)"),
          SM(r"\s+Average number of grid points per atom\s*\.\.\.\s+(?P<orca_x_avg_nb_grid_pts_per_atom_final>[-+0-9.eEdD]+)")
          ]
       ),
       # Final SCF total Energy:
       SM(name = 'Total Energy',
          startReStr = r"\s*TOTAL SCF ENERGY\s*",
          sections = ["x_orca_total_energy"],
          subMatchers = [
          SM(r"\s*Total Energy\s+:\s+(?P<orca_x_total_energy__hartree>[-+0-9.eEdD]+)"),
          # Energy Components:
          SM(name = 'Energy Components',
             startReStr = r"\s*Components:\s*",
             sections = ["x_orca_energy_componets"],
             subMatchers = [
             SM(r"\s*Nuclear Repulsion\s*:\s+(?P<orca_x_nuc_repulsion__hartree>[-+0-9.eEdD]+)"),
             SM(r"\s*Electronic Energy\s*:\s+(?P<orca_x_elec_energy__hartree>[-+0-9.eEdD]+)"),
             SM(r"\s*One Electron Energy:\s+(?P<orca_x_one_elec_energy__hartree>[-+0-9.eEdD]+)"),
             SM(r"\s*Two Electron Energy:\s+(?P<orca_x_two_elec_energy__hartree>[-+0-9.eEdD]+)"),
             # Virial Components:
             SM(r"\s*Potential Energy\s*:\s+(?P<orca_x_potential_energy__hartree>[-+0-9.eEdD]+)"),
             SM(r"\s*Kinetic Energy\s*:\s+(?P<orca_x_kinetc_energy__hartree>[-+0-9.eEdD]+)"),
             SM(r"\s*Virial Ratio\s*:\s+(?P<orca_x_virial_ratio>[-+0-9.eEdD]+)"),
             # DFT Components:
             SM(r"\s*N\(Alpha\)\s*:\s+(?P<orca_x_nb_elect_alpha_channel>[-+0-9.eEdD]+)"),
             SM(r"\s*N\(Beta\)\s*:\s+(?P<orca_x_nb_elect_beta_channel>[-+0-9.eEdD]+)"),
             SM(r"\s*N\(Total\)\s*:\s+(?P<orca_x_nb_elect_total>[-+0-9.eEdD]+)"),
             SM(r"\s*E\(X\)\s*:\s+(?P<orca_x_exchange_energy__hartree>[-+0-9.eEdD]+)"),
             SM(r"\s*E\(C\)\s*:\s+(?P<orca_x_correlation_energy__hartree>[-+0-9.eEdD]+)"),
             SM(r"\s*E\(XC\)\s*:\s+(?P<orca_x_exchange_correlation_energy__hartree>[-+0-9.eEdD]+)")
             ]
          ),
          # Final SCF convergence:
          SM(r"\s*Last Energy change\s+\.\.\.\s+(?P<orca_x_last_energy_change__hartree>[-+0-9.eEdD]+)\s+Tolerance :\s*(?P<orca_x_last_energy_change_tolerance__hartree>[-+0-9.eEdD]+)"),
          SM(r"\s*Last MAX-Density change\s+\.\.\.\s+(?P<orca_x_last_max_density_change__hartree>[-+0-9.eEdD]+)\s+Tolerance :\s*(?P<orca_x_last_max_density_tolerance__hartree>[-+0-9.eEdD]+)"),
          SM(r"\s*Last RMS-Density change\s+\.\.\.\s+(?P<orca_x_last_rms_density_change__hartree>[-+0-9.eEdD]+)\s+Tolerance :\s*(?P<orca_x_last_rms_density_tolerance__hartree>[-+0-9.eEdD]+)")
          ]
       ),
       # Orbitals Energies:
       SM(name = 'Orbital Energies',
          startReStr = r"\s*ORBITAL ENERGIES\s*",
          sections = ["x_orca_orbital_energies"],
          subMatchers = [
          SM(r"\s*(?P<orca_x_orbital_nb>[0-9]+)\s+(?P<orca_x_orbital_occupation_nb>[-+0-9]+)\s+(?P<orca_x_orbital_energy__hartree>)", repeats = True)
          ]
       ),
       # Mulliken population analysis:
       SM(name = 'Mulliken population analysis',
          startReStr = r"\s*\* MULLIKEN POPULATION ANALYSIS \*\s*",
          sections = ["x_orca_mulliken_analysis"],
          subMatchers = [
          SM(r"\s*(?P<orca_x_atom_nb>[0-9]+)\s+(?P<orca_x_atom_species>[a-zA-Z]+):\s+(?P<orca_x_mulliken_atom_charge>)", repeats = True),
          SM(r"\s*Sum of atomic charges:\s*(?P<orca_x_mulliken_total_charge>[-+0-9.eEdD]+)"),
          # Mulliken reduced orbital charges (mroc):
          SM(r"\s*(?P<orca_x_atom_nb_mroc>[0-9]+)\s+(?P<orca_x_atom_species_mroc>[a-zA-Z]+)(?P<orca_x_atom_orbital_mroc>[-+0-9a-zA-Z]+)\s*:\s*(?P<orca_x_mulliken_partial_orbital_charge_mroc>[-+0-9.eEdD]+)", repeats = True)
         ]
       ),
       # Time table:
       SM(name = 'timings',
          startReStr = r"TIMINGS",
          sections = ["x_orca_timings"],
          subMatchers = [
          SM(r"\s*Total SCF time:\s+(?P<orca_x_total_days_time>[0-9]+) days (?P<orca_x_total_hours_time>[0-9]+) hours (?P<orca_x_total_mins_time>[0-9]+) min (?P<orca_x_total_secs_time>[0-9]+) sec"),
          SM(r"\s*Total time\s*\.\.\.\.\s*(?P<orca_x_final_time>[0-9.]+) sec"),
          SM(r"\s*Sum of individual times\s*\.\.\.\.\s*(?P<orca_x_sum_individual_times>[0-9.]+) sec"),
          SM(r"\s*Fock matrix formation\s*\.\.\.\.\s*(?P<orca_x_fock_matrix_formation>[0-9.]+) sec"),
          SM(r"\s*Coulomb formation\s*\.\.\.\.\s*(?P<orca_x_coulomb_formation>[0-9.]+) sec"),
          SM(r"\s*Split-RI-J\s*\.\.\.\.\s*(?P<orca_x_split_rj>[0-9.]+) sec"),
          SM(r"\s*XC integration\s*\.\.\.\.\s*(?P<orca_x_xc_integration>[0-9.]+) sec"),
          SM(r"\s*Basis function eval\.\s*\.\.\.\.\s*(?P<orca_x_basis_fn_evaluation>[0-9.]+) sec"),
          SM(r"\s*Density eval\.\s*\.\.\.\.\s*(?P<orca_x_density_evaluation>[0-9.]+) sec"),
          SM(r"\s*XC-Functional eval\.\s*\.\.\.\.\s*(?P<orca_x_xc_functional_evaluation>[0-9.]+) sec"),
          SM(r"\s*XC-Potential eval\.\s*\.\.\.\.\s*(?P<orca_x_potential_evaluation>[0-9.]+) sec"),
          SM(r"\s*Diagonalization\s*\.\.\.\.\s*(?P<orca_x_diagonalization>[0-9.]+) sec"),
          SM(r"\s*Density matrix formation\s*\.\.\.\.\s*(?P<orca_x_density_matrix_formation>[0-9.]+) sec"),
          SM(r"\s*Population analysis\s*\.\.\.\.\s*(?P<orca_x_population_analysis>[0-9.]+) sec"),
          SM(r"\s*Initial guess\s*\.\.\.\.\s*(?P<orca_x_initial_guess>[0-9.]+) sec"),
          SM(r"\s*Orbital Transformation\s*\.\.\.\.\s*(?P<orca_x_orbital_transformation>[0-9.]+) sec"),
          SM(r"\s*Orbital Orthonormalization\s*\.\.\.\.\s*(?P<orca_x_orbital_orthonormalization>[0-9.]+) sec"),
          SM(r"\s*DIIS solution\s*\.\.\.\.\s*(?P<orca_x_diis_solution>[0-9.]+) sec"),
          SM(r"\s*Grid generation\s*\.\.\.\.\s*(?P<orca_x_grid_generation>[0-9.]+) sec")
          ]
       ),
       # Here new stuff:
#       SM(name = '',
#          startReStr = r"",
#          sections = [],
#          subMatchers = [
#          SM(),
#          SM(),
#          SM()
#          ]
#       )
     ]
  )
#
# c) SimpleMatcher for Calculation Including Atomic Optimizations:
# ****************************************************************
#

mainFileDescription = build_OrcaMainFileSimpleMatcher()



#
# This is the rest of Fawzi's code:
# *********************************
#
# loading metadata from nomad-meta-info/meta_info/nomad_meta_info/fhi_aims.nomadmetainfo.json

parserInfo = {
  "name": "sample_parser",
  "version": "1.0"
}

metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/orca.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)

if __name__ == "__main__":
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo, superContext=OrcaContext())
