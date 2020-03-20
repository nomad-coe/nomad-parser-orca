# Copyright 2016-2018 Sebastian Alarcon Villaseca, Sebastián Alarcón Villaseca, Fawzi Mohamed, Micael Oliveira, Ankit Kariryaa, Danio Brambila
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# -*- coding: utf-8 -*
from builtins import object
import numpy as np
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
        self.CalculationGIndex = None

    def initialize_values(self):
        """allows to reset values if the same superContext is used to parse different files"""
        self.CalculationGIndex = None

    def startedParsing(self, path, parser):
        """called when parsing starts"""
        self.parser = parser

        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv

        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()

    def onClose_x_orca_atom_positions(self, backend, gIndex, value):
            x = value["x_orca_atom_positions_x"]
       	    y = value["x_orca_atom_positions_y"]
       	    z = value["x_orca_atom_positions_z"]
            pos = np.zeros((len(x),3), dtype=float)
            pos[:, 0] = x
       	    pos[:, 1] = y
       	    pos[:, 2] = z
            backend .addArrayValues("atom_positions", pos)
            backend.addValue("atom_labels", value["x_orca_atom_labels"])

    def onClose_x_orca_final_geometry(self, backend, gIndex, value):
            x = value["x_orca_atom_positions_x_geo_opt"]
            y = value["x_orca_atom_positions_y_geo_opt"]
            z = value["x_orca_atom_positions_z_geo_opt"]
            pos = np.zeros((len(x),3), dtype=float)
            pos[:,0] = x
            pos[:,1] = y
            pos[:,2] = z
            backend.addArrayValues("geometry_optimization_converged", pos)
            backend.addValue("atom_labels", value["x_orca_atom_labels_geo_opt"])

    def onClose_x_orca_section_functionals(self, backend, gIndex, section):
        functional_names = section["x_orca_XC_functional_type"]

        if functional_names == None: functional = "HF" #default method is Hartree-Fock
        else: functional = functional_names[-1]

        if functional:
            functionalMap = {

		#***************************************
		# Local functionals
		#***************************************
		# HFS: Hartree-Fock Slater
                "HF":     ["HF_X"],
                "HFS":    ["HF_X"],
		# XAlpha: The famous old Slater Xa theory
		"XAlpha": ["HF_X", "LDA_C_XALPHA"],
		# LSD: Local spin density (VWN-5A form)
		"LSD":    ["LDA_X", "LDA_C_VWN_RPA", "VWN5_RPA"],
                "LDA":    ["LDA_X", "LDA_C_VWN_RPA", "VWN5_RPA"],
		# VWN5: Local spin density (VWN-5)
		"VWN":    ["LDA_X", "LDA_C_VWN_RPA", "VWN5_RPA"],
                "VWN5":   ["LDA_X", "LDA_C_VWN_RPA", "VWN5_RPA"],
		# VWN3: Local spin density (VWN-3)
		"VWN3":   ["LDA_X", "LDA_C_VWN_3"],
		# PWLDA: Local spin density (PW-LDA)
		"PWLDA":  ["LDA_X", "LDA_C_PW"],
		#***************************************
		# ’Pure’ GGA functionals
		#***************************************
		# BNULL: Becke ’88 exchange, no corr.
		"BNULL":  ["GGA_X_B88"],
		# BVWN: Becke ’88 exchane, VWN-5 corr.
		"BVWN":   ["GGA_X_B88", "LDA_C_VWN_RPA", "VWN5_RPA"],
		# BP: Becke X-Perdew 86 correlation
		"BP":     ["LDA_X", "GGA_X_B88", "LDA_C_VWN", "GGA_C_P86"],
		"BP86":   ["LDA_X", "GGA_X_B88", "LDA_C_VWN", "GGA_C_P86"],
		# PW91: Perdew-Wang GGA-II ’91 func.
		"PW91":   ["GGA_X_PW91", "GGA_C_PW91"],
		# PWP: Perdew-Wang ’91 exchange and Perdew ’86 correlation
                "PWP":    ["GGA_X_PW91", "GGA_C_P86"],
		# mPWPW: Modified PW with PW correlation
		"mPWPW":  ["GGA_X_MPW91", "GGA_C_PW91"],
		# mPWLYP: same with LYP correlation
		"mPWLYP": ["GGA_X_MPW91", "GGA_C_LYP"],
		# BLYP: Becke X with LYP correlation
		"BLYP":   ["GGA_X_B88", "GGA_C_LYP"],
		# GP: Gill ’96 X, Perdew ’86 corr.
		"GP":     ["GGA_X_G96", "GGA_C_P86"],
		#GLYP: Gill ’96 X with LYP correlation
		"GLYP":   ["GGA_X_G96", "GGA_C_LYP"],
		# PBE: Perdew-Burke-Ernzerhof
		"PBE":    ["GGA_X_PBE", "GGA_C_PBE"],
		# revPBE: Revised PBE (exchange scaling)
		"revPBE": ["GGA_X_PBE_R", "GGA_C_PBE"],
		# RPBE: Revised PBE (functional form of X)
		"RPBE":   ["GGA_X_RPBE", "GGA_C_PBE"],
		# PWP: PW91 exchange + P86 correlation
		"PW91":   ["GGA_X_PW91", "GGA_C_P86"],
		# OLYP: Handy’s ‘optimal’ exchange and Lee-Yang-Parr correlation (the optimized exchange and LYP)
		"OLYP":   ["GGA_X_OPTX", "GGA_C_LYP", "GGA_XC_OPWLYP_D"],
		# OPBE: the optimized exchange and PBE
		"OPBE":   ["GGA_XC_OPBE_D"],
		#XLYP: the Xu/Goddard exchange and LYP
		"XLYP":   ["GGA_XC_XLYP"],
		# B97-D: Grimme’s GGA including D2 dispersion correction
		"B97-D":  ["GGA_XC_B97_D"],
                "B97-D3": ["GGA_XC_B97_D"],
		#***************************************
		# Meta GGA and hybrid meta-GGA functals.
		#***************************************
		# TPSS: the TPPS functional
		"TPSS":   ["LDA_X", "MGGA_X_TPSS", "LDA_C_PW", "MGGA_C_TPSS"],
		# TPSSh: The hybrid version of TPSS (10% HF exchange)
		"TPSSh":  ["HYB_MGGA_XC_TPSSH"],
		# TPSS0: A 25% exchange version of TPSSh that yields improved energetics
                # compared to TPSSh but is otherwise not well tested
		# "TPSS0":  ["HYB_MGGA_XC_TPSSH"],
		# M06L: The Minnesota M06-L meta-GGA functional
		"M06L":   ["MGGA_X_M06_L", "MGGA_C_M06_L"],
		# M06: The M06 hybrid meta-GGA (27% HF exchange)
		"M06":    ["MGGA_X_M06", "MGGA_C_M06"],
		# M062X: The M06-2X version with 54% HF exchange
		"M062X":  ["MGGA_X_M06_2X", "MGGA_C_M06_2X"],
		#***************************************
		# Hybrid functionals
		#***************************************
		# B1LYP: One parameter Hybrid of BLYP
		"B1LYP":     ["HYB_GGA_XC_B1LYP"],
		# B3LYP: Three parameter Hybrid of BLYP
                "B3LYP":     ["HYB_GGA_XC_B3LYP"],
		# B1P: Analogous with Perdew exchange (The one-parameter hybrid version of BP86)
		"B1P":       ["HYB_GGA_XC_B1PW91"],
		# B3P: Analogous with Perdew exchange (The three-parameter hybrid version of BP86)
		"B3P":       ["HYB_GGA_XC_B3P86"],
		# B3PW: The three-parameter hybrid version of PW91
		"B3PW":      ["HYB_GGA_XC_B3PW91"],
		# G1LYP: 1 par. analog with Gill 96 X
		"G1LYP":     ["GGA_X_G96", "GGA_C_OP_G96"],
		# G3LYP: 3 par. analog with Gill 96 X
		# G1P: similar with P correlation
		# G3P # similar with P correlation
		# PBE0: 1 parameter version of PBE (PBEH, PBE0)
		"PBE0":      ["HYB_GGA_XC_PBEH"],
		# PW1PW: One-parameter hybrid version of PW91
		"PW1PW":     ["HYB_GGA_XC_B1PW91"],
                "PWP91_1":   ["HYB_GGA_XC_B1PW91"],
		"PWP1":      ["HYB_GGA_XC_B1PW91"],
		# mPW1PW: 1 parameter version of mPWPW
		"mPW1PW":    ["HYB_GGA_XC_MPW1PW"],
		# mPW1LYP: One-parameter hybrid version of mPWLYP
		"mPW1LYP":   ["HYB_GGA_XC_MPWLYP1M"],
		# O3LYP: 3 parameter version of OLYP
		"O3LYP":     ["HYB_GGA_XC_O3LYP"],
		# X3LYP: 3 parameter version of XLYP
		"X3LYP":     ["HYB_GGA_XC_X3LYP"],
		# PW6B95: Hybrid functional by Truhlar [93]
		"PW6B95":    ["HYB_MGGA_XC_PW6B95"],
		# B97: Becke’s original hybrid
		"B97":       ["HYB_GGA_XC_B97"],
		# BHANDHLYP: Half-and-half Becke hybrid functional
		"BHANDHLYP": ["HYB_GGA_XC_BHANDHLYP"],
		#***************************************
		# Double Hybrid functionals (mix in MP2)
		#***************************************
		# B2PLYP: Grimme’s 2006 double hybrid
		"B2PLYP":    ["HYB_GGA_XC_B2PLYP"],
		# mPW2PLYP: Schwabe/Grimme improved double hybrid
		"mB2PLYP":   ["HYB_GGA_XC_MB2PLYP"],
		# PWPB95: New Grimme double hybrid [94]
		"PWPB95":    ["HYB_MGGA_XC_MPW1B95"]
            }
        nomadNames = functionalMap.get(functional)
        if not nomadNames:
            raise Exception("Unhandled xc functional %s found" % functional)
        for name in nomadNames:
            s = backend.openSection("section_XC_functionals")
            backend.addValue('XC_functional_name', name)
            backend.closeSection("section_XC_functionals", s)

    def onClose_section_method(self, backend, gIndex, value):
        if value["electronic_structure_method"][-1] == "TDDFT":
            gi = backend.openSection("section_method_to_method_refs")
            backend.addValue("method_to_method_ref", (gIndex - 1))
            backend.addValue("method_to_method_kind", 'starting_point')
            backend.closeSection("section_method_to_method_refs", gi)

            gi = backend.openSection("section_calculation_to_calculation_refs")
            backend.addValue("calculation_to_calculation_ref", (self.CalculationGIndex - 1))
            backend.addValue("calculation_to_calculation_kind", 'starting_point')
            backend.closeSection("section_calculation_to_calculation_refs", gi)

    def onOpen_section_single_configuration_calculation(self, backend, gIndex, value):
        self.CalculationGIndex = gIndex

    def onClose_section_eigenvalues(self, backend, gIndex, value):
        occupations = np.array(value["x_orca_orbital_occupation_nb"])
        eigenvalues = np.array(value["x_orca_orbital_energy"])

        number_of_eigenvalues = len(eigenvalues)
        backend.addValue("number_of_eigenvalues", number_of_eigenvalues)
        backend.addValue("number_of_eigenvalues_kpoints", 1)

        occupations = np.array(value["x_orca_orbital_occupation_nb"])
        eigenvalues = np.array(value["x_orca_orbital_energy"])

        backend.addArrayValues("eigenvalues_occupation", occupations.reshape([1, 1, number_of_eigenvalues]))

        backend.addArrayValues("eigenvalues_values", eigenvalues.reshape([1, 1, number_of_eigenvalues]))

    def onClose_section_excited_states(self, backend, gIndex, value):
        number_of_excited_states = len(value["x_orca_excitation_energy"])
        backend.addValue("number_of_excited_states", number_of_excited_states)

        backend.addArrayValues("excitation_energies", np.array(value["x_orca_excitation_energy"]))

        backend.addArrayValues("oscillator_strengths", np.array(value["x_orca_oscillator_strength"]))

        x = value["x_orca_transition_dipole_moment_x"]
        y = value["x_orca_transition_dipole_moment_y"]
        z = value["x_orca_transition_dipole_moment_z"]
        tdm = np.zeros((number_of_excited_states,3), dtype=float)
        tdm[:, 0] = x
        tdm[:, 1] = y
        tdm[:, 2] = z
        backend.addArrayValues("transition_dipole_moments", tdm)

    def onClose_section_run(self, backend, gIndex, value):
        versionPieces = [value[x] for x in ["x_orca_program_version", "x_orca_program_svn", "x_orca_program_compilation_date", "x_orca_program_compilation_time"]]
        version = []
        for v in versionPieces:
            if v:
                version.append(v[0])
        if version:
            backend.addValue("program_version", " ".join(version))

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
        startReStr = r"\s*\* O  R  C  A \*\s*",
        forwardMatch = True,
        sections = ["section_run"],
        fixedStartValues = {"program_name": "ORCA", "program_basis_set_type": "Gaussians" },
        subMatchers = [
            SM(name = 'ProgramHeader',
               startReStr = r"\s*\* O   ?R   ?C   ?A \*\s*",
               subMatchers = [
                   SM(r"\s*Program Version\s*(?P<x_orca_program_version>[0-9a-zA-Z_.].*)"),
                   SM(r"\s*\(SVN:\s*\$(?P<x_orca_program_svn>[^$]+)\$\)\s*"),
                   SM(r"\s*\(\$Date\:\s*(?P<x_orca_program_compilation_date>[0-9a-zA-Z].+?)\s*\$\)")
                   ]),
            buildSinglePointMatcher(),
            buildGeoOptMatcher(),
            buildMp2Matcher(),
            buildCIMatcher(),
            buildTDDFTMatcher()
            ])

def buildSinglePointMatcher():
#
# b) SimpleMatcher for Single Point Calculation:
# **********************************************
#
  return SM(name = 'SinglePointCalculation',
     startReStr = r"\s*\* Single Point Calculation \*\s",
     sections = ["section_single_configuration_calculation"],
     subMatchers = buildSinglePointSubMatchers()
  )

def buildSinglePointSubMatchers():
    return [
       # Get atomic positions:
       SM(name = 'Atomic Coordinates',
          startReStr = r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*",
          sections = ["section_system", "x_orca_atom_positions"],
          subMatchers = [
          SM(r"\s+(?P<x_orca_atom_labels>[a-zA-Z]+)\s+(?P<x_orca_atom_positions_x__angstrom>[-+0-9.]+)\s+(?P<x_orca_atom_positions_y__angstrom>[-+0-9.]+)\s+(?P<x_orca_atom_positions_z__angstrom>[-+0-9.]+)", repeats = True)
          ]
       ),
       # Get basis set information:
       SM(name = 'Basis set information',
          startReStr = r"BASIS SET INFORMATION\s*",
          sections = ["section_basis_set"],
          subMatchers = [
          # Atom labels and basis set:
          SM(r"\s*Group\s+[0-9]+\s+Type\s+(?P<x_orca_basis_set_atom_labels>[a-zA-Z]+)\s+:\s+(?P<x_orca_basis_set>[0-9a-z]+)\s+contracted\s+to\s+(?P<x_orca_basis_set_contracted>[0-9a-z]+)\s+pattern\s+\{[0-9/]+\}", repeats = True),
          # Auxiliary basis set information:
          SM(name = 'Auxiliary basis set information',
             startReStr = r"AUXILIARY BASIS SET INFORMATION\s*",
             sections = ["section_basis_set"],
             subMatchers = [
             SM(r"\s*Group\s+[0-9]+\s+Type\s+(?P<x_orca_basis_set_atom_labels>[a-zA-Z]+)\s+:\s+(?P<x_orca_auxiliary_basis_set>[0-9a-z]+)\s+contracted\s+to\s+(?P<x_orca_auxiliary_basis_set_contracted>[0-9a-z]+)\s+pattern\s+\{[0-9/]+\}", repeats = True)
             ]
          )
          ]
       ),
       # Basis set statistics and startup info:
       SM(name = 'Basis set statistics and startup info',
          startReStr = r"\s*BASIS SET STATISTICS AND STARTUP INFO\s*",
          sections = ["section_basis_set"],
          subMatchers = [
          # Gaussian basis set:
          SM(name = 'Gaussian basis set',
             startReStr = r"\s*Gaussian basis set:\s*",
             subMatchers = [
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<x_orca_nb_of_primitive_gaussian_shells>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of primitive gaussian functions\s+\.\.\.\s+(?P<x_orca_nb_of_primitive_gaussian_functions>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of contracted shells\s+\.\.\.\s+(?P<x_orca_nb_of_contracted_shells>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of contracted basis functions\s+\.\.\.\s+(?P<x_orca_nb_of_contracted_basis_functions>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+Highest angular momentum\s+\.\.\.\s+(?P<x_orca_highest_angular_moment>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+Maximum contraction depth\s+\.\.\.\s+(?P<x_orca_maximum_contraction_depth>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<x_orca_nb_primitive_gaussian_shells>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             ]
          ),
          # Gaussian auxiliary basis set:
          SM(name = 'Gaussian auxiliary basis set',
             startReStr = r"\s*Auxiliary gaussian basis set:\s*",
             subMatchers = [
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<x_orca_nb_of_primitive_gaussian_shells_aux>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of primitive gaussian functions\s+\.\.\.\s+(?P<x_orca_nb_of_primitive_gaussian_functions_aux>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of contracted shells\s+\.\.\.\s+(?P<x_orca_nb_of_contracted_shells_aux>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of contracted aux-basis functions\s+\.\.\.\s+(?P<x_orca_nb_of_contracted_basis_functions_aux>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+Highest angular momentum\s+\.\.\.\s+(?P<x_orca_highest_angular_moment_aux>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+Maximum contraction depth\s+\.\.\.\s+(?P<x_orca_maximum_contraction_depth_aux>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             SM(r"\s+# of primitive gaussian shells\s+\.\.\.\s+(?P<x_orca_nb_primitive_gaussian_shells_aux>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             ]
          )
          ]
       ),
       # SCF Settings:
       SM(name = 'Orca SCF settings',
          startReStr = r"\s*ORCA SCF\s*",
          sections = ["section_method", "x_orca_section_functionals"],
          fixedStartValues = {
              "electronic_structure_method": "DFT"
          },
          subMatchers = [
          # A - For HF methods:
          SM(r"\s+Ab initio Hamiltonian\s+Method\s+\.\.\.\s+(?P<x_orca_XC_functional_type>[-+0-9a-zA-Z()]+)"),
          # B - For DFT methods:
          SM(r"\s+Density Functional\s+Method\s+\.\.\.\s+(?P<x_orca_XC_functional_type>[a-zA-Z()]+)"),
          SM(r"\s+Exchange Functional\s+Exchange\s+\.\.\.\s+(?P<x_orca_exchange_functional>[a-zA-Z0-9]+)"),
          SM(r"\s+X-Alpha parameter\s+XAlpha\s+\.\.\.\s+(?P<x_orca_xalpha_param>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Becke's b parameter\s+XBeta\s+\.\.\.\s+(?P<x_orca_beckes_beta_param>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Correlation Functional\s+Correlation\s+\.\.\.\s+(?P<x_orca_correl_functional>[a-zA-Z0-9]+)"),
          SM(r"\s+LDA part of GGA corr\.\s+LDAOpt\s+\.\.\.\s+(?P<x_orca_lda_part_of_gga_corr>[a-zA-Z-+0-9]+)"),
          # Common settings:
          SM(r"\s+Scalar relativistic method\s+\.\.\.\s+(?P<x_orca_scalar_relativistic_method>[a-zA-Z-+0-9]+)"),
          SM(r"\s+Speed of light used\s+Velit\s+\.\.\.\s+(?P<x_orca_speed_of_light_used>[0-9.]+)"),
          SM(r"\s+Hartree-Fock type\s+HFTyps+\.\.\.\s+(?P<x_orca_hf_type>[a-zA-Z]+)"),
          SM(r"\s+Total Charge\s+Charge\s+\.\.\.\s+(?P<x_orca_total_charge>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Multiplicity\s+Mult\s+\.\.\.\s+(?P<x_orca_multiplicity>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Number of Electrons\s+NEL\s+\.\.\.\s+(?P<x_orca_nelectrons>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Nuclear Repulsion\s+ENuc\s+\.\.\.\s+(?P<x_orca_nuclear_repulsion__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          # Convergence Tolerance:
          SM(r"\s+Convergence Check Mode ConvCheckMode\s*\.\.\.\s+(?P<x_orca_convergence_check_mode>[a-zA-Z-+0-9.eEdD_-]+)"),
          SM(r"\s+Energy Change\s+TolE\s*\.\.\.\s+(?P<x_orca_energy_change_tolerance__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+1-El\. energy change\s*\.\.\.\s+(?P<x_orca_1_elect_energy_change__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          # DFT Grid generation:
          SM(r"\s+General Integration Accuracy\s+IntAcc\s*\.\.\.\s+(?P<x_orca_gral_integ_accuracy>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Radial Grid Type\s+RadialGrid\s*\.\.\.\s+(?P<x_orca_radial_grid_type>[a-zA-Z-_]+)"),
          SM(r"\s+Angular Grid \(max\. acc\.\)\s+AngularGrid\s*\.\.\.\s+(?P<x_orca_angular_grid>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Angular grid pruning method\s+GridPruning\s*\.\.\.\s+(?P<x_orca_grid_pruning_method>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Weight generation scheme\s*WeightScheme\s*\.\.\.\s+(?P<x_orca_weight_gener_scheme>[a-zA-Z0-9]+)"),
          SM(r"\s+Basis function cutoff\s+BFCut\s*\.\.\.\s+(?P<x_orca_basis_fn_cutoff>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Integration weight cutoff\s+WCut\s*\.\.\.\s+(?P<x_orca_integr_weight_cutoff>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+# of grid points \(after initial pruning\)\s+WCut\s*\.\.\.\s+(?P<x_orca_nb_grid_pts_after_initial_pruning>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+# of grid points \(after weights\+screening\)\s*\.\.\.\s+(?P<x_orca_nb_grid_pts_after_weights_screening>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Total number of grid points\s*\.\.\.\s+(?P<x_orca_total_nb_grid_pts>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Total number of batches\s*\.\.\.\s+(?P<x_orca_total_nb_batches>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Average number of points per batch\s*\.\.\.\s+(?P<x_orca_avg_nb_points_per_batch>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Average number of grid points per atom\s*\.\.\.\s+(?P<x_orca_avg_nb_grid_pts_per_atom>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)")
          ]
       ),
       # SCF iterations:
       SM(name = 'Orca SCF iterations',
          startReStr = r"\s*SCF ITERATIONS\s*",
          subMatchers = [
          SM(r"\s*(?P<x_orca_iteration_nb>[0-9]+)\s+(?P<energy_total_scf_iteration__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)",
             sections = ["section_scf_iteration"],
             repeats = True),
          ]
       ),
#      *****************************************************
#      *                     SUCCESS                       *
#      *           SCF CONVERGED AFTER  XY CYCLES          *
#      *****************************************************
       SM(name = 'Final step after convergence',
          startReStr = r"Setting up the final grid:",
          sections = ["section_scf_iteration", "section_basis_set"],
          subMatchers = [
          # Final DFT Grid generation:
          SM(r"\s+General Integration Accuracy\s+IntAcc\s*\.\.\.\s+(?P<x_orca_gral_integ_accuracy_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Radial Grid Type\s+RadialGrid\s*\.\.\.\s+(?P<x_orca_radial_grid_type_final>[a-zA-Z-_]+)"),
          SM(r"\s+Angular Grid \(max\. acc\.\)\s+AngularGrid\s*\.\.\.\s+(?P<x_orca_angular_grid_final>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Angular grid pruning method\s+GridPruning\s*\.\.\.\s+(?P<x_orca_grid_pruning_method_final>[a-zA-Z-+0-9.eEdD-_]+)"),
          SM(r"\s+Weight generation scheme\s*WeightScheme\s*\.\.\.\s+(?P<x_orca_weight_gener_scheme_final>[a-zA-Z0-9]+)"),
          SM(r"\s+Basis function cutoff\s+BFCut\s*\.\.\.\s+(?P<x_orca_basis_fn_cutoff_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Integration weight cutoff\s+WCut\s*\.\.\.\s+(?P<x_orca_integr_weight_cutoff_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+# of grid points \(after initial pruning\)\s+WCut\s*\.\.\.\s+(?P<x_orca_nb_grid_pts_after_initial_pruning_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+# of grid points \(after weights\+screening\)\s*\.\.\.\s+(?P<x_orca_nb_grid_pts_after_weights_screening_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Total number of grid points\s*\.\.\.\s+(?P<x_orca_total_nb_grid_pts_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Total number of batches\s*\.\.\.\s+(?P<x_orca_total_nb_batches_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Average number of points per batch\s*\.\.\.\s+(?P<x_orca_avg_nb_points_per_batch_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s+Average number of grid points per atom\s*\.\.\.\s+(?P<x_orca_avg_nb_grid_pts_per_atom_final>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          # Final SCF total Energy:
          SM(name = 'Total Energy',
             startReStr = r"\s*TOTAL SCF ENERGY\s*",
             sections = [],
             subMatchers = [
             SM(r"\s*Total Energy\s+:\s+(?P<energy_total__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             # Energy Components:
             SM(name = 'Energy Components',
                startReStr = r"\s*Components:\s*",
                subMatchers = [
                SM(r"\s*Nuclear Repulsion\s*:\s+(?P<x_orca_nuc_repulsion__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*Electronic Energy\s*:\s+(?P<x_orca_elec_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*One Electron Energy:\s+(?P<x_orca_one_elec_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*Two Electron Energy:\s+(?P<x_orca_two_elec_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                # Virial Components:
                SM(r"\s*Potential Energy\s*:\s+(?P<x_orca_potential_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*Kinetic Energy\s*:\s+(?P<x_orca_kinetc_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*Virial Ratio\s*:\s+(?P<x_orca_virial_ratio>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                # DFT Components:
                SM(r"\s*N\(Alpha\)\s*:\s+(?P<x_orca_nb_elect_alpha_channel>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*N\(Beta\)\s*:\s+(?P<x_orca_nb_elect_beta_channel>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*N\(Total\)\s*:\s+(?P<x_orca_nb_elect_total>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*E\(X\)\s*:\s+(?P<x_orca_exchange_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*E\(C\)\s*:\s+(?P<x_orca_correlation_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
                SM(r"\s*E\(XC\)\s*:\s+(?P<x_orca_exchange_correlation_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)")
                ]
             )
             ]
          ),
          # Final SCF convergence:
          SM(r"\s*Last Energy change\s+\.\.\.\s+(?P<x_orca_last_energy_change__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+Tolerance :\s*(?P<x_orca_last_energy_change_tolerance__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Last MAX-Density change\s+\.\.\.\s+(?P<x_orca_last_max_density_change__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+Tolerance :\s*(?P<x_orca_last_max_density_tolerance__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Last RMS-Density change\s+\.\.\.\s+(?P<x_orca_last_rms_density_change__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+Tolerance :\s*(?P<x_orca_last_rms_density_tolerance__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)")
          ]
       ),
       # Orbitals Energies and Mulliken population analysis:
       SM(name = 'Orbital Energies',
          startReStr = r"\s*ORBITAL ENERGIES\s*",
          sections = ["section_dos", "section_eigenvalues"],
          subMatchers = [
          SM(r"\s*(?P<x_orca_orbital_nb>[0-9]+)\s+(?P<x_orca_orbital_occupation_nb>[-+.0-9]+)\s+(?P<x_orca_orbital_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)", repeats = True),
          # Mulliken population analysis:
          SM(name = 'Mulliken population analysis',
             startReStr = r"\s*\* MULLIKEN POPULATION ANALYSIS \*\s*",
             sections = [],
             subMatchers = [
             SM(r"\s*(?P<x_orca_atom_nb>[0-9]+)\s+(?P<x_orca_atom_species>[a-zA-Z]+):\s+(?P<x_orca_mulliken_atom_charge>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)", repeats = True),
             SM(r"\s*Sum of atomic charges:\s*(?P<x_orca_mulliken_total_charge>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
             # Mulliken reduced orbital charges (mroc):
             SM(r"\s*(?P<x_orca_atom_nb_mroc>[0-9]+)\s+(?P<x_orca_atom_species_mroc>[a-zA-Z]+)(?P<x_orca_atom_orbital_mroc>[-+0-9a-zA-Z]+)\s*:\s*(?P<x_orca_mulliken_partial_orbital_charge_mroc>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)", repeats = True)
             ]
          )
          ]
       ),
       # Time table:
       SM(name = 'timings',
          startReStr = r"\s*TIMINGS\s*",
          subMatchers = [
          SM(r"\s*Total SCF time:\s+(?P<x_orca_total_days_time>[0-9]+) days (?P<x_orca_total_hours_time>[0-9]+) hours (?P<x_orca_total_mins_time>[0-9]+) min (?P<x_orca_total_secs_time>[0-9]+) sec"),
          SM(r"\s*Total time\s*\.\.\.\.\s*(?P<x_orca_final_time>[0-9.]+) sec"),
          SM(r"\s*Sum of individual times\s*\.\.\.\.\s*(?P<x_orca_sum_individual_times>[0-9.]+) sec"),
          SM(r"\s*Fock matrix formation\s*\.\.\.\.\s*(?P<x_orca_fock_matrix_formation>[0-9.]+) sec"),
          SM(r"\s*Coulomb formation\s*\.\.\.\.\s*(?P<x_orca_coulomb_formation>[0-9.]+) sec"),
          SM(r"\s*Split-RI-J\s*\.\.\.\.\s*(?P<x_orca_split_rj>[0-9.]+) sec"),
          SM(r"\s*XC integration\s*\.\.\.\.\s*(?P<x_orca_xc_integration>[0-9.]+) sec"),
          SM(r"\s*Basis function eval\.\s*\.\.\.\.\s*(?P<x_orca_basis_fn_evaluation>[0-9.]+) sec"),
          SM(r"\s*Density eval\.\s*\.\.\.\.\s*(?P<x_orca_density_evaluation>[0-9.]+) sec"),
          SM(r"\s*XC-Functional eval\.\s*\.\.\.\.\s*(?P<x_orca_xc_functional_evaluation>[0-9.]+) sec"),
          SM(r"\s*XC-Potential eval\.\s*\.\.\.\.\s*(?P<x_orca_potential_evaluation>[0-9.]+) sec"),
          SM(r"\s*Diagonalization\s*\.\.\.\.\s*(?P<x_orca_diagonalization>[0-9.]+) sec"),
          SM(r"\s*Density matrix formation\s*\.\.\.\.\s*(?P<x_orca_density_matrix_formation>[0-9.]+) sec"),
          SM(r"\s*Population analysis\s*\.\.\.\.\s*(?P<x_orca_population_analysis>[0-9.]+) sec"),
          SM(r"\s*Initial guess\s*\.\.\.\.\s*(?P<x_orca_initial_guess>[0-9.]+) sec"),
          SM(r"\s*Orbital Transformation\s*\.\.\.\.\s*(?P<x_orca_orbital_transformation>[0-9.]+) sec"),
          SM(r"\s*Orbital Orthonormalization\s*\.\.\.\.\s*(?P<x_orca_orbital_orthonormalization>[0-9.]+) sec"),
          SM(r"\s*DIIS solution\s*\.\.\.\.\s*(?P<x_orca_diis_solution>[0-9.]+) sec"),
          SM(r"\s*Grid generation\s*\.\.\.\.\s*(?P<x_orca_grid_generation>[0-9.]+) sec")
          ]
       )]

def buildGeoOptMatcher():
#
# c) SimpleMatcher for geometry optimization:
# *******************************************
#
       return SM(name = 'Geometry optimization',
          startReStr = r"\s*\* Geometry Optimization Run \*\s*",
          sections = ["section_sampling_method"],
          subMatchers = [
          # Geometry optimization settings:
          SM(r"Update method\s*(?P<x_orca_update_method>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_update_method_name>[a-zA-Z]+)"),
          SM(r"Choice of coordinates\s*(?P<x_orca_coords_choice>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_coords_choice_name>[a-zA-Z ]+)"),
          SM(r"Initial Hessian\s*(?P<x_orca_initial_hessian>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_initial_hessian_name>[a-zA-Z'\( \)]+)"),
          SM(r"Energy Change\s*(?P<x_orca_energy_change_tol>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_energy_change_tol_value__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) Eh"),
          SM(r"Max\. Gradient\s*(?P<x_orca_max_gradient_tol>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_max_gradient_tol_value__hartree_bohr_1>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) Eh/bohr"),
          SM(r"RMS Gradient\s*(?P<x_orca_rms_gradient_tol>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_rms_gradient_tol_value__hartree_bohr_1>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) Eh/bohr"),
          SM(r"Max\. Displacement\s*(?P<x_orca_max_displacement_tol>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_max_displacement_tol_value__bohr>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) bohr"),
          SM(r"RMS Displacement\s*(?P<x_orca_rms_displacement_tol>[a-zA-Z]+)\s+\.\.\.\.\s+(?P<x_orca_rms_displacement_tol_value__bohr>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) bohr"),
          SM(name="geoOptCycle",
             startReStr=r"\s*\*\s*GEOMETRY OPTIMIZATION CYCLE\s*(?P<x_orca_geo_opt_cycle>[0-9]+)\s*\*\s*",
             sections = ["section_single_configuration_calculation"],
             subMatchers = buildSinglePointSubMatchers()
          ),
          # Final geometry:
          # CARTESIAN COORDINATES (ANGSTROEM)
          SM(name = 'final geometry',
             startReStr = r"\s*\*\*\*\s*THE OPTIMIZATION HAS CONVERGED\s*\*\*\*\s*",
             sections = ["x_orca_final_geometry"],
             subMatchers = [
             SM(r"\s+(?P<x_orca_atom_labels_geo_opt>[a-zA-Z]+)\s+(?P<x_orca_atom_positions_x_geo_opt__angstrom>[-+0-9.]+)\s+(?P<x_orca_atom_positions_y_geo_opt__angstrom>[-+0-9.]+)\s+(?P<x_orca_atom_positions_z_geo_opt__angstrom>[-+0-9.]+)", repeats = True)
             ]
          ),
       # *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***
       #
          SM(name = 'Final step after geometry convergence',
             startReStr = r"Setting up the final grid:",
             sections = ["section_single_configuration_calculation"],
             subMatchers = buildSinglePointSubMatchers()
         )
        ]
       )

def buildMp2Matcher():
#
# d) Post-processing calculations:
# ********************************
#
       # MP2 Calculation (post-proc):
    return SM(name = 'mp2',
          startReStr = r"\s*ORCA MP2 CALCULATION\s*",
              sections = ["section_method", "section_single_configuration_calculation"],
          fixedStartValues = {
              "electronic_structure_method": "MP2"
          },
          subMatchers = [
          SM(r"\s*Dimension of the basis\s*\.\.\.\s*(?P<x_orca_mp2_basis_dimension>[0-9.]+)"),
          SM(r"\s*Overall scaling of the MP2 energy\s*\.\.\.\s*(?P<x_orca_scaling_mp2_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) Eh"),
          SM(r"\s*Dimension of the aux-basis\s*\.\.\.\s*(?P<x_orca_mp2_aux_basis_dimension>[0-9.]+)"),
          SM(r"\s*RI-MP2 CORRELATION ENERGY:\s*(?P<energy_method_current__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) Eh"),
          SM(r"\s*MP2 TOTAL ENERGY:\s*(?P<energy_total__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) Eh")
          ]
       )
def buildCIMatcher():
       # Driven CI (post-proc):
    return SM(name = 'CI',
          startReStr = r"\s*ORCA-MATRIX DRIVEN CI\s*",
          sections = ["section_method", "section_single_configuration_calculation"],
#          fixedStartValues = {
#              "electronic_structure_method": "CI"
#          },
          subMatchers = [
          SM(r"\s*Correlation treatment\s*\.\.\.\s*(?P<electronic_structure_method>[()a-zA-Z0-9.]+)"),
          SM(r"\s*Single excitations\s*\.\.\.\s*(?P<x_orca_single_excitations_on_off>[a-zA-Z]+)"),
          SM(r"\s*Orbital optimization\s*\.\.\.\s*(?P<x_orca_orbital_opt_on_off>[a-zA-Z]+)"),
          SM(r"\s*Calculation of Z vector\s*\.\.\.\s*(?P<x_orca_z_vector_calc_on_off>[a-zA-Z]+)"),
          SM(r"\s*Calculation of Brueckner orbitals\s*\.\.\.\s*(?P<x_orca_Brueckner_orbitals_calc_on_off>[a-zA-Z]+)"),
          SM(r"\s*Perturbative triple excitations\s*\.\.\.\s*(?P<x_orca_perturbative_triple_excitations_on_off>[a-zA-Z]+)"),
          SM(r"\s*Calculation of F12 correction\s*\.\.\.\s*(?P<x_orca_f12_correction_on_off>[a-zA-Z]+)"),
          SM(r"\s*Frozen core treatment\s*\.\.\.\s*(?P<x_orca_frozen_core_treatment>[0-9.a-zA-Z( )]+)"),
          SM(r"\s*Reference Wavefunction\s*\.\.\.\s*(?P<x_orca_reference_wave_function>[0-9.a-zA-Z]+)"),
          SM(r"\s*Number of AO's\s*\.\.\.\s*(?P<x_orca_nb_of_atomic_orbitals>[0-9]+)"),
          SM(r"\s*Number of electrons\s*\.\.\.\s*(?P<x_orca_nb_of_electrons>[0-9]+)"),
          SM(r"\s*Number of correlated electrons\s*\.\.\.\s*(?P<x_orca_nb_of_correlated_electrons>[0-9]+)"),
          SM(r"\s*Integral transformation\s*\.\.\.\s*(?P<x_orca_integral_transformation>[a-zA-Z( )]+)"),
          SM(r"\s*K\(C\) Formation\s*\.\.\.\s*(?P<x_orca_K_C_formation>[+-a-zA-Z( )]+)"),
          SM(r"\s*Convergence tolerance \(max\. residuum\)\s*\.\.\.\s*(?P<x_orca_convergence_tol_max_residuum>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Level shift for amplitude update\s*\.\.\.\s*(?P<x_orca_level_shift_amplitude_update>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          # Partial Coulomb Transformation:
          SM(r"\s*Transformation type\s*\.\.\.\s*(?P<x_orca_coulomb_transformation_type>[a-zA-Z ()]+)"),
          SM(r"\s*Dimension of the basis\s*\.\.\.\s*(?P<x_orca_coulomb_transformation_dimension_basis>[0-9]+)"),
          SM(r"\s*Number of internal alpha-MOs\s*\.\.\.\s*(?P<x_orca_nb_internal_alpha_mol_orbitals>[-+0-9( )]+)"),
          SM(r"\s*Number of internal beta-MOs\s*\.\.\.\s*(?P<x_orca_nb_internal_beta_mol_orbitals>[-+0-9( )]+)"),
          SM(r"\s*Pair cutoff\s*\.\.\.\s*(?P<x_orca_pair_cutoff__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?) Eh"),
          SM(r"\s*AO-integral source\s*\.\.\.\s*(?P<x_orca_atomic_orbital_integral_source>[a-zA-Z]+)"),
          SM(r"\s*Integral package used\s*\.\.\.\s*(?P<x_orca_integral_package_used>[a-zA-Z]+)"),
          SM(r"\s*Number of Alpha-MO pairs included\s*\.\.\.\s*(?P<x_orca_nb_alpha_pairs_included>[0-9]+)"),
          SM(r"\s*Number of Beta-MO pairs included\s*\.\.\.\s*(?P<x_orca_nb_beta_pairs_included>[0-9]+)"),
          # Spin-unrestricted guess:
          SM(r"\s*EMP2\(aa\)=\s*(?P<x_orca_mp2_energy_spin_aa__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*EMP2\(bb\)=\s*(?P<x_orca_mp2_energy_spin_bb__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*EMP2\(ab\)=\s*(?P<x_orca_mp2_energy_spin_ab__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*E\(0\)\s*\.\.\.\s*(?P<x_orca_mp2_initial_guess__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*E\(MP2\)\s*\.\.\.\s*(?P<x_orca_mp2_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Initial E\(tot\)\s*\.\.\.\s*(?P<x_orca_mp2_total_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*<T\|T>\s*\.\.\.\s*(?P<x_orca_T_and_T_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Number of pairs included\s*\.\.\.\s*(?P<x_orca_total_nb_pairs_included>[0-9]+)"),
          # iterations (e.g.:UHF COUPLED CLUSTER ITERATIONS):
          SM(r"\s*(?P<x_orca_ci_iteration_nb>[0-9]+)\s+(?P<x_orca_ci_total_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+(?P<x_orca_ci_correl_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+(?P<x_orca_ci_deltaE_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+(?P<x_orca_ci_residual_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+(?P<x_orca_ci_iteration_time>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s+(?P<x_orca_ci_half_s_and_s_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)", repeats = True),
          # Final Coupled Cluster Energies:
          SM(r"\s*E\(CORR\)\s*\.\.\.\s*(?P<x_orca_ccsd_correlation_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*E\(TOT\)\s*\.\.\.\s*(?P<x_orca_ccsd_total_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Singles norm <S\|S>\*\*1/2\s*\.\.\.\s*(?P<x_orca_single_norm_half_ss__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*T1 diagnostic\s*\.\.\.\s*(?P<x_orca_t1_diagnostic__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Triples Correction (T)\s*\.\.\.\s*(?P<x_orca_ccsdt_total_triples_correction__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*alpha-alpha-alpha\s*\.\.\.\s*(?P<x_orca_ccsdt_aaa_triples_contribution__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*alpha-alpha-beta\s*\.\.\.\s*(?P<x_orca_ccsdt_aab_triples_contribution__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*alpha-beta -beta\s*\.\.\.\s*(?P<x_orca_ccsdt_aba_triples_contribution__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*beta -beta -beta\s*\.\.\.\s*(?P<x_orca_ccsdt_bbb_triples_contribution__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*Final correlation energy\s*\.\.\.\s*(?P<x_orca_ccsdt_final_corr_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*E\(CCSD\)\s*\.\.\.\s*(?P<x_orca_ccsd_final_energy__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"),
          SM(r"\s*E\(CCSD\(T\)\)\s*\.\.\.\s*(?P<energy_total__hartree>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)")
          ]
       )

def buildTDDFTMatcher():
    # TDDFT Calculation (post-proc):
    return SM(name='tddft',
              startReStr=r"\s*ORCA TD-DFT/TDA CALCULATION\s*",
              sections=["section_single_configuration_calculation", "section_method", "section_excited_states"],
              fixedStartValues={
                  "electronic_structure_method": "TDDFT"
              },
              subMatchers=[
              SM(r"\s*ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\s*"),
              SM(r"\s*[0-9]+\s+(?P<x_orca_excitation_energy__inversecm>[-+.0-9]+)\s+[-+.0-9]+\s+(?P<x_orca_oscillator_strength>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)"
                 r"\s+[-+.0-9]+\s+(?P<x_orca_transition_dipole_moment_x__bohr>[-+.0-9]+)\s+(?P<x_orca_transition_dipole_moment_y__bohr>[-+.0-9]+)"
                 r"\s+(?P<x_orca_transition_dipole_moment_z__bohr>[-+.0-9]+)\s*", repeats = True),
              SM(r"\s*ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\s*")
              ]
        )

       # Here new stuff:
#       SM(name = '',
#          startReStr = r"",
#          sections = [],
#          subMatchers = [
#          SM(),
#          SM(),
#          SM()
#          ]
#       ),
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

class OrcaParser():
   """ A proper class envolop for running this parser from within python. """
   def __init__(self, backend, **kwargs):
       self.backend_factory = backend

   def parse(self, mainfile):
       from unittest.mock import patch
       logging.info('orca parser started')
       logging.getLogger('nomadcore').setLevel(logging.WARNING)
       backend = self.backend_factory("orca.nomadmetainfo.json")
       with patch.object(sys, 'argv', ['<exe>', '--uri', 'nmd://uri', mainfile]):
           mainFunction(
               mainFileDescription,
               None,
               parserInfo,
               superContext=OrcaContext(),
               superBackend=backend)

       return backend
