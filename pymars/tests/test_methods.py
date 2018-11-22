""" Tests the reduction methods implemented by pyMARS """

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import cantera as ct

import drgep
import pfa
import drg
import sensitivity_analysis

ROOT_DIR =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def relative_location( file):
	file_path = os.path.join(ROOT_DIR, file)
	return file_path

def testDRGEPSA():

	# Original model
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	# Conditions for reduction
	conditions =  relative_location("example_files/example_input_file.txt")
	error = 5.0
	target_species = ["CH4","O2"]
	retained_species = ["CH4","O2","N2","H2O","CO2"]
	final_error = [0]
	epsilon_star = .5

	# Run DRGEP
	reduced_model = drgep.run_drgep(solution_object, conditions, error, target_species, retained_species, path_to_original, final_error)

	# Expected answer	
	path_to_answer = relative_location("pymars/tests/drgep_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same
	assert len(reduced_model.species()) == len(expected_model.species())
	assert len(reduced_model.reactions()) == len(expected_model.reactions())

	# Run SA	
	reduced_model = sensitivity_analysis.run_sa(solution_object, reduced_model, epsilon_star, final_error, conditions, target_species, retained_species, error)

	# Get expected model	
	path_to_answer = relative_location("pymars/tests/sa_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same	
	assert reduced_model.species_names == expected_model.species_names
	assert len(reduced_model.reactions()) == len(expected_model.reactions())


def testPFA():

	# Original model
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	# Conditions for reduction	
	conditions =  relative_location("example_files/example_input_file.txt")
	error = 5.0
	target_species = ["CH4","O2"]
	retained_species = ["CH4","O2","N2","H2O","CO2"]
	final_error = [0]
	epsilon_star = .5

	# Run PFA
	reduced_model = pfa.run_pfa(solution_object, conditions, error, target_species, retained_species, path_to_original, final_error)

	# Expected answer	
	path_to_answer = relative_location("pymars/tests/pfa_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same
	assert reduced_model.species_names == expected_model.species_names
	assert len(reduced_model.reactions()) == len(expected_model.reactions())

def testDRG():

	# Original model
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	# Conditions for reduction	
	conditions =  relative_location("example_files/example_input_file.txt")
	error = 5.0
	target_species = ["CH4","O2"]
	retained_species = ["CH4","O2","N2","H2O","CO2"]
	final_error = [0]
	epsilon_star = .5

	# Run DRG
	reduced_model = drg.run_drg(solution_object, conditions, error, target_species, retained_species, path_to_original, final_error)

	# Expected answer	
	path_to_answer = relative_location("pymars/tests/drg_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same
	assert reduced_model.species_names == expected_model.species_names
	assert len(reduced_model.reactions()) == len(expected_model.reactions())
