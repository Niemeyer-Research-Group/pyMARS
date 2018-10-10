""" Tests the reduction methods implemented by pyMARS """

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import cantera as ct

import drgep
import sensitivity_analysis


def testDRGEPSA():

	# Original model
	path_to_original = os.path.dirname(os.path.abspath(__file__)) + "/../../example_files/gri30.cti" 
	solution_object = ct.Solution(path_to_original)

	# Conditions for reduction	
	conditions = os.path.dirname(os.path.abspath(__file__)) + "/../../example_files/example_input_file.txt"
	error = 5.0
	target_species = ["CH4","O2"]
	retained_species = ["CH4","O2","N2","H2O","CO2"]
	final_error = [0]
	epsilon_star = .5

	# Run DRGEP
	reduced_model = drgep.run_drgep(solution_object, conditions, error, target_species, retained_species, path_to_original, final_error)

	# Expected answer	
	path_to_answer = os.path.dirname(os.path.abspath(__file__)) + "/drgep_gri30.cti"
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same
	assert len(reduced_model.species()) == len(expected_model.species())
	assert len(reduced_model.reactions()) == len(expected_model.reactions())

	# Run SA	
	reduced_model = sensitivity_analysis.run_sa(solution_object, reduced_model, epsilon_star, final_error, conditions, target_species, retained_species, error)

	# Get expected model	
	path_to_answer = os.path.dirname(os.path.abspath(__file__)) + "/sa_gri30.cti"
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same	
	assert len(reduced_model.species()) == len(expected_model.species())
	assert len(reduced_model.reactions()) == len(expected_model.reactions())

