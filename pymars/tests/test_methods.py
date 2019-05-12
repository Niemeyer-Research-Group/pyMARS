""" Tests the reduction methods implemented by pyMARS """

import sys
import os
import pkg_resources

import cantera as ct

from .. import drgep
from .. import pfa
from .. import drg
from .. import sensitivity_analysis

def relative_location(file):
	file_path = os.path.join(file)
	return pkg_resources.resource_filename(__name__, file_path)


def testDRGEPSA():

	# Original model
	solution_object = ct.Solution("gri30.cti")

	# Conditions for reduction
	conditions =  relative_location("example_input_file.txt")
	error = 5.0
	target_species = ["CH4","O2"]
	retained_species = ["CH4","O2","N2","H2O","CO2"]
	final_error = [0]
	epsilon_star = .5

	# Run DRGEP
	result = drgep.run_drgep(solution_object, conditions, error, target_species, 
						 	 retained_species, "gri30.cti", final_error, epsilon_star
							 )
	reduced_model = result[0]
	limbo = result[1]

	# Expected answer	
	path_to_answer = relative_location("drgep_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same
	assert len(reduced_model.species()) == len(expected_model.species())
	assert len(reduced_model.reactions()) == len(expected_model.reactions())

	# Run SA	
	reduced_model = sensitivity_analysis.run_sa(solution_object, reduced_model, final_error, 
												conditions, target_species, retained_species, 
												error, limbo
												)

	# Get expected model	
	path_to_answer = relative_location("sa_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same	
	assert reduced_model.species_names == expected_model.species_names
	assert len(reduced_model.reactions()) == len(expected_model.reactions())


def testPFA():

	# Original model
	solution_object = ct.Solution("gri30.cti")

	# Conditions for reduction	
	conditions =  relative_location("example_input_file.txt")
	error = 5.0
	target_species = ["CH4","O2"]
	retained_species = ["CH4","O2","N2","H2O","CO2"]
	final_error = [0]
	epsilon_star = .5

	# Run PFA
	reduced_model = pfa.run_pfa(solution_object, conditions, error, target_species, 
								retained_species, "gri30.cti", final_error
								)[0]

	# Expected answer	
	path_to_answer = relative_location("pfa_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same
	assert reduced_model.species_names == expected_model.species_names
	assert len(reduced_model.reactions()) == len(expected_model.reactions())

def testDRGASA():

	# Original model
	solution_object = ct.Solution("gri30.cti")

	# Conditions for reduction	
	conditions =  relative_location("example_input_file.txt")
	error = 5.0
	target_species = ["CH4","O2"]
	retained_species = ["CH4","O2","N2","H2O","CO2"]
	final_error = [0]
	epsilon_star = .5

	# Run DRG
	result = drg.run_drg(solution_object, conditions, error, target_species, 
						 retained_species, "gri30.cti", final_error, epsilon_star
						 )
	reduced_model = result[0]
	limbo = result[1]

	# Expected answer	
	path_to_answer = relative_location("drg_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same
	assert reduced_model.species_names == expected_model.species_names
	assert len(reduced_model.reactions()) == len(expected_model.reactions())
	
	# Run SA	
	reduced_model = sensitivity_analysis.run_sa(solution_object, reduced_model, final_error, 
												conditions, target_species, retained_species, 
												error, limbo
												)

	# Get expected model	
	path_to_answer = relative_location("drgasa_gri30.cti")
	expected_model = ct.Solution(path_to_answer)

	# Make sure models are the same	
	assert reduced_model.species_names == expected_model.species_names
	assert len(reduced_model.reactions()) == len(expected_model.reactions())
