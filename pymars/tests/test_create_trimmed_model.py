""" Tests the create trimmed model unit used by pyMARS """

import os
import pkg_resources

import pytest
import cantera as ct

from ..create_trimmed_model import trim

def relative_location(file):
	file_path = os.path.join(file)
	return pkg_resources.resource_filename(__name__, file_path)


def test_GRI_minus_three():
	"""Tests removal of three species from GRI Mech 3.0.
	"""

	# Original model to remove things from
	solution_object = ct.Solution('gri30.cti')

	# Create exclusion list for test case
	exclusion_list = ["CH4", "O2", "N2"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "gri30.cti")

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 50
	expected_reactions_num = 237

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 3

	# Make sure removed species are not included
	assert "CH4" not in reduced_model.species_names
	assert "O2" not in reduced_model.species_names
	assert "N2" not in reduced_model.species_names

def test_GRI_minus_zero():
	"""Tests removal of zero species from GRI Mech 3.0.
	"""

	# Original model to remove things from
	solution_object = ct.Solution("gri30.cti")

	# Create exclusion list for test case
	exclusion_list = []

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "gri30.cti")

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 53
	expected_reactions_num = 325

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species())

	# Make sure target is not removed 
	assert "CH4" in reduced_model.species_names


def test_artificial_minus_one():
	"""Test removing one species from artificial model.
	"""

	# Original model to remove things from
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["H"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 3
	expected_reactions_num = 1

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 1

	# Make sure removed species are not included
	assert "H" not in reduced_model.species_names


@pytest.mark.xfail
def testArtRemoveAll():
	"""Test removing all four species in an artificial model.

	Fails because Cantera will not produce a Solution with no species/reactions.
	"""

	# Original model to remove things from
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["H", "H2", "O2", "H2O"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 0
	expected_reactions_num = 0

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 4

	# Make sure removed species are not included
	assert "H" not in reduced_model.species_names

def testArtRemoveInvalid():
	"""Test removing species not present in model.
	"""

	# Original model to remove things from
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["CH4"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 4
	expected_reactions_num = 2

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species())

def testArtRemoveInvalidAnd1():
	"""Test removing mixture of species both in and not in artificial model.
	"""

	# Original model to remove things from
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["H", "CH4"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 3
	expected_reactions_num = 1

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 1

	# Make sure removed species are not included
	assert "H" not in reduced_model.species_names

def test_GRI_minus_10():
	"""Test removing 10 species from GRI Mech 3.0
	"""

	# Original model to remove things from
	solution_object = ct.Solution("gri30.cti")

	# Create exclusion list for test case
	exclusion_list = ["CH4", "O2", "N2", "H", "OH", "H2O", "CH2", "CH3", "CO", "AR"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "gri30.cti")

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 43
	expected_reactions_num = 14

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 10

	# Make sure removed species are not included
	for sp in exclusion_list:
		assert sp not in reduced_model.species_names
