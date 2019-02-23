""" Tests the create trimmed model unit used by pyMARS """

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import pytest
import cantera as ct

from create_trimmed_model import trim

ROOT_DIR =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def relative_location(file):
	file_path = os.path.join(ROOT_DIR, file)
	return file_path

########
# Input: GRI Mech 3.0 (53 species), removal list of three species
# Output: Reduced mechanism with 50 species and none of the species in the removal list
########
def testGRIMinus3():

	# Original model to remove things from
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["CH4", "O2", "N2"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "gri30.cti")[1]

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 50
	expected_reactions_num = 237

	# Print stuff
	print("-- Trim 3 from GRI3.0 --")
	print("Trimmed num species: " + str(reduced_model_num_species))
	print("Trimmed num reactions: " + str(reduced_model_num_reactions))
	print("Trimmed model species: ")
	print(reduced_model.species_names)

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 3

	# Make sure removed species are not included
	assert "CH4" not in reduced_model.species_names
	assert "O2" not in reduced_model.species_names
	assert "N2" not in reduced_model.species_names

########
# Input: GRI Mech 3.0 (53 species), empty removal list
# Output: Reduced mechanism with 53 species, identical to the input 
########
def testGRIMinus0():

	# Original model to remove things from
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = []

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "gri30.cti")[1]

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 53
	expected_reactions_num = 325

	# Print stuff
	print("-- Trim 0 from GRI3.0 --")
	print("Trimmed num species: " + str(reduced_model_num_species))
	print("Trimmed num reactions: " + str(reduced_model_num_reactions))
	print("Trimmed model species: ")
	print(reduced_model.species_names)

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species())

	# Make sure target is not removed 
	assert "CH4" in reduced_model.species_names

########
# Input: Artificial mechanism (4 species), removal list of one species
# Output: Reduced mechanism with 3 species and not the species in the removal list
########
def testArtMinus1():

	# Original model to remove things from
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["H"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")[1]

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 3
	expected_reactions_num = 1

	# Print stuff
	print("-- Trim 1 from Artificial Model --")
	print("Trimmed num species: " + str(reduced_model_num_species))
	print("Trimmed num reactions: " + str(reduced_model_num_reactions))
	print("Trimmed model species: ")
	print(reduced_model.species_names)

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 1

	# Make sure removed species are not included
	assert "H" not in reduced_model.species_names

########
# Input: Artificial mechanism (4 species), removal list of all 4 valid speices
# Output: Reduced mechanism with 0 species
#   *This test fails because Cantera will not produce an empty solution object
########
@pytest.mark.xfail
def testArtRemoveAll():

	# Original model to remove things from
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["H", "H2", "O2", "H2O"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")[1]

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 0
	expected_reactions_num = 0

	# Print stuff
	print("-- Trim all from Artificial Model --")
	print("Trimmed num species: " + str(reduced_model_num_species))
	print("Trimmed num reactions: " + str(reduced_model_num_reactions))
	print("Trimmed model species: ")
	print(reduced_model.species_names)

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 4

	# Make sure removed species are not included
	assert "H" not in reduced_model.species_names

########
# Input: Artificial mechanism (4 species), removal list of one species that is not in the model
# Output: Reduced mechanism with 4 species
########
def testArtRemoveInvalid():

	# Original model to remove things from
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["CH4"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")[1]

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 4
	expected_reactions_num = 2

	# Print stuff
	print("-- Trim 1 invalid from Artificial Model --")
	print("Trimmed num species: " + str(reduced_model_num_species))
	print("Trimmed num reactions: " + str(reduced_model_num_reactions))
	print("Trimmed model species: ")
	print(reduced_model.species_names)

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species())

########
# Input: Artificial mechanism (4 species), removal list of H and CH4 
# Output: Reduced mechanism with 3 species and not the species in the removal list
########
def testArtRemoveInvalidAnd1():

	# Original model to remove things from
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["H", "CH4"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "a-m.cti")[1]

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 3
	expected_reactions_num = 1

	# Print stuff
	print("-- Trim 1 from Artificial Model with invalid --")
	print("Trimmed num species: " + str(reduced_model_num_species))
	print("Trimmed num reactions: " + str(reduced_model_num_reactions))
	print("Trimmed model species: ")
	print(reduced_model.species_names)

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 1

	# Make sure removed species are not included
	assert "H" not in reduced_model.species_names

########
# Input: Random inputs 
# Output: Catch bad input error
# *No error is thrown, so this test fails
########
@pytest.mark.xfail
def testBadInput():

	# Run trim unit
	reduced_model = trim("solution_object", "exclusion_list", "a-m.cti")[1]

########
# Input: GRI Mech 3.0 (53 species), removal list of 10 species
# Output: Reduced mechanism with 43 species and none of the species in the removal list
########
def testGRIMinus10():

	# Original model to remove things from
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	# Create exclusion list for test case
	exclusion_list = ["CH4", "O2", "N2", "H", "OH", "H2O", "CH2", "CH3", "CO", "AR"]

	# Run trim unit
	reduced_model = trim(solution_object, exclusion_list, "gri30.cti")[1]

	#Get number of species/reactions in reduced model	
	reduced_model_num_species = len(reduced_model.species())
	reduced_model_num_reactions = len(reduced_model.reactions())

	# Expected answer	
	expected_species_num = 43
	expected_reactions_num = 14

	# Print stuff
	print("-- Trim 10 from GRI3.0 --")
	print("Trimmed num species: " + str(reduced_model_num_species))
	print("Trimmed num reactions: " + str(reduced_model_num_reactions))
	print("Trimmed model species: ")
	print(reduced_model.species_names)

	# Make sure number matches what is expected
	assert reduced_model_num_species == expected_species_num
	assert reduced_model_num_reactions == expected_reactions_num
	assert reduced_model_num_species == len(solution_object.species()) - 10

	# Make sure removed species are not included
	assert "CH4" not in reduced_model.species_names
	assert "O2" not in reduced_model.species_names
	assert "N2" not in reduced_model.species_names
	assert "H" not in reduced_model.species_names
	assert "OH" not in reduced_model.species_names
	assert "H2O" not in reduced_model.species_names
	assert "CH2" not in reduced_model.species_names
	assert "CH3" not in reduced_model.species_names
	assert "CO" not in reduced_model.species_names
	assert "AR" not in reduced_model.species_names

testGRIMinus10()
