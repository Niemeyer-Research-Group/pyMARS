""" Tests the create output file unit used by pyMARS """

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import pytest
import cantera as ct

from create_trimmed_model import trim
from soln2cti import write

ROOT_DIR =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def relative_location(file):
	file_path = os.path.join(ROOT_DIR, file)
	return file_path

########
# Input: Solution representing the GRI 3.0 Mech
# Output: A file storing the information for the GRI 3.0 Mech
########
def testGRIwrite():
	
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	write(solution_object)

	output_path = "pym_gri30.cti"
	new_solution_object = ct.Solution(output_path)

	assert solution_object.species_names == new_solution_object.species_names
	assert len(solution_object.reactions()) == len(new_solution_object.reactions())

########
# Input: Solution representing the Artificial Mechanism
# Output: A file storing the information for the Artificial Mechanism
########
def testArtWrite():
	
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	write(solution_object)

	output_path = "pym_gas.cti"
	new_solution_object = ct.Solution(output_path)

	assert solution_object.species_names == new_solution_object.species_names
	assert len(solution_object.reactions()) == len(new_solution_object.reactions())

 
########
# Input: Solution representing the GRI 3.0 Mech slightly reduced
# Output: A file storing the information for the GRI 3.0 Mech slightly reduced
########
def testGRIwriteRed():
	
	path_to_original = relative_location("example_files/gri30.cti")
	solution_object = ct.Solution(path_to_original)

	exclusion_list = ["CH4", "O2", "N2"]
	solution_object = trim(solution_object, exclusion_list, "gri30.cti")[1]
	
	write(solution_object)

	output_path = "pym_trimmed_gri30.cti"
	new_solution_object = ct.Solution(output_path)

	path_to_original = relative_location("pymars/tests/eout_gri30.cti")
	solution_object = ct.Solution(path_to_original)

	assert solution_object.species_names == new_solution_object.species_names
	assert len(solution_object.reactions()) == len(new_solution_object.reactions())

########
# Input: Solution representing the Artificial Mechanism slightly reduced
# Output: A file storing the information for the Artificial Mech slightly reduced
########
def testArtWriteRed():
	
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	exclusion_list = ["H", "O2"]
	solution_object = trim(solution_object, exclusion_list, "gas.cti")[1]
	
	write(solution_object)

	output_path = "pym_trimmed_gas.cti"
	new_solution_object = ct.Solution(output_path)

	path_to_original = relative_location("pymars/tests/eout_artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)

	assert solution_object.species_names == new_solution_object.species_names
	assert len(solution_object.reactions()) == len(new_solution_object.reactions())


########
# Input: Garbage input
# Output: Program catches the error
#  *Note that this fails because the funciton does not handle invalid input in any way
########
@pytest.mark.xfail
def testBadInput():
	
	write("garbage")
