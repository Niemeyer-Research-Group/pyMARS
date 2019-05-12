""" Tests the reduction methods implemented by pyMARS """

import sys
import os
import pytest
import pkg_resources

import cantera as ct

from .. import pfa

PA = {"H":1.5, "H2":2, "H2O":3, "O2":1};
CA = {"H":2, "H2":3, "H2O":5, "O2":4};
PAB = {"H_H2":4, "H_H2O":5, "H_O2":1, "H2_H2O":1, "H2_H":2, "H2_O2":3, "H2O_H2":2, 
	   "H2O_H":3, "H2O_O2":4, "O2_H2":5, "O2_H2O":1,"O2_H":2
	   }
CAB = {"H_H2":1, "H_H2O":2, "H_O2":3, "H2_H2O":4, "H2_H":5, "H2_O2":1, "H2O_H2":5, 
	   "H2O_H":4, "H2O_O2":3, "O2_H2":2, "O2_H2O":1,"O2_H":2
	   }

def relative_location(file):
	file_path = os.path.join(file)
	return pkg_resources.resource_filename(__name__, file_path)


def testGoodInput():
	# Original model
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_1(solution_object,PA,CA,PAB,CAB)
	print(var[0]["H_H2"])
	print(var[0]["H_H2O"])
	print(var[0]["H_O2"])
	assert var[0]["H_H2"] == 2.0
	assert var[0]["H_H2O"] == 2.5
	assert var[0]["H_O2"] == 0.5
	
def testNegativePA():
	PA["H"] = -1.5
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_1(solution_object,PA,CA,PAB,CAB)
	print(var[0])
	assert var[0]["H_H2"] == 2.0
	assert var[0]["H_H2O"] == 2.5
	assert var[0]["H_O2"] == 0.5

#@pytest.mark.xfail	
def testBothNegative():
	PA["H"] = -2
	CA["H"] = -1.5
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_1(solution_object,PA,CA,PAB,CAB)
	print(var[0])
	assert var[0]["H_H2"] == -2.6666666666666665
	assert var[0]["H_H2O"] == -3.3333333333333335
	assert var[0]["H_O2"] == -0.6666666666666666

def testDivideByZero():
	PA["H"] = 0
	CA["H"] = 0
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_1(solution_object,PA,CA,PAB,CAB)
	print(var[0])
	assert var[0]["H_H2"] == 0
	assert var[0]["H_H2O"] == 0
	assert var[0]["H_O2"] == 0
	
@pytest.mark.xfail	
def testPAisChar():
	#PA "H" value will be a char
	PA["H"] = 'f'
	CA["H"] = 2
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_1(solution_object,PA,CA,PAB,CAB)
	print(var[0])

@pytest.mark.xfail	
def testPAMissing():
	#PA will not have the speceis "H" defined.
	PA = {"H2":2, "H2O":3, "O2":1}
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_1(solution_object,PA,CA,PAB,CAB)
	print(var[0])

@pytest.mark.xfail	
def testEmptyDictionary():
	PA = {}
	CA = {}
	PAB = {}
	CAB = {}
	path_to_original = relative_location("artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_1(solution_object,PA,CA,PAB,CAB)
	print(var[0])
