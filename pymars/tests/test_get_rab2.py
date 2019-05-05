""" Tests the reduction methods implemented by pyMARS """

import sys
import os
import pytest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import cantera as ct

import pfa

ROOT_DIR =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rAB_p1 = {"H_H2":4, "H_H2O":5, "H_O2":1, "H2_H2O":1, "H2_H":2, "H2_O2":3, "H2O_H2":2, "H2O_H":3, "H2O_O2":4, "O2_H2":5, "O2_H2O":1,"O2_H":2}
rAB_c1 = {"H_H2":1, "H_H2O":2, "H_O2":3, "H2_H2O":4, "H2_H":5, "H2_O2":1, "H2O_H2":5, "H2O_H":4, "H2O_O2":3, "O2_H2":2, "O2_H2O":1,"O2_H":2}

def relative_location( file):
	file_path = os.path.join(ROOT_DIR, file)
	return file_path

#get_rAB2(new_solu,rab_p1,rab_c1)

def testGoodInput():

	# Original model
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_2(solution_object,rAB_p1,rAB_c1);
	print(var[0]["H_H2"]);
	print(var[0]["H_H2O"]);
	print(var[0]["H_O2"]);
	assert var[0]["H_H2"] == 15
	assert var[0]["H_H2O"] == 5
	assert var[0]["H_O2"] == 32
	
def testGoodInputZero():
	rAB_p1["H_O2"] = 0;
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_2(solution_object,rAB_p1,rAB_c1);
	print(var[0]["H_H2"]);
	print(var[0]["H_H2O"]);
	print(var[0]["H_O2"]);
	assert var[0]["H_H2"] == 10
	assert var[0]["H_H2O"] == 4
	assert var[0]["H_O2"] == 32

def testGoodInputNegative():
	rAB_p1["H_O2"] = -3;
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_2(solution_object,rAB_p1,rAB_c1);
	print(var[0]["H_H2"]);
	print(var[0]["H_H2O"]);
	print(var[0]["H_O2"]);
	assert var[0]["H_H2"] == -5
	assert var[0]["H_H2O"] == 1
	assert var[0]["H_O2"] == 32

@pytest.mark.xfail	
def testrAB_hasChar():
	rAB_p1["H_O2"] = 'f';
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_2(solution_object,rAB_p1,rAB_c1);

@pytest.mark.xfail	
def testIncompleteDictionary():	
	rAB_p1 = {"H_H2":4, "H_H2O":5, "H2_H2O":1, "H2_H":2, "H2_O2":3, "H2O_H2":2, "H2O_H":3, "H2O_O2":4, "O2_H2":5, "O2_H2O":1,"O2_H":2}
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_2(solution_object,rAB_p1,rAB_c1);

@pytest.mark.xfail	
def testEmptyDictionary():
	rAB_p1 = {}
	rAB_c1 = {}
	path_to_original = relative_location("pymars/tests/artificial-mechanism.cti")
	solution_object = ct.Solution(path_to_original)
	var = pfa.get_rAB_2(solution_object,rAB_p1,rAB_c1);

