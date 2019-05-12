""" Tests the reduction methods implemented by pyMARS """

import sys
import os
import pkg_resources
import pytest

import cantera as ct
from .. import readin_initial_conditions

def relative_location(file):
	file_path = os.path.join(file)
	return pkg_resources.resource_filename(__name__, file_path)


#@pytest.mark.xfail
def testGoodInput1():
	#test reading in example_input_file
	conditions =  relative_location("example_input_file.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	print("---------------------------------------------------------------");	
	print("condObject[1].pressure: " + condObject[1].pressure)
	print("condObject[1].temperature: " + condObject[1].temperature)
	print("condObject[1].species: ", condObject[1].species)
	print("condObject[1].fuel: " + condObject[1].fuel)
	print("condObject[1].oxid: " + condObject[1].oxid)
	
	assert condObject[0].pressure == " 1.0\n"
	assert condObject[0].temperature == " 1000\n"
	assert condObject[0].species == {'CH4': '1.0', 'N2': '3.76', 'O2': '1.0'}
	assert condObject[0].fuel == "CH4:1.0"
	assert condObject[0].oxid == 'O2:1.0,N2:3.76'
	
	assert condObject[1].pressure == " 1.0\n"
	assert condObject[1].temperature == " 1200\n"
	assert condObject[1].species == {'CH4': '1.0', 'N2': '3.76', 'O2': '1.0'}
	assert condObject[1].fuel == "CH4:1.0"
	assert condObject[1].oxid == 'O2:1.0,N2:3.76'

#@pytest.mark.xfail
def testGoodInput2():
	#test reading in example_input_artificial
	conditions =  relative_location("example_input_artificial.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	
	assert condObject[0].pressure == " 1.0\n"
	assert condObject[0].temperature == " 500\n"
	assert condObject[0].species == {'H': '1.0', 'O2': '1.0'}
	assert condObject[0].fuel == "H:1.0"
	assert condObject[0].oxid == 'O2:1.0'

#@pytest.mark.xfail
def testRearrangedInput():
	#test reading in re-arranged example_input_file
	conditions =  relative_location("inputfiles/example_input_file_new_order.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	print("---------------------------------------------------------------");	
	print("condObject[1].pressure: " + condObject[1].pressure)
	print("condObject[1].temperature: " + condObject[1].temperature)
	print("condObject[1].species: ", condObject[1].species)
	print("condObject[1].fuel: " + condObject[1].fuel)
	print("condObject[1].oxid: " + condObject[1].oxid)
	
	assert condObject[0].pressure == " 1.0\n"
	assert condObject[0].temperature == " 1000\n"
	assert condObject[0].species == {'CH4': '1.0', 'N2': '3.76', 'O2': '1.0'}
	assert condObject[0].fuel == "CH4:1.0"
	assert condObject[0].oxid == 'O2:1.0,N2:3.76'
	
	assert condObject[1].pressure == " 1.0\n"
	assert condObject[1].temperature == " 1200\n"
	assert condObject[1].species == {'CH4': '1.0', 'N2': '3.76', 'O2': '1.0'}
	assert condObject[1].fuel == "CH4:1.0"
	assert condObject[1].oxid == 'O2:1.0,N2:3.76'

#@pytest.mark.xfail
def testIllogicalSpeciesNames():
	#test reading file in correct order and format, but illogical values
	conditions =  relative_location("inputfiles/example_input_file_illogical.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	
	assert condObject[0].pressure == " hello\n"
	assert condObject[0].temperature == " I\n"
	assert condObject[0].species == {'everthing': '5.0', 'oxidizer': '1.0', 'peanutbutter': '3.7'}
	assert condObject[0].fuel == "everthing:5.0"
	assert condObject[0].oxid == 'oxidizer:1.0,peanutbutter:3.7'
	

@pytest.mark.xfail
def testStringForFloats():
	#test string values when expecting float
	conditions =  relative_location("inputfiles/example_input_file_bad_float.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	
	assert condObject[0].pressure == " 1.0\n"
	assert condObject[0].temperature == " 1000\n"
	assert condObject[0].species == {'CH4': 'youtube', 'N2': '4.5.6', 'O2': 'gello'}
	assert condObject[0].fuel == "CH4:youtube"
	assert condObject[0].oxid == 'O2:gello,N2:4.5.6'
	
@pytest.mark.xfail
def testNoValuesAfterSpecies():
	conditions =  relative_location("inputfiles/example_input_file_no_species_value.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	
	assert condObject[1].pressure == " 1.0\n"
	assert condObject[1].temperature == " 1000\n"
	assert condObject[1].species == {'CH4', 'N2', 'O2'}
	assert condObject[1].fuel == "CH4"
	assert condObject[1].oxid == 'O2,N2'

def testNegativeSpeciesValues():
	#test reading in negative species values
	conditions =  relative_location("inputfiles/example_input_file_negative_species_values.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	
	assert condObject[0].pressure == " 1.0\n"
	assert condObject[0].temperature == " 1000\n"
	assert condObject[0].species == {'CH4': '-1.0', 'N2': '-3.0', 'O2': '-2.0'}
	assert condObject[0].fuel == "CH4:-1.0"
	assert condObject[0].oxid == 'O2:-2.0,N2:-3.0'

@pytest.mark.xfail
def testNoOpeningKeyword():
	#test file with no "CONV"
	#should not read anything.
	conditions =  relative_location("inputfiles/example_input_file_no_open_keyword.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	
	assert condObject[0].pressure == " 1.0\n"
	assert condObject[0].temperature == " 1000\n"
	assert condObject[0].species == {'CH4': '1.0', 'N2': '3.76', 'O2': '1.0'}
	assert condObject[0].fuel == "CH4:1.0"
	assert condObject[0].oxid == 'O2:1.0,N2:3.76'

@pytest.mark.xfail
def testNoClosingKeyword():
	#test file with no "END"
	#Will read, but will never close.
	conditions =  relative_location("inputfiles/example_input_file_no_close_keyword.txt")
	condObject = readin_initial_conditions.readin_conditions(conditions)
	print("condObject[0].pressure: " + condObject[0].pressure)
	print("condObject[0].temperature: " + condObject[0].temperature)
	print("condObject[0].species: ", condObject[0].species)
	print("condObject[0].fuel: " + condObject[0].fuel)
	print("condObject[0].oxid: " + condObject[0].oxid)
	
	assert condObject[0].pressure == " 1.0\n"
	assert condObject[0].temperature == " 1000\n"
	assert condObject[0].species == {'CH4': '1.0', 'N2': '3.76', 'O2': '1.0'}
	assert condObject[0].fuel == "CH4:1.0"
	assert condObject[0].oxid == 'O2:1.0,N2:3.76'
#testNoClosingKeyword()
