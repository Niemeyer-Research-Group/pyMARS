""" Tests the dictionary lowest model unit used by pyMARS """

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import pytest
import cantera as ct

from sensitivity_analysis import dic_lowest

ROOT_DIR =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def relative_location(file):
    file_path = os.path.join(ROOT_DIR, file)
    return file_path

"""
 NORMAL INPUT DICTIONARY CONTENTS:
    KEY ::    species to remove.
                expressed as a string ("H2O", "N2", etc)
    VALUE ::  error that is introduced by removing just that species.
                expressed as a float, typically between 1 and 15, 
                but with boundary 0 to infinity
""" 

########
# Input: dictionary values that vary in type
# Output: program breaks 
########
@pytest.mark.xfail
def test_bad_input():
    print("---dic_lowest <test 1> input :: dictionary w/ nonsense values") 

    badDic = {"AB": False, "BC": "Darth Vader", "CD": 4387, "DE": 'Z'}
    output = dic_lowest(badDic)     

    print("---dic_lowest <test 1> output :: " + output) 


########
# Input: empty dictionary 
# Output: "error" -- input is caught 
########
def test_empty_input():
    print("---dic_lowest <test 2> input :: empty dictionary")

    emptyDic = {}
    output = dic_lowest(emptyDic)

    print("---dic_lowest <test 2> output :: " + output) 

    assert output == "error"

########
# Input: large lowest value 
# Output: "error" -- incorrect
########
@pytest.mark.xfail
def test_huge_error_value():
    print("---dic_lowest <test 3> input :: huge error value") 

    # lowest value is above upper threshold
    dic = {"AB": 9999999999999999, "BC": 99999999999999, "CD": 10000000000}
    output = dic_lowest(dic)

    print("---dic_lowest <test 3> output :: " + output) 

    assert output == "CD" 


########
# Input: good input where smallest value occurs at end of dictionary 
# Output: lowest value, which is at the end of the dictionary 
########
def test_good_input_1():
    print("---dic_lowest <test 4> input :: good input 1. should be DE")

    dic = {"AB": 10, "BC": 9, "CD": 8, "DE": 5} 
    output = dic_lowest(dic)

    print("---dic_lowest <test 4> output :: " + output) 

    assert output == "DE"

########
# Input: good input where smallest value occurs in the middle of the dictionary
# Output: lowest value, which is in the middle of the dictionary 
########
def test_good_input_2():
    print("---dic_lowest <test 5> input :: good input 2. should be H2O")

    dic = {"AB": 10, "H2O": 1, "CD": 8, "DE": 5} 
    output = dic_lowest(dic)

    print("---dic_lowest <test 5> output :: " + output) 

    assert output == "H2O"

########
# Input: good input where smallest value occurs at beginning of dictionary 
# Output: lowest value, which is at the beginning of the dictionary 
########
def test_good_input_3():
    print("---dic_lowest <test 6> input :: good input 3. should be Keyus")

    dic = {"Keyus": 1, "BC": 9, "CD": 8, "DE": 5} 
    output = dic_lowest(dic)

    print("---dic_lowest <test 6> output :: " + output) 

    assert output == "Keyus"
