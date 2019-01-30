""" Tests the dictionary lowest model unit used by pyMARS """

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

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

# TEST 1 :: nonsense dictionary values passed in =====================
def test_bad_input():
    print("---dic_lowest <test 1> input :: dictionary w/ nonsense values") 

    badDic = {"AB": False, "BC": "Darth Vader", "CD": 4387, "DE": 'Z'}
    output = dic_lowest(badDic)     

    print("---dic_lowest <test 1> output :: " + output) 

# not running test currently because it destroys execution


# TEST 2 :: empty dictionary passed in  ==============================
def test_empty_input():
    print("---dic_lowest <test 2> input :: empty dictionary")

    emptyDic = {}
    output = dic_lowest(emptyDic)

    print("---dic_lowest <test 2> output :: " + output) 

# run test 2
test_empty_input()


# TEST 3 :: Huge error value ========================================= 
def test_huge_error_value():
    print("---dic_lowest <test 3> input :: huge error value") 

    # lowest value is above upper threshold
    dic = {"AB": 9999999999999999, "BC": 99999999999999, "CD": 10000000000}
    output = dic_lowest(dic)

    print("---dic_lowest <test 3> output :: " + output) 

# run test 3
test_huge_error_value()


# TEST 4 :: Good input 1 ============================================== 
def test_good_input_1():
    print("---dic_lowest <test 4> input :: good input 1. should be DE")

    dic = {"AB": 10, "BC": 9, "CD": 8, "DE": 5} 
    output = dic_lowest(dic)

    print("---dic_lowest <test 4> output :: " + output) 

# run test 4
test_good_input_1()


# TEST 5 :: Good input 2 ============================================== 
def test_good_input_2():
    print("---dic_lowest <test 5> input :: good input 2. should be H2O")

    dic = {"AB": 10, "H2O": 1, "CD": 8, "DE": 5} 
    output = dic_lowest(dic)

    print("---dic_lowest <test 5> output :: " + output) 

# run test 5
test_good_input_2()


# TEST 6 :: Good input 3 ============================================== 
def test_good_input_3():
    print("---dic_lowest <test 6> input :: good input 3. should be Caius")

    dic = {"Caius": 1, "BC": 9, "CD": 8, "DE": 5} 
    output = dic_lowest(dic)

    print("---dic_lowest <test 6> output :: " + output) 

# run test 6
test_good_input_3()

