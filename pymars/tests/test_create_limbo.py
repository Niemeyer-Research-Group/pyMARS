""" Tests the create limbo  model unit used by pyMARS """

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import cantera as ct

from sensitivity_analysis import create_limbo 

ROOT_DIR =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def relative_location(file):
    file_path = os.path.join(ROOT_DIR, file)
    return file_path

"""
 NORMAL INPUT:

    REDUCED_MODEL :: model reduced by the previous reduction
                     cantera solutions object
                     this can be tested with the cti from the example_files folder. either with gri30 (which has
                        53 species) or with artificial mech which has 4 species. 

    EP_STAR       :: epsilon star value
                     float between 0 and 1
                      
    DRGEP_COEFFS  :: dictionary containing strings for species names with accompanying floats between 0 and 1
                     these floats will be checked against the epsilon star value
                     the strings need to be species within the model (cantera solutions object)
                     
    SAFE          :: array of species names (strings) 
                     species in the safe array should not be removed from the reduced model (put in limbo)


GRI30.CTI CONTAINS:

    H2  H   O   O2   OH   H2O   C   CH
    and 45 others...


LIMBO DEFINTION:
    limbo is a list of species that could potentially be removed
    to be in limbo, two primary conditions must be met:
        1. value specified in the drgep dictionary is less than the epsilon star value
        2. the element in not in the "safe" array

""" 

# TEST 1 :: SAFE 1  =====================
def test_safe_1():

    # reduced model
    path_to_reduced = relative_location("example_files/gri30.cti")
    solution_object = ct.Solution(path_to_reduced)

    # drgep dictionary
    dic = {"H2": 0.1, "H": 0.15, "O": 0.2, "O2": 0.25, "OH": 0.3, "H2O": 0.35, "C": 0.4, "CH": 0.45}

    # input conditions
    print("---create_limbo <test 1> input :: ")

    # unchanging conditions
    print("---  reduced_model: gri30.cti H2  H  O  O2  OH  H2O  C  CH ...(45)")
    print("---  drgep_coeffs: ranging 0.1 to 0.45 by increments of .5 for testing purposes.") 

    # variable conditions
    print("---  ep_star: 0.5")
    ep_star = 0.5
    print("---  safe: H, O, OH")
    safe = ["H", "O", "OH"]

    output = create_limbo(solution_object, ep_star, dic, safe)

    print("---create_limbo <test 1> output :: ") 
    print(output)

# run test 1
test_safe_1()

# TEST 2 :: SAFE 2  =====================
def test_safe_2():

    # reduced model
    path_to_reduced = relative_location("example_files/gri30.cti")
    solution_object = ct.Solution(path_to_reduced)

    # drgep dictionary
    dic = {"H2": 0.1, "H": 0.15, "O": 0.2, "O2": 0.25, "OH": 0.3, "H2O": 0.35, "C": 0.4, "CH": 0.45}

    # input conditions
    print("---create_limbo <test 2> input :: ")

    # unchanging conditions
    print("---  reduced_model: gri30.cti H2  H  O  O2  OH  H2O  C  CH ...(45)")
    print("---  drgep_coeffs: ranging 0.1 to 0.45 by increments of .5 for testing purposes.") 

    # variable conditions
    print("---  ep_star: 0.5")
    ep_star = 0.5
    print("---  safe: H2, O2, H2O")
    safe = ["H2", "O2", "H2O"]

    output = create_limbo(solution_object, ep_star, dic, safe)

    print("---create_limbo <test 2> output :: ") 
    print(output)

# run test 2
test_safe_2()

# TEST 3 :: EPSILON STAR 1 =====================
def test_epstar_1():

    # reduced model
    path_to_reduced = relative_location("example_files/gri30.cti")
    solution_object = ct.Solution(path_to_reduced)

    # drgep dictionary
    dic = {"H2": 0.1, "H": 0.15, "O": 0.2, "O2": 0.25, "OH": 0.3, "H2O": 0.35, "C": 0.4, "CH": 0.45}

    # input conditions
    print("---create_limbo <test 3> input :: ")

    # unchanging conditions
    print("---  reduced_model: gri30.cti H2  H  O  O2  OH  H2O  C  CH ...(45)")
    print("---  drgep_coeffs: ranging 0.1 to 0.45 by increments of .5 for testing purposes.") 

    # variable conditions
    print("---  ep_star: 0.05")
    ep_star = 0.05
    print("---  safe: empty list")
    safe = []

    output = create_limbo(solution_object, ep_star, dic, safe)

    print("---create_limbo <test 3> output :: ") 
    print(output)

# run test 3
test_epstar_1()

# TEST 4 :: EPSILON STAR 2 =====================
def test_epstar_2():

    # reduced model
    path_to_reduced = relative_location("example_files/gri30.cti")
    solution_object = ct.Solution(path_to_reduced)

    # drgep dictionary
    dic = {"H2": 0.1, "H": 0.15, "O": 0.2, "O2": 0.25, "OH": 0.3, "H2O": 0.35, "C": 0.4, "CH": 0.45}

    # input conditions
    print("---create_limbo <test 4> input :: ")

    # unchanging conditions
    print("---  reduced_model: gri30.cti H2  H  O  O2  OH  H2O  C  CH ...(45)")
    print("---  drgep_coeffs: ranging 0.1 to 0.45 by increments of .5 for testing purposes.") 

    # variable conditions
    print("---  ep_star: 1")
    ep_star = 1
    print("---  safe: empty list")
    safe = []

    output = create_limbo(solution_object, ep_star, dic, safe)

    print("---create_limbo <test 4> output :: ") 
    print(output)

# run test 4
test_epstar_2()

# TEST 5 :: BAD INPUT 1 =====================
def test_bad_input_1():

    # reduced model
    path_to_reduced = relative_location("example_files/gri30.cti")
    solution_object = ct.Solution(path_to_reduced)

    # drgep dictionary
    dic = {"H2": 0.1, "H": 0.15, "O": 0.2, "O2": 0.25, "OH": 0.3, "H2O": 0.35, "C": 0.4, "CH": 0.45}

    # input conditions
    print("---create_limbo <test 5> input :: ")

    # unchanging conditions
    print("---  reduced_model: gri30.cti H2  H  O  O2  OH  H2O  C  CH ...(45)")
    print("---  drgep_coeffs: bad input. array instead of dictionary") 

    # variable conditions
    print("---  ep_star: 0.5")
    ep_star = 0.5
    print("---  safe: H2, O2, H2O", "O", "OH", "C", "CH")
    safe = ["H2", "O2", "H2O", "O", "OH", "C", "CH"]

    output = create_limbo(solution_object, ep_star, safe, safe)

    print("---create_limbo <test 5> output :: ") 
    print(output)

# run test 5
# test_bad_input_1()
#    breaks stuff

# TEST 6 :: BAD INPUT 2 =====================
def test_bad_input_2():

    # reduced model
    path_to_reduced = relative_location("example_files/gri30.cti")
    solution_object = ct.Solution(path_to_reduced)

    # drgep dictionary
    dic = {"H2": 0.1, "H": 0.15, "O": 0.2, "O2": 0.25, "OH": 0.3, "H2O": 0.35, "C": 0.4, "CH": 0.45}

    # input conditions
    print("---create_limbo <test 6> input :: ")

    # unchanging conditions
    print("---  reduced_model: gri30.cti H2  H  O  O2  OH  H2O  C  CH ...(45)")
    print("---  drgep_coeffs: ranging 0.1 to 0.45 by increments of .5 for testing purposes.") 

    # variable conditions
    print("---  ep_star: 0.5")
    ep_star = 0.5
    print("---  safe: bad input. dictionary instead of array") 
    safe = ["H2", "O2", "H2O", "O", "OH", "C", "CH"]

    output = create_limbo(solution_object, ep_star, dic, dic)

    print("---create_limbo <test 6> output :: ") 
    print(output)

# run test 6
test_bad_input_2()
