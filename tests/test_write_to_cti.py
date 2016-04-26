#tests the cti output file
import cantera as ct
import os

def test_solution_object():
    initial_solution_object=ct.Solution('gri30.cti')
    data_file='gri30.cti'
