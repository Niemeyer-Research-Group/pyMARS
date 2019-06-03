"""Tests for drg module"""
import sys
import os
import pkg_resources

import pytest
import numpy as np
import networkx as nx
import cantera as ct

from ..sampling import SamplingInputs
from ..pfa import run_pfa

# Taken from http://stackoverflow.com/a/22726782/1569494
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from contextlib import contextmanager
    import shutil
    import tempfile
    import errno

    @contextmanager
    def TemporaryDirectory():
        name = tempfile.mkdtemp()
        try:
            yield name
        finally:
            try:
                shutil.rmtree(name)
            except OSError as e:
                # Reraise unless ENOENT: No such file or directory
                # (ok if directory has already been deleted)
                if e.errno != errno.ENOENT:
                    raise


def relative_location(file):
    file_path = os.path.join(file)
    return pkg_resources.resource_filename(__name__, file_path)


class TestCreatePFAMatrix:
    def test_artificial_mech(self):
        # species: H, H2, H2O, O2
        PA = np.array({'H':1.5, 'H2':2, 'H2O':3, 'O2':1}
        CA = {'H':2, 'H2':3, 'H2O':5, 'O2':4}
        PAB = {
            'H_H2':4, 'H_H2O':5, 'H_O2':1, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
            'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
            }
        CAB = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
            'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
            }

        # Original model
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)

        assert var[0]['H_H2'] == 2.0
        assert var[0]['H_H2O'] == 2.5
        assert var[0]['H_O2'] == 0.5
        
    def testNegativePA(self):
        PA = {'H':1.5, 'H2':2, 'H2O':3, 'O2':1}
        CA = {'H':2, 'H2':3, 'H2O':5, 'O2':4}
        PAB = {
            'H_H2':4, 'H_H2O':5, 'H_O2':1, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
            'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
            }
        CAB = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
            'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
            }

        PA['H'] = -1.5
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)
        print(var[0])
        assert var[0]['H_H2'] == 2.0
        assert var[0]['H_H2O'] == 2.5
        assert var[0]['H_O2'] == 0.5

    def testDivideByZero(self):
        PA = {'H':1.5, 'H2':2, 'H2O':3, 'O2':1}
        CA = {'H':2, 'H2':3, 'H2O':5, 'O2':4}
        PAB = {
            'H_H2':4, 'H_H2O':5, 'H_O2':1, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
            'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
            }
        CAB = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
            'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
            }

        PA['H'] = 0
        CA['H'] = 0
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)
        print(var[0])
        assert var[0]['H_H2'] == 0
        assert var[0]['H_H2O'] == 0
        assert var[0]['H_O2'] == 0
    
    def testGoodInput(self):
        rAB_p1 = {
            'H_H2':4, 'H_H2O':5, 'H_O2':1, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
            'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
            }
        rAB_c1 = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
            'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
            }
        # Original model
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_2(solution_object, rAB_p1, rAB_c1)
        print(var[0]['H_H2'])
        print(var[0]['H_H2O'])
        print(var[0]['H_O2'])
        assert var[0]['H_H2'] == 15
        assert var[0]['H_H2O'] == 5
        assert var[0]['H_O2'] == 32
        
    def testGoodInputZero(self):
        rAB_p1 = {
            'H_H2':4, 'H_H2O':5, 'H_O2':1, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
            'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
            }
        rAB_c1 = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
            'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
            }
        rAB_p1['H_O2'] = 0;
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_2(solution_object,rAB_p1,rAB_c1)
        print(var[0]['H_H2'])
        print(var[0]['H_H2O'])
        print(var[0]['H_O2'])
        assert var[0]['H_H2'] == 10
        assert var[0]['H_H2O'] == 4
        assert var[0]['H_O2'] == 32

    def testGoodInputNegative(self):
        rAB_p1 = {
            'H_H2':4, 'H_H2O':5, 'H_O2':1, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
            'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
            }
        rAB_c1 = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
            'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
            }
        rAB_p1['H_O2'] = -3;
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_2(solution_object,rAB_p1,rAB_c1)
        print(var[0]['H_H2'])
        print(var[0]['H_H2O'])
        print(var[0]['H_O2'])
        assert var[0]['H_H2'] == -5
        assert var[0]['H_H2O'] == 1
        assert var[0]['H_O2'] == 32


class TestRunPFA:
    def testPFA(self):
        # Original model
        solution_object = ct.Solution('gri30.cti')

        # Conditions for reduction	
        conditions = relative_location('example_input_file.txt')
        error = 5.0
        target_species = ["CH4","O2"]
        retained_species = ["CH4","O2","N2","H2O","CO2"]
        final_error = [0]

        # Run PFA
        reduced_model = run_pfa(
            solution_object, conditions, error, target_species, 
            retained_species, "gri30.cti", final_error
            )

        # Expected answer	
        path_to_answer = relative_location("pfa_gri30.cti")
        expected_model = ct.Solution(path_to_answer)

        # Make sure models are the same
        assert check_equal(reduced_model.model.species_names, expected_model.species_names)
        assert reduced_model.model.n_reactions == expected_model.n_reactions
