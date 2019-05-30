"""Tests for drg module"""
import sys
import os
import pkg_resources

import pytest
import numpy as np
import networkx as nx
import cantera as ct

from ..sampling import SamplingInputs
from ..pfa import run_pfa, get_rAB_1, get_rAB_2

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


class TestGetRAB1:
    def testGoodInput(self):
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

        # Original model
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)
        print(var[0]['H_H2'])
        print(var[0]['H_H2O'])
        print(var[0]['H_O2'])
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

    #@pytest.mark.xfail	
    def testBothNegative(self):
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

        PA['H'] = -2
        CA['H'] = -1.5
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)
        print(var[0])
        assert var[0]['H_H2'] == -2.6666666666666665
        assert var[0]['H_H2O'] == -3.3333333333333335
        assert var[0]['H_O2'] == -0.6666666666666666

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
        
    @pytest.mark.xfail	
    def testPAisChar(self):
        """PA 'H' value will be a char"""
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
        PA['H'] = 'f'
        CA['H'] = 2
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)
        print(var[0])

    @pytest.mark.xfail	
    def testPAMissing(self):
        """PA will not have the speceis 'H' defined."""
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
        PA = {'H2':2, 'H2O':3, 'O2':1}
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)
        print(var[0])

    @pytest.mark.xfail	
    def testEmptyDictionary(self):
        PA = {}
        CA = {}
        PAB = {}
        CAB = {}
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_1(solution_object,PA,CA,PAB,CAB)
        print(var[0])


class TestGetRAB2:    
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

    @pytest.mark.xfail	
    def testrAB_hasChar(self):
        rAB_p1 = {
            'H_H2':4, 'H_H2O':5, 'H_O2':1, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
		    'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
		    }
        rAB_c1 = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
		    'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
		    }
        rAB_p1['H_O2'] = 'f';
        path_to_original = relative_location('pymars/test/artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_2(solution_object,rAB_p1,rAB_c1)

    @pytest.mark.xfail	
    def testIncompleteDictionary(self):	
        rAB_c1 = {
            'H_H2':1, 'H_H2O':2, 'H_O2':3, 'H2_H2O':4, 'H2_H':5, 'H2_O2':1, 'H2O_H2':5, 
		    'H2O_H':4, 'H2O_O2':3, 'O2_H2':2, 'O2_H2O':1,'O2_H':2
		    }
        rAB_p1 = {
            'H_H2':4, 'H_H2O':5, 'H2_H2O':1, 'H2_H':2, 'H2_O2':3, 'H2O_H2':2, 
            'H2O_H':3, 'H2O_O2':4, 'O2_H2':5, 'O2_H2O':1,'O2_H':2
            }
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_2(solution_object,rAB_p1,rAB_c1)

    @pytest.mark.xfail	
    def testEmptyDictionary(self):
        rAB_p1 = {}
        rAB_c1 = {}
        path_to_original = relative_location('artificial-mechanism.cti')
        solution_object = ct.Solution(path_to_original)
        var = get_rAB_2(solution_object,rAB_p1,rAB_c1)


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
        assert check_equal(reduced_model.species_names, expected_model.species_names)
        assert reduced_model.n_reactions == expected_model.n_reactions
