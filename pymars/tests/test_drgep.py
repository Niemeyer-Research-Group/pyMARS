"""Tests for drgep module"""

import sys
import os
import pkg_resources

import pytest
import numpy as np
import cantera as ct
import networkx as nx

from ..sampling import data_files, InputIgnition, InputLaminarFlame
from ..drgep import create_drgep_matrix, graph_search_drgep, get_importance_coeffs
from ..drgep import reduce_drgep, run_drgep, trim

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


def check_equal(list1, list2):
    """Check whether two lists have the same contents (regardless of order).

    Taken from https://stackoverflow.com/a/12813909

    Parameters
    ----------
    list1 : list
        First list, containing all of a particular type
    list2: list
        Second list, containing all of a particular type

    Returns
    -------
    bool
        ``True`` if lists are equal

    """
    return len(list1) == len(list2) and sorted(list1) == sorted(list2)


@pytest.mark.skip(reason='old tests for removed function')
class TestTrimDRGEP:

    ##################
    # Input: The GRI model with a dictionary threshold combo to save 2.  
    # Output: An exclusion list of 51 species
    ##################
    def testGRItrim51(self):

        # Original model
        solution_object = ct.Solution("gri30.cti")

        # Dictionary, retained species, and threshold value to input
        max_dict = {"CH4": 1.0, "N2": .25, "O2": .75}
        threshold_value = .5
        retained_species = []
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 51

        # Print stuff out
        print(" -- Exclude all but CH4 and O2 -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "CH4" not in exclusion_list
        assert "O2" not in exclusion_list

    ##################
    # Input: The artificial model with a dictionary threshold combo to remove 1.  
    # Output: An exclusion list of 1 species
    ##################
    def testArtTrim1(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {"H": 1.0, "H2": .2, "O2": 1, "H2O": 1}
        threshold_value = .5
        retained_species = []
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 1

        # Print stuff out
        print(" -- Exclude only H2 -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list

    ##################
    # Input: The artificial model with a dictionary threshold combo to remove 2.  
    # Output: An exclusion list of 2 species
    ##################
    def testArtTrim2(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {"H": 1.0, "H2": .2, "O2": .45, "H2O": 1}
        threshold_value = .5
        retained_species = []
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 2

        # Print stuff out
        print(" -- Exclude H2 and O2 -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list
        assert "O2" in exclusion_list

    ##################
    # Input: The artificial model with a dictionary threshold combo to remove 2. 
    #        However, the retianed species list saves one of them from being removed. 
    # Output: An exclusion list of 1 species
    ##################
    def testArtRetain(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {"H": 1.0, "H2": .2, "O2": .45, "H2O": 1}
        threshold_value = .5
        retained_species = ["O2"]
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 1

        # Print stuff out
        print(" -- Exclude H2, retain O2 -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list

    ##################
    # Input: The artificial model with a dictionary threshold combo to remove 4.  
    # Output: An exclusion list of 4 species
    ##################
    def testArtRemoveAll(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {"H": 0, "H2": 0, "O2": 0, "H2O": 0}
        threshold_value = .1
        retained_species = []
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 4

        # Print stuff out
        print(" -- Exclude all -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list
        assert "H2O" in exclusion_list
        assert "H" in exclusion_list
        assert "O2" in exclusion_list

    ##################
    # Input: The artificial model with a dictionary value equal to the threshold.  
    # Output: The species is removed.
    ##################
    def testArtRemoveEqual(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {"H": 1.0, "H2": .2, "O2": 1, "H2O": 1}
        threshold_value = .2
        retained_species = []
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 1

        # Print stuff out
        print(" -- OIC == Threshold -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list

    ##################
    # Input: Dictionary is empty.  
    # Output: Nothing.
    ##################
    def testEmptyDict(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {}
        threshold_value = .5
        retained_species = ["O2"]
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 3

        # Print stuff out
        print(" -- Empty Dict -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list
        assert "H2O" in exclusion_list
        assert "H" in exclusion_list


    ##################
    # Input: Dictionary is mixed types.  
    # Output: Doesn't break, still works.
    ##################
    def testInvalidDict(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {"LK": False, 3: 7.5, "H2O": True}
        threshold_value = .5
        retained_species = ["O2"]
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 2

        # Print stuff out
        print(" -- Invalid Input -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list
        assert "H" in exclusion_list

    ##################
    # Input: Threshold is greater than 1.  
    # Output: Everything is removed.
    ##################
    def testThresholdGreaterThan1(self):

        # Original model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Dictionary, retained species, and threshold value to input
        max_dict = {"H": 1.0, "H2": .2, "O2": 1, "H2O": 1}
        threshold_value = 5
        retained_species = []
        done = [False]	

        # Run unit
        exclusion_list = trim_drgep(max_dict, solution_object, threshold_value, retained_species, done)

        expected_list_len = 4

        # Print stuff out
        print(" -- Threshold Greater Than 1 -- ")
        print("Exclusion list length: " + str(len(exclusion_list)))
        print(" ")
        print(exclusion_list)

        # Assert list is correct length and doesn't include keepers
        assert len(exclusion_list) == expected_list_len
        assert "H2" in exclusion_list
        assert "H2O" in exclusion_list
        assert "H" in exclusion_list
        assert "O2" in exclusion_list


class TestGraphSearchDRGEP:
    def test_DRGEP_GraphSearchSimple(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'
            ])

        graph.add_weighted_edges_from([
            ('A','B', 1), ('B','C',1), ('C','D',1),
            ('D','E',1), ('E','F',1)
            ])

        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])

        #temporary solution
        max_dic = graph_search_drgep(subgraph, 'A')

        assert 'A' in max_dic
        assert max_dic['A']==1
        assert max_dic['B']==1
        assert max_dic['C']==1
        assert max_dic['D']==1
        assert max_dic['E']==1
        assert max_dic['F']==1

    def test_DRGEP_GraphSearchDouble(self):
        """Starting from A, nodes A-F should be found doubling each time
        """
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 
                            'L', 'M', 'N', 'O'
                            ])

        graph.add_weighted_edges_from([('A','B', 2), ('B','C',2), ('C','D',2),
                                    ('D','E',2), ('E','F',2)
                                    ])

        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])


        #temporary solution
        max_dic = graph_search_drgep(subgraph, 'A')

        assert 'A' in max_dic
        assert max_dic['A']==1
        assert max_dic['B']==2
        assert max_dic['C']==4
        assert max_dic['D']==8
        assert max_dic['E']==16
        assert max_dic['F']==32

    def test_DRGEP_GraphSearchHalve(self):
        """starting from A, nodes A-F should be found decreasing by half each time
        """
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 
                            'L', 'M', 'N', 'O'
                            ])

        graph.add_weighted_edges_from([('A','B', .5), ('B','C',.5), ('C','D',.5), 
                                    ('D','E',.5), ('E','F',.5)
                                    ])


        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])


        #temporary solution
        max_dic = graph_search_drgep(subgraph, 'A')

        assert 'A' in max_dic
        assert max_dic['A']==1
        assert max_dic['B']==.5
        assert max_dic['C']==.25
        assert max_dic['D']==.125
        assert max_dic['E']==.0625
        assert max_dic['F']==.03125

    def test_DRGEP_GraphSearchFromG(self):
        """starting from G, nodes G and L should be the only nodes found
        """
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 
                            'L', 'M', 'N', 'O'
                            ])

        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                    ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                    ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                    ('E','G',0), ('G','I',0), ('G','M',0),
                                    ('G','L',1.0), ('E','H',1.0), ('H','J',0)
                                    ])

        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])


        #temporary solution
        max_dic = graph_search_drgep(subgraph, 'G')

        assert 'G' in max_dic
        assert max_dic['G']
        assert [n in max_dic for n in ['G', 'L']]
        assert [n not in max_dic for n in ['A','B','C','D','E','F','H','I', 'J', 'K', 'L','M', 'N']]
        assert max_dic['G']==1
        assert max_dic['L']==1

    def test_DRGEP_GraphSearchComplicated(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        graph.add_weighted_edges_from([('A','B', .1), ('C','D',.2), ('B','C',.3), ('D','E',.6),
                                    ('E','F',.5), ('F','A',.4), ('A','E',.7), ('C','G',.01), 
                                    ('G','A',.01), ('E','B',.05)
                                    ])


        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])


        #temporary solution
        max_dic = graph_search_drgep(subgraph, 'A')
        
        assert 'A' in max_dic
        assert [n in max_dic for n in ['A','B','C','D','E','F','G']]
        assert max_dic['A']==1
        assert max_dic['B']==.1
        assert max_dic['C']==0.03
        assert max_dic['D']==.006
        assert max_dic['E']==.7
        assert max_dic['F']==.35
        assert max_dic['G']==.0003

    def test_DRGEP_GraphSearch3Inputs(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'
            ])

        graph.add_weighted_edges_from([
            #('A','F', 0), ('A','N',0),
            ('C','F',1.0), ('A','C',1.0), ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

        target_species= ['A', 'C', 'D']

        max_dic = graph_search_drgep(graph, target_species)
    
        assert 'A' in max_dic
        assert 'C' in max_dic
        assert 'D' in max_dic
        assert [n in max_dic for n in [
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I','J','K', 'L','M','N','O',
            ]]
        assert [n not in max_dic for n in ['B']]

    def test_DRGEP_GraphSearchComplicated2Inputs(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G'])
        graph.add_weighted_edges_from([
            ('A','B', .5),('G','C',.9),('C','D',.5),('D','E',.5),
            ('E','F',.5),('B','F',0.0002),('A','C',.0225)
            ])

        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])

        #temporary solution
        max_dic = graph_search_drgep(subgraph, ['A','G'])
        
        assert 'A' in max_dic
        assert [n in max_dic for n in ['A','B','C','D','E','F','G']]
        assert max_dic['A']==1
        assert max_dic['B']==.5
        assert max_dic['C']==0.9
        assert max_dic['D']==.45
        assert max_dic['E']==.225
        assert max_dic['F']==.1125
        assert max_dic['G']==1
    
    @pytest.mark.xfail
    def test_DRGEP_graphsearchEmptyGraph(self):
        graph = nx.DiGraph()
        max_dic=[]
        assert 'A' not in max_dic
        max_dic= graphsearch(graph, 'A') #this fails can't figure out what the real input to compare to would be
        assert 'A' not in max_dic

    @pytest.mark.xfail
    def test_DRGEP_GraphshearchwithatargetThatsnotinGraph(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'
            ])

        graph.add_weighted_edges_from([
            ('A','F', 0), ('C','F',1.0), ('A','C',1.0),
            ('A','N',0), ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            ('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), ('H','J',0)
            ])

        max_dic = graph_search_drgep(graph, 'Z')
        assert 'Z' in max_dic

    def test_DRGEP_Graphsearchforinfinteloops(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E'])
        
        graph.add_weighted_edges_from([
            ('A', 'B', 1),('B', 'C', 1), ('C', 'D', 1), ('D','E',1), ('E', 'A', 1)
            ])
        
        max_dic= graph_search_drgep(graph, 'A')
        for n in max_dic:
            print (n)
        assert 'A' in max_dic
        assert [n in max_dic for n in ['A', 'C', 'D', 'B', 'E']]
        
    @pytest.mark.xfail
    def test_DRGEP_GraphShearchWithATargetThatsNotInGraphAndOneThatIs(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'
            ])

        graph.add_weighted_edges_from([
            ('A','F', 0), ('C','F',1.0), ('A','C',1.0),
            ('A','N',0), ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            ('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), ('H','J',0)
            ])

        max_dic = graph_search_drgep(graph, ['B','Z'])
        assert 'B' in max_dic

    def test_DRGEP_GraphsearchwithListofLength1(self):
        graph = nx.DiGraph()
        graph.add_node('A')

        max_dic = graph_search_drgep(graph, 'A')
        assert 'A' in max_dic
        assert len(max_dic)==1
        assert max_dic['A']==1

    def test_DRGEP_GraphSearchWithTwoOfTheSameNodeInTheGraph(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([
            'A', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'
            ])
    
        graph.add_weighted_edges_from([
            ('A','F', 0), ('C','F',1.0), ('A','C',1.0),
            ('A','N',0), ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            ('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), ('H','J',0)
            ])

        max_dic = graph_search_drgep(graph, 'A')
        assert 'A' in max_dic
        assert [n in max_dic for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in max_dic for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]

    def test_DRGEP_GraphSearchWithTwoOfTheSameItemInTheTargetList(self):
        graph = nx.DiGraph()
        graph.add_nodes_from([
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'
            ])
    
        graph.add_weighted_edges_from([
            ('A','F', 0), ('C','F',1.0), ('A','C',1.0),
            ('A','N',0), ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            ('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), ('H','J',0)
            ])

        max_dic = graph_search_drgep(graph, ['A','A'])
        assert 'A' in max_dic
        assert [n in max_dic for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in max_dic for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]


class TestGetImportanceCoeffs:
    def test_one_condition_multiple_targets(self):
        """Tests correct coefficients with one condition but multiple targets"""
        matrices = [
            np.array([[0, 0.6, 0], [1.0, 0, 1.0], [0, 0.1, 0]])
            ]
        coefficients = get_importance_coeffs(['A', 'B', 'C'], ['A', 'C'], matrices)
        assert coefficients['A'] == 1.0
        assert coefficients['B'] == 0.6
        assert coefficients['C'] == 1.0

        matrices = [
            np.array([[0, 0.1, 0], [1.0, 0, 1.0], [0, 0.6, 0]])
            ]
        coefficients = get_importance_coeffs(['A', 'B', 'C'], ['A', 'C'], matrices)
        assert coefficients['A'] == 1.0
        assert coefficients['B'] == 0.6
        assert coefficients['C'] == 1.0

    def test_multiple_conditions(self):
        """Tests correct coefficients with multiple conditions"""
        matrices = [
            np.array([[0, 0.6], [1.0, 0.0]]),
            np.array([[0, 0.0], [1.0, 0.0]]),
            np.array([[0, 0.1], [1.0, 0.0]]),
            ]
        coefficients = get_importance_coeffs(['A', 'B'], ['A'], matrices)
        assert coefficients['A'] == 1.0
        assert coefficients['B'] == 0.6


class TestRunDRGEP:
    def test_gri_reduction(self):
        """Tests driver run_drgep method"""
        model_file = 'gri30.cti'

        # Conditions for reduction
        ignition_conditions = [
            InputIgnition(
                kind='constant volume', pressure=1.0, temperature=1000.0, equivalence_ratio=1.0,
                fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}
                ),
            InputIgnition(
                kind='constant volume', pressure=1.0, temperature=1200.0, equivalence_ratio=1.0,
                fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}
                ),
                ]

        data_files['output_ignition'] = relative_location(
            os.path.join('assets', 'example_ignition_output.txt')
            )
        data_files['data_ignition'] = relative_location(
            os.path.join('assets', 'example_ignition_data.dat')
            )
        error = 5.0

        # Run DRG
        with TemporaryDirectory() as temp_dir:
            reduced_model = run_drgep(
                model_file, error, ['CH4', 'O2'], ['N2'], ignition_conditions,
                num_threads=1, path=temp_dir
                )

        # Expected answer
        expected_model = ct.Solution(relative_location(os.path.join('assets', 'drgep_gri30.cti')))
        
        # Make sure models are the same
        assert check_equal(reduced_model.model.species_names, expected_model.species_names)
        assert reduced_model.model.n_reactions == expected_model.n_reactions
        assert round(reduced_model.error, 2) == 3.22
    
    def test_dri_reductionflame(self):
        """Tests driver run_drgep method with only flame input data"""
        model_file = 'gri30.cti'

        flame_conditions = [
            InputLaminarFlame(
                kind='constant pressure', pressure=1.0, temperature=1000.0, equivalence_ratio=1.0,
                fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}, width=0.01
                ),
            InputLaminarFlame(
                kind='constant pressure', pressure=1.0, temperature=1200.0, equivalence_ratio=1.0,
                fuel={'CH4': 1.0}, oxidizer={'O2': 1.0, 'N2': 3.76}, width=0.1
                ),
        ]
        
        data_files['output_flame'] = relative_location(
            os.path.join('assets', 'example_flame_output.txt')
            )
        data_files['data_flame'] = relative_location(
            os.path.join('assets', 'example_flame_data.dat')
            )
        error = 5.0

        # Run DRG
        with TemporaryDirectory() as temp_dir:
            reduced_model = run_drgep(
                model_file, error, ['CH4', 'O2'], ['N2'],  flame_conditions=flame_conditions,
                num_threads=1, path=temp_dir
                )

        # Expected answer
        expected_model = ct.Solution(relative_location(os.path.join('assets', 'drgep_gri30.cti')))
        
        # Make sure models are the same
        #assert check_equal(reduced_model.model.species_names, expected_model.species_names)
        #assert reduced_model.model.n_reactions == expected_model.n_reactions
        #assert round(reduced_model.error, 2) == 3.22
