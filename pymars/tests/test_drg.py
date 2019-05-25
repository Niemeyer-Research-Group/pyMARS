"""Tests for drg module"""
import sys
import os
import pkg_resources

import pytest

import numpy as np
import networkx as nx
import cantera as ct

from ..drg import graph_search, create_drg_matrix, run_drg

def relative_location(file):
	file_path = os.path.join(file)
	return pkg_resources.resource_filename(__name__, file_path)

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


class TestCreateDRGMatrix:
    """Tests for create_drg_matrix method"""
    
    @pytest.mark.skip
    def testArtificial(self):
        """Uses artificial mechanism to test"""
        # Load model
        path_to_original = relative_location("artificial-mechanism.cti")
        solution_object = ct.Solution(path_to_original)

        # Set up simulation data from givin initial coniditions
        conditions = relative_location("example_input_artificial.txt") # Locate conditions file
        conditions_array = readin_conditions(conditions) # Load conditions
        sim_array = helper.setup_simulations(conditions_array, solution_object) # Set up simulation
        helper.simulate(sim_array) # Run autoignition simulation

        rate_edge_data = drg.get_rates_drg(sim_array, solution_object) # Run unit

        print(rate_edge_data)

        # Pull out timestep one denomenator and numerator dicts
        ic_one = rate_edge_data[list(rate_edge_data.keys())[0]]
        tstep_one = ic_one[list(ic_one.keys())[0]]
        denoms = tstep_one[0]
        numers = tstep_one[1]

        # Expected values for denomenators
        expected_denoms = {}
        expected_denoms["H2O"] = 1.9573216e-13
        expected_denoms["H2"] = .00025854374
        expected_denoms["O2"] = 9.7866081e-14
        expected_denoms["H"] = .00051708749

        assert np.isclose(expected_denoms["H2O"],denoms["H2O"],abs_tol=1.0e-17)
        assert np.isclose(expected_denoms["H2"],denoms["H2"],abs_tol=1.0e-10)
        assert np.isclose(expected_denoms["O2"],denoms["O2"],abs_tol=1.0e-18)
        assert np.isclose(expected_denoms["H"],denoms["H"],abs_tol=1.0e-10)

        expected_numers = {}
        expected_numers["H2O_H2"] = 1.9573216e-13
        expected_numers["H2O_O2"] = 1.9573216e-13
        expected_numers["H2_O2"] = 1.9573216e-13
        expected_numers["H2_H2O"] = 1.9573216e-13
        expected_numers["O2_H2"] = 9.7866081e-14
        expected_numers["O2_H2O"] = 9.7866081e-14
        expected_numers["H2_H"] = .00025854374
        expected_numers["H_H2"] = .00051708749
        
        assert np.isclose(expected_numers["H2O_H2"],numers["H2O_H2"],abs_tol=1.0e-17)
        assert np.isclose(expected_numers["H2O_O2"],numers["H2O_O2"],abs_tol=1.0e-17)
        assert np.isclose(expected_numers["H2_O2"],numers["H2_O2"],abs_tol=1.0e-17)
        assert np.isclose(expected_numers["H2_H2O"],numers["H2_H2O"],abs_tol=1.0e-17)
        assert np.isclose(expected_numers["O2_H2"],numers["O2_H2"],abs_tol=1.0e-18)
        assert np.isclose(expected_numers["O2_H2O"],numers["O2_H2O"],abs_tol=1.0e-18)
        assert np.isclose(expected_numers["H2_H"],numers["H2_H"],abs_tol=1.0e-18)
        assert np.isclose(expected_numers["H_H2"],numers["H_H2"],abs_tol=1.0e-18)


class TestGraphSearch:
    """Tests for graph_search method"""
    #generate test graph
    #starting from A, nodes A,E,C,F,D,I,H,O should be the only nodes found
    def testGraphSearchOneInput(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                ('E','G',0), ('G','I',0), ('G','M',0),
                                ('G','L',1.0), ('E','H',1.0), ('H','J',0)])

        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])


        #temporary solution
        essential_nodes = graph_search(subgraph, 'A')
    
        assert 'A' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in essential_nodes for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]

    #generate test graph
    #starting from A, nodes A,E,C,F,D,I,H,O should be the only nodes found
    def testGraphSearchOneInput2(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                ('E','G',0), ('G','I',0), ('G','M',0),
                                ('G','L',1.0), ('E','H',1.0), ('H','J',0)])

        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])
        
        #temporary solution
        essential_nodes = graph_search(subgraph, 'G')
    
        assert 'G' in essential_nodes
        assert [n not in essential_nodes for n in ['A','B', 'C', 'D', 'J', 'K', 'I', 'O', 'F', 'E', 'H','M', 'N']]
        assert [n in essential_nodes for n in [ 'G', 'L']]

    def testGraphSearch3Inputs(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                ('E','G',0), ('G','I',0), ('G','M',0),
                                ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


        target_species= ['A','C','D']

        essential_nodes = graph_search(graph, target_species)
    
        assert 'A' in essential_nodes
        assert 'C' in essential_nodes
        assert 'D' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in essential_nodes for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]
    
    @pytest.mark.xfail
    def testgraphsearchEmptyGraph (self):
        graph = nx.DiGraph()
        essential_nodes=[]
        assert 'A' not in essential_nodes
        essential_nodes= graph_search(graph, 'A') #this fails can't figure out what the real input to compare to would be
        assert 'A' not in essential_nodes

    @pytest.mark.xfail
    def testGraphshearchwithatargetThatsnotinGraph(self):	
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                ('E','G',0), ('G','I',0), ('G','M',0),
                                ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


        essential_nodes = graph_search(graph, 'Z')
        assert 'Z' in essential_nodes

    def testGraphsearchforinfinteloops(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E'])
        
        graph.add_weighted_edges_from([('A', 'B', 1),('B', 'C', 1), ('C', 'D', 1), ('D','E',1), ('E', 'A', 1)])
        
        essential_nodes= graph_search(graph, 'A')
        for n in essential_nodes:
            print (n)
        assert 'A' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'B', 'E']]
        
    @pytest.mark.xfail
    def testGraphShearchWithATargetThatsNotInGraphAndOneThatIs(self):	
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                ('E','G',0), ('G','I',0), ('G','M',0),
                                ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


        essential_nodes = graph_search(graph, ['B','Z'])
        assert 'B' in essential_nodes

    def testGraphsearchwithListofLength1(self):
        graph = nx.DiGraph()
        graph.add_node('A')


        essential_nodes = graph_search(graph, 'A')
        assert 'A' in essential_nodes
        assert len(essential_nodes)==1

    def testGraphSearchWithTwoOfTheSameItemInTheGraph(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])
    
        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                ('E','G',0), ('G','I',0), ('G','M',0),
                                ('G','L',1.0), ('E','H',1.0), ('H','J',0)])

        essential_nodes = graph_search(graph, 'A')
        assert 'A' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in essential_nodes for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]

    def testGraphSearchWithTwoOfTheSameItemInTheTargetList(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])
    
        graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                                ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                                ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                                ('E','G',0), ('G','I',0), ('G','M',0),
                                ('G','L',1.0), ('E','H',1.0), ('H','J',0)])

        essential_nodes = graph_search(graph, ['A','A'])
        assert 'A' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in essential_nodes for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]
        

class TestRunDRG:
    """Tests driver run_drg method"""
    # Original model
    model_file = "gri30.cti"
    solution_object = ct.Solution(model_file)

	# Conditions for reduction
    conditions = relative_location(os.path.join('inputfiles', 'example_input_file.yaml'))
    error = 5.0
    target_species = ["CH4", "O2"]
    retained_species = ["CH4","O2","N2","H2O","CO2"]

	# Run DRG
    reduced_model = run_drg(
        model_file, conditions, error, target_species, retained_species, 
        )

	# Expected answer
    path_to_answer = relative_location("drg_gri30.cti")
    expected_model = ct.Solution(path_to_answer)
    
    # Make sure models are the same
    assert reduced_model.model.species_names == expected_model.species_names
    assert len(reduced_model.model.reactions()) == len(expected_model.reactions())
