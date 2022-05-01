"""Tests for pfa module"""
import sys
import os
import pkg_resources

import pytest

import numpy as np
import networkx as nx
import cantera as ct

from ..sampling import data_files, InputIgnition, InputLaminarFlame
from ..pfa import graph_search, create_pfa_matrix, run_pfa, reduce_pfa

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


class TestCreatePFAMatrix:
    """Tests for create_pfa_matrix method"""

    def test_qss_artificial(self):
        """Test using four species artificial model with QSS species from 2006 DRG paper.
        
        # R \approx F / 1e3
        """
        R1 = ct.Reaction.fromCti('''reaction('F => R', [1.0, 0.0, 0.0])''')
        R2 = ct.Reaction.fromCti('''reaction('R => P', [1.0e3, 0.0, 0.0])''')
        R3 = ct.Reaction.fromCti('''reaction('R => Pp', [1.0, 0.0, 0.0])''')

        F = ct.Species('F', 'H:1')
        R = ct.Species('R', 'H:1')
        P = ct.Species('P', 'H:1')
        Pp = ct.Species('Pp', 'H:1')
        for sp in [F, R, P, Pp]:
            sp.thermo = ct.ConstantCp(
                300, 1000, 101325, (300, 1.0, 1.0, 1.0)
                )
        model = ct.Solution(
            thermo='IdealGas', kinetics='GasKinetics',
            species=[F, R, P, Pp], reactions=[R1, R2, R3]
            )
        mass_fracs = [1., 1./1.e3, 0., 0.]
        state = 1000, ct.one_atm, mass_fracs
        matrix = create_pfa_matrix(state, model)
        denom = max(mass_fracs[0], (1e3 + 1) * mass_fracs[1])
        correct = np.array([
            [0, 1.0, 1e3*mass_fracs[1] / denom, mass_fracs[1] / denom],
            [mass_fracs[0] / denom, 0, 1e3*mass_fracs[1] / denom, mass_fracs[1] / denom],
            [mass_fracs[0] / denom, 1.0, 0, 0],
            [mass_fracs[0] / denom, 1.0, 0, 0]
            ])
        assert np.allclose(correct, matrix, rtol=1e-3)

    def test_pe_artificial(self):
        """Test using three species artificial model with PE reactions from 2006 DRG paper.
        # Confirmed by hand.
        """
        R1 = ct.Reaction.fromCti('''reaction('F <=> R', [1.0e3, 0.0, 0.0])''')
        R2 = ct.Reaction.fromCti('''reaction('R <=> P', [1.0, 0.0, 0.0])''')

        F = ct.Species('F', 'H:1')
        R = ct.Species('R', 'H:1')
        P = ct.Species('P', 'H:1')

        for sp in [F, R, P]:
            sp.thermo = ct.ConstantCp(
                300, 1000, 101325, (300, 1.0, 1.0, 1.0)
                )
        model = ct.Solution(
            thermo='IdealGas', kinetics='GasKinetics',
            species=[F, R, P], reactions=[R1, R2]
            )
        conc_R = 0.1
        conc_F = ((1 + 1e-3)*conc_R - (1/2e3))/(1 - (1/2e3))
        conc_P = 1.0 - (conc_R + conc_F)
        state = 1000, ct.one_atm, [conc_F, conc_R, conc_P]
        matrix = create_pfa_matrix(state, model)

        correct = np.array([
            [0, 1.0, 1],
            [.5, 0, 1],
            [.5, 1.0, 0],
            ])
        assert np.allclose(correct, matrix, rtol=1e-3)
    
    def test_dormant_modes(self):
        """Test using three species artificial model with dormant modes from 2006 DRG paper.
        # Confirmed by hand. 
        """
        R1 = ct.Reaction.fromCti('''reaction('A <=> B', [1.0, 0.0, 0.0])''')
        R2 = ct.Reaction.fromCti('''reaction('B <=> C', [1.0e-3, 0.0, 0.0])''')

        A = ct.Species('A', 'H:1')
        B = ct.Species('B', 'H:1')
        C = ct.Species('C', 'H:1')

        for sp in [A, B, C]:
            sp.thermo = ct.ConstantCp(
                300, 1000, 101325, (300, 1.0, 1.0, 1.0)
                )
        model = ct.Solution(
            thermo='IdealGas', kinetics='GasKinetics',
            species=[A, B, C], reactions=[R1, R2]
            )
        state = 1000, ct.one_atm, [1.0, 2.0, 1.0]
        matrix = create_pfa_matrix(state, model)

        correct = np.array([
            [0, 1.0, 0],
            [.999000999, 0, .000999000999],
            [0, 1.0, 0],
            ])
        assert np.allclose(correct, matrix, rtol=1e-3)

        conc_A = 1.370536
        conc_B = 1.370480
        conc_C = 1.258985
        state = 1000, ct.one_atm, [conc_A, conc_B, conc_C]
        matrix = create_pfa_matrix(state, model)

        correct = np.array([
            [0, 1.0, 1],
            [.50226468, 0, 1],
            [.50226468, 1, 0],
            ])
        assert np.allclose(correct, matrix, rtol=1e-3)


class TestGraphSearch:
    """Tests for graph_search method"""
    #generate test graph
    #starting from A, nodes A,E,C,F,D,I,H,O should be the only nodes found
    def testGraphSearchOneInput(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(
            ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )

        graph.add_weighted_edges_from([
            #('A','F', 0), ('A','N',0), 
            # ('C','F',1.0), ('A','C',1.0),
            ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

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
        graph.add_nodes_from(
            ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )

        graph.add_weighted_edges_from([
            #('A','F', 0), ('A','N',0),
            # ('C','F',1.0), ('A','C',1.0),
            ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

        subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])
        
        #temporary solution
        essential_nodes = graph_search(subgraph, 'G')
    
        assert 'G' in essential_nodes
        for n in ['A','B', 'C', 'D', 'J', 'K', 'I', 'O', 'F', 'E', 'H', 'M', 'N']:
            assert n not in essential_nodes
        assert [n in essential_nodes for n in [ 'G', 'L']]

    def testGraphSearch3Inputs(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(
            ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )

        graph.add_weighted_edges_from(
            [ ('C','F', 1), ('A','C', 1),
            #('A','F', 0), ('A','N', 0), 
             ('N','C', 1), ('C','D', 1),
             ('D','I', 1), ('I','O', 1), ('A','E', 1),
             #('E','G', 0), ('G','I', 0), ('G','M', 0),
             ('G','L', 1), ('E','H', 1), 
             #('H','J', 0)
             ])

        target_species= ['A', 'C', 'D']

        essential_nodes = graph_search(graph, target_species)
    
        assert 'A' in essential_nodes
        assert 'C' in essential_nodes
        assert 'D' in essential_nodes
        for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']:
            assert n in essential_nodes
        for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']:
            assert n not in essential_nodes
    
    def testgraphsearch_no_targets (self):
        graph = nx.DiGraph()
        graph.add_nodes_from(
            ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )

        graph.add_weighted_edges_from([
            #('A','F', 0), ('A','N',0),
            ('C','F',1.0), ('A','C',1.0),
            ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

        essential_nodes = graph_search(graph, [])
        assert not essential_nodes

    @pytest.mark.xfail
    def testGraphshearchwithatargetThatsnotinGraph(self):	
        graph = nx.DiGraph()
        graph.add_nodes_from(
            ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )

        graph.add_weighted_edges_from([
            #('A','F', 0), ('A','N',0),
            ('C','F',1.0), ('A','C',1.0),
            ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

        essential_nodes = graph_search(graph, 'Z')
        assert 'Z' in essential_nodes

    def testGraphsearchforinfinteloops(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(['A', 'B', 'C', 'D', 'E'])
        
        graph.add_weighted_edges_from(
            [('A', 'B', 1), ('B', 'C', 1), ('C', 'D', 1), ('D', 'E',1), ('E', 'A', 1)]
            )
        
        essential_nodes= graph_search(graph, 'A')
        assert 'A' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'B', 'E']]
        
    @pytest.mark.xfail
    def testGraphShearchWithATargetThatsNotInGraphAndOneThatIs(self):	
        graph = nx.DiGraph()
        graph.add_nodes_from(
            ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )

        graph.add_weighted_edges_from([
            #('A','F', 0), ('A','N',0),
            ('C','F',1.0), ('A','C',1.0),
            ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

        essential_nodes = graph_search(graph, ['B', 'Z'])
        assert 'B' in essential_nodes

    def testGraphsearchwithListofLength1(self):
        graph = nx.DiGraph()
        graph.add_node('A')


        essential_nodes = graph_search(graph, 'A')
        assert 'A' in essential_nodes
        assert len(essential_nodes) == 1

    def testGraphSearchWithTwoOfTheSameItemInTheGraph(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(
            ['A', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )
    
        graph.add_weighted_edges_from([
            #('A','F',0), ('A','N',0),
            ('C','F',1.0), ('A','C',1.0), 
            ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

        essential_nodes = graph_search(graph, 'A')
        assert 'A' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in essential_nodes for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]

    def testGraphSearchWithTwoOfTheSameItemInTheTargetList(self):
        graph = nx.DiGraph()
        graph.add_nodes_from(
            ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O']
            )
    
        graph.add_weighted_edges_from([
            #('A','F', 0), ('A','N',0),
            ('C','F',1.0), ('A','C',1.0),
            ('N','C',1.0), ('C','D',1.0),
            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
            #('E','G',0), ('G','I',0), ('G','M',0),
            ('G','L',1.0), ('E','H',1.0), 
            #('H','J',0)
            ])

        essential_nodes = graph_search(graph, ['A','A'])
        assert 'A' in essential_nodes
        assert [n in essential_nodes for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
        assert [n not in essential_nodes for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]


class TestReducePFA:
    def test_gri_reduction_multiple_thresholds(self):
        """Tests reduce_pfa method with multiple thresholds"""
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

        data = np.genfromtxt(
            relative_location(os.path.join('assets', 'example_ignition_data.dat')), 
            delimiter=','
            )

        model = ct.Solution(model_file)
        matrices = []
        for state in data:
            matrices.append(create_pfa_matrix((state[0], state[1], state[2:]), model))
        
        with TemporaryDirectory() as temp_dir:
            reduced_model = reduce_pfa(
                model_file, ['CH4', 'O2'], ['N2'], 0.14, matrices, 
                np.array([1.066766136745876281e+00, 4.334773545084597696e-02]), ignition_conditions=ignition_conditions,
                previous_model=None, threshold_upper=None, num_threads=1, path=temp_dir
                )
        
        expected_species = [
            'H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'C', 'CH', 'CH2', 'CH2(S)', 
            'CH3', 'CH4', 'CO', 'CO2', 'HCO', 'CH2O', 'CH2OH', 'CH3O', 'CH3OH', 'C2H2', 'C2H3',
            'C2H4', 'C2H5', 'C2H6', 'HCCO', 'CH2CO', 'N', 'NH', 'NH2', 'NNH', 'NO', 'N2O',
            'HNO', 'CN', 'HCN', 'H2CN', 'HCNN', 'HCNO', 'HOCN', 'HNCO', 'NCO', 'N2', 'CH2CHO'
            ]
        assert check_equal(reduced_model.model.species_names, expected_species)
        assert reduced_model.model.n_reactions == 281
        assert round(reduced_model.error, 2) == .14


class TestRunPFA:
    def test_gri_reduction(self):
        """Tests driver run_pfa method"""
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

        # Run PFA
        with TemporaryDirectory() as temp_dir:
            reduced_model = run_pfa(
                model_file, error,  ['CH4', 'O2'], ['N2'], ignition_conditions=ignition_conditions,
                num_threads=1, path=temp_dir
                )

        # Expected answer
        expected_model = ct.Solution(relative_location(os.path.join('assets', 'pfa_gri30.cti')))
        
        # Make sure models are the same
        assert check_equal(reduced_model.model.species_names, expected_model.species_names)
        assert reduced_model.model.n_reactions == expected_model.n_reactions
        assert round(reduced_model.error, 2) == 3.64
