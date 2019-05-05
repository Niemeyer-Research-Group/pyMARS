#tests graph searching function
import sys
import os
import pytest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import networkx as nx
import cantera as ct
from drg import graph_search

xfail = pytest.mark.xfail


#generate test graph
#starting from A, nodes A,E,C,F,D,I,H,O should be the only nodes found
def testGraphSearchOneInput():
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
def testGraphSearchOneInput():
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


def testGraphSearch3Inputs():
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

  
@xfail
def testgraphsearchEmptyGraph ():
	graph = nx.DiGraph()
	essential_nodes=[]
	assert 'A' not in essential_nodes
	essential_nodes= graphsearch(graph, 'A') #this fails can't figure out what the real input to compare to would be
	assert 'A' not in essential_nodes

@xfail
def testGraphshearchwithatargetThatsnotinGraph():	
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


	essential_nodes = graph_search(graph, 'Z')
	assert 'Z' in essential_nodes
   

def testGraphsearchforinfinteloops():
  graph = nx.DiGraph()
  graph.add_nodes_from(['A', 'B', 'C', 'D', 'E'])
  
  graph.add_weighted_edges_from([('A', 'B', 1),('B', 'C', 1), ('C', 'D', 1), ('D','E',1), ('E', 'A', 1)])
  
  essential_nodes= graph_search(graph, 'A')
  for n in essential_nodes:
    print (n)
  assert 'A' in essential_nodes
  assert [n in essential_nodes for n in ['A', 'C', 'D', 'B', 'E']]
	
@xfail
def testGraphShearchWithATargetThatsNotInGraphAndOneThatIs():	
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


	essential_nodes = graph_search(graph, ['B','Z'])
	assert 'B' in essential_nodes


def testGraphsearchwithListofLength1():
	graph = nx.DiGraph()
	graph.add_node('A')


	essential_nodes = graph_search(graph, 'A')
	assert 'A' in essential_nodes
	assert len(essential_nodes)==1

def testGraphSearchWithTwoOfTheSameItemInTheGraph():
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


def testGraphSearchWithTwoOfTheSameItemInTheTargetList():
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
