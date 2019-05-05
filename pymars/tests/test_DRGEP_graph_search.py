#test_DRGEP_s graph searching function
import sys
import os
import pytest
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")

import networkx as nx
import cantera as ct
from drgep import graph_search_drgep

xfail = pytest.mark.xfail

#generate test_DRGEP_ graph
#starting from A, nodes A,E,C,F,D,I,H,O should be the only nodes found
def test_DRGEP_GraphSearchSimple():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','B', 1), ('B','C',1), ('C','D',1),('D','E',1),
       	                     ('E','F',1)])

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


#generate test_DRGEP_ graph
#starting from A, nodes A,-F should be found doubling each time
def test_DRGEP_GraphSearchDouble():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','B', 2), ('B','C',2), ('C','D',2),('D','E',2),
       	                     ('E','F',2)])

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


#generate test_DRGEP_ graph
#starting from A, nodes A,-F should be found decreasing by half each time
def test_DRGEP_GraphSearchHalve():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','B', .5), ('B','C',.5), ('C','D',.5),('D','E',.5),
       	                     ('E','F',.5)])


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


#generate test_DRGEP_ graph
#starting from G, nodes G and L should be the only nodes found
def test_DRGEP_GraphSearchFromG():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])

	subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])


	#temporary solution
	max_dic = graph_search_drgep(subgraph, 'G')

	assert 'G' in max_dic
	assert max_dic['G']
	assert [n in max_dic for n in ['G', 'L']]
	assert [n not in max_dic for n in ['A','B','C','D','E','F','H','I', 'J', 'K', 'L','M', 'N']]
	assert max_dic['G']==1
	assert max_dic['L']==1

def test_DRGEP_GraphSearchComplicated():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G'])
	graph.add_weighted_edges_from([('A','B', .1), ('C','D',.2), ('B','C',.3), ('D','E',.6), ('E','F',.5), ('F','A',.4), ('A','E',.7), ('C','G',.01), ('G','A',.01), ('E','B',.05)])


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

def test_DRGEP_GraphSearch3Inputs():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


	target_species= ['A','C','D']

	max_dic = graph_search_drgep(graph, target_species)
  
	assert 'A' in max_dic
	assert 'C' in max_dic
	assert 'D' in max_dic
	assert [n in max_dic for n in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I','J','K','L','M','N','O',]]
	assert [n not in max_dic for n in ['B']]

def test_DRGEP_GraphSearchComplicated2Inputs():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G'])
	graph.add_weighted_edges_from([('A','B', .5),('G','C',.9),('C','D',.5),('D','E',.5),('E','F',.5),('B','F',0.0002),('A','C',.0225)])


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
@xfail
def test_DRGEP_graphsearchEmptyGraph ():
	graph = nx.DiGraph()
	max_dic=[]
	assert 'A' not in max_dic
	max_dic= graphsearch(graph, 'A') #this fails can't figure out what the real input to compare to would be
	assert 'A' not in max_dic

@xfail
def test_DRGEP_GraphshearchwithatargetThatsnotinGraph():	
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


	max_dic = graph_search_drgep(graph, 'Z')
	assert 'Z' in max_dic
   

def test_DRGEP_Graphsearchforinfinteloops():
  graph = nx.DiGraph()
  graph.add_nodes_from(['A', 'B', 'C', 'D', 'E'])
  
  graph.add_weighted_edges_from([('A', 'B', 1),('B', 'C', 1), ('C', 'D', 1), ('D','E',1), ('E', 'A', 1)])
  
  max_dic= graph_search_drgep(graph, 'A')
  for n in max_dic:
    print (n)
  assert 'A' in max_dic
  assert [n in max_dic for n in ['A', 'C', 'D', 'B', 'E']]
	
@xfail
def test_DRGEP_GraphShearchWithATargetThatsNotInGraphAndOneThatIs():	
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])


	max_dic = graph_search_drgep(graph, ['B','Z'])
	assert 'B' in max_dic


def test_DRGEP_GraphsearchwithListofLength1():
	graph = nx.DiGraph()
	graph.add_node('A')


	max_dic = graph_search_drgep(graph, 'A')
	assert 'A' in max_dic
	assert len(max_dic)==1
	assert max_dic['A']==1

def test_DRGEP_GraphSearchWithTwoOfTheSameNodeInTheGraph():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])
  
	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])

	max_dic = graph_search_drgep(graph, 'A')
	assert 'A' in max_dic
	assert [n in max_dic for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
	assert [n not in max_dic for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]


def test_DRGEP_GraphSearchWithTwoOfTheSameItemInTheTargetList():
	graph = nx.DiGraph()
	graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])
  
	graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
       	                     ('A','N',0), ('N','C',1.0), ('C','D',1.0),
       	                     ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
       	                     ('E','G',0), ('G','I',0), ('G','M',0),
       	                     ('G','L',1.0), ('E','H',1.0), ('H','J',0)])

	max_dic = graph_search_drgep(graph, ['A','A'])
	assert 'A' in max_dic
	assert [n in max_dic for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
	assert [n not in max_dic for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]
