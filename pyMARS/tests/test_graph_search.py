#tests graph searching function
import networkx as nx
from pyMARS import graph_search
import cantera as ct




#generate test graph
#starting from A, nodes A,E,C,F,D,I,H,O should be the only nodes found
graph = nx.DiGraph()
graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                            ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                            ('E','G',0), ('G','I',0), ('G','M',0),
                            ('G','L',1.0), ('E','H',1.0), ('H','J',0)])
original_edges = graph.edges_iter()

subgraph = nx.DiGraph([(u,v,d) for u,v,d in graph.edges(data=True) if d['weight'] > 0])


#temporary solution
soln = ct.Solution('pym_gri30.cti')
essential_nodes = graph_search.graph_search(soln, subgraph, 'A')

assert 'A' in essential_nodes
assert [n in essential_nodes for n in ['A', 'C', 'D', 'I', 'O', 'F', 'E', 'H']]
assert [n not in essential_nodes for n in ['B', 'G', 'J', 'K', 'L', 'M', 'N']]
