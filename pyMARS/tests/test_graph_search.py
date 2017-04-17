#tests graph searching function
import networkx as nx
from pyMARS import graph_search





#generate test graph
#starting from A, nodes E,C,F,D,I,H,O should be the only nodes found
graph = nx.DiGraph()
graph.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])

graph.add_weighted_edges_from([('A','F', 0), ('C','F',1.0), ('A','C',1.0),
                            ('A','N',0), ('N','C',1.0), ('C','D',1.0),
                            ('D','I',1.0), ('I','O',1.0), ('A','E',1.0),
                            ('E','G',0), ('G','I',0), ('G','M',0),
                            ('G','L',1.0), ('E','H',1.0), ('H','J',1.0)])
