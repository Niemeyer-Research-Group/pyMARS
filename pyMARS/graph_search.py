import networkx as nx
import cantera as ct
import h5py

def graph_search(solution_object, nx_graph):
    """Search nodal graph and generate list of species to exclude

    :param solution_object:
        Original Cantera solution Object

    :param nx_graph:
        networkx graph object

    :returns exclusion_list:
        String containing names of species to exclude
    """
