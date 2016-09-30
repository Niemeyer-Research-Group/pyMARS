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

    target_species = str(raw_input('Enter target starting species: '))
    essential_nodes = list(nx.dfs_preorder_nodes(nx_graph, target_species))

    exclusion_list = list()
    ex_list = []

    for species in solution_object.species():
        if species.name not in essential_nodes:
            exclusion_list.append(species.name)
            ex_list.append(species.name)

    exclusion_list_string = '\''
    for species in exclusion_list:
        exclusion_list_string += str(species) + ','
    exclusion_list_string = exclusion_list_string.rstrip(',')
