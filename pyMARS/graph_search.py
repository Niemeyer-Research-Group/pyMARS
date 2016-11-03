import networkx as nx
import cantera as ct
import h5py

def graph_search(solution_object, nx_graph, target_species):
    """Search nodal graph and generate list of species to exclude

    Parameters
    ----------
    solution_object : obj
        Original Cantera solution Object
    nx_graph : obj
        networkx graph object of solution

    Returns
    -------
    exclusion_list : str
        String containing names of species to exclude
    """

    essential_nodes = list(nx.dfs_preorder_nodes(nx_graph, target_species))

    exclusion_list = []
    for species in solution_object.species():
        if species.name not in essential_nodes:
            exclusion_list.append(species.name)
    return exclusion_list
