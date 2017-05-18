import networkx as nx
import cantera as ct
import h5py

def graph_search(nx_graph, target_species):
    """Search nodal graph and generate list of species to exclude

    Parameters
    ----------
    nx_graph : obj
        networkx graph object of solution\
    target_species : list
        List of target species to search from

    Returns
    -------
    essential_nodes : str
        String containing names of essential species
    """

    if len(target_species) > 1:
        essential_nodes = list()
        for target in target_species:
            essential = list(nx.dfs_preorder_nodes(nx_graph, target))
            for sp in essential:
                if sp not in essential_nodes:
                    essential_nodes.append(sp)
    else:
        essential_nodes = list(nx.dfs_preorder_nodes(nx_graph, target_species[0]))

    return essential_nodes
