import networkx as nx
import cantera as ct
import h5py
from dijkstra import ss_dijkstra_path_length_modified

def graph_search_drgep(nx_graph, target_species, threshold, species,done):
    """Search nodal graph and generate list of species to exclude

    Parameters
    ----------
    nx_graph : obj
        networkx graph object of solution\
    target_species : list
        List of target species to search from
    threshold: float
	The threshold for when to cut species out of the solution.
    species: list
	A list of the species in the model.
    done: singleton
	A singleton boolean value that represents wether the graph has more species on it that have not been cut out or not.  

    Returns
    -------
    essential_nodes : str
        String containing names of essential species
    """

    essential_nodes = []
    for target in target_species:
        dic = ss_dijkstra_path_length_modified(nx_graph, target) #Get dictionary of each values maximum path
        for sp in species: 
	    if sp.name in dic: #For all species in the dictionary, if the greatest path to them is greater than the threshold, add them to the essential species list.
                if dic[sp.name] > threshold and sp not in essential_nodes:
                    essential_nodes.append(sp)
	done[0] = True
	for sp in species:
	    if sp.name in dic: #If more could be cut out in the future, set done to false.  
      	        if dic[sp.name] > threshold:
	            done[0] = False
    return essential_nodes
