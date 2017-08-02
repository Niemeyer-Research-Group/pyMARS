import networkx as nx
import cantera as ct
import h5py
from dijkstra import ss_dijkstra_path_length_modified

def graph_search_drgep(nx_graph, target_species):
    """Search nodal graph and generate a dictionary of the greatest paths to all species from one of the targets.  

    Parameters
    ----------
    nx_graph : obj
        networkx graph object of solution\
    target_species : list
        List of target species to search from
    
    Returns
    -------
    max_dic : dictionary
        Values of the greatest possible path to each species from one of the targets on this graph keyed by species name.  
    """

    max_dic = {} #A dictionary holding the maximum path to each species.  
    for target in target_species: 
        dic = ss_dijkstra_path_length_modified(nx_graph, target) #Get dictionary of each values maximum path
        for sp in dic:    #If the species are not in the max dictionary or the new value for that species is greater than the one in the max dictionary, add it to the max dictionary.  
	    if sp not in max_dic:                                                                                                                                         
                max_dic[sp] = dic[sp]
	    elif max_dic[sp] < dic[sp]:
                max_dic[sp] = dic[sp]
    return max_dic
