import networkx as nx
import numpy as np
import h5py
from collections import Counter
from graph_search_drgep import graph_search_drgep
import time as tm

def make_graph_drgep(solution_object, total_edge_data, target_species):
    """ Use the Direct Relation Graph (DRG) method to build a nodal graph of
        species and their edge-weights above a certain threshold value

    Parameters
    ----------
    solution_object : obj
        Cantera Solution object
    hdf5_file : str
        data file containing individual reaction production rates
    target_species: list of strings
	a list containing the names of the target species.

    Returns
    -------
    max_dic: Dictionary
        A Dictionary keyed by species name with values that represent that species importance to the system.  
    """
    start_time = tm.time()
    #initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    graph = nx.DiGraph() #Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.

    #add nodes for all species to graph
    #for species in species_objects: 
    #    graph.add_node(species.name)
    #ri_total = {}
    #ri_partial = {}
    #error_list = {}
    #core_species = []
    
    max_dic = {} #Dictionary holding the maximum values for the iteration
    #calculate edge weights based on list received from get_rate_data
    #initial condition
    for ic in total_edge_data.iterkeys(): #For each initial condition
        for species in species_objects: #Make graph
            graph.add_node(species.name)
        #timestep
        for tstep in total_edge_data[ic].iterkeys(): #Set edge values for the graph
            numerator = total_edge_data[ic][tstep][2]
            denominator = total_edge_data[ic][tstep][1]
            #each species
            for edge in numerator:
                try:
                    edge_name = edge.split('_', 1)
                    species_a_name = edge_name[0]
                    species_b_name = edge_name[1]
                    #dge weight between two species
	            if denominator[species_a_name] != 0:
                        weight = abs(float(numerator[edge])/float(denominator[species_a_name]))
                        if graph.has_edge(species_a_name, species_b_name):
                            old_weight = graph[species_a_name][species_b_name]['weight']
                            if weight > old_weight and weight <= 1:
                                graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
			    elif weight > 1:
			    	print "Error.  Edge weights should not be greater than one."
			    	exit()
                        elif weight <= 1:
                            graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
		        elif weight > 1:
		    		print "Error.  Edge weights should not be greater than one."
		     		exit()
                except IndexError:
                    print edge
                    continue
        
        dic = graph_search_drgep(graph, target_species) #Search graph for max values to each species based on targets
        for sp in dic: #Add to max dictionary if it is new or greater than the value already there. 
            if sp not in max_dic:
                max_dic[sp] = dic[sp]
            elif dic[sp] > max_dic[sp]:
                max_dic[sp] = dic[sp]
        graph.clear() #Reset graph
    return max_dic
 

def run_drgep(max_dic, solution_object, threshold_value, keeper_list, done):
    """ Use the dictionary created by the drgep method to determine what should
        be cut out of the model at a specific threshold value.         

    Parameters
    ----------
    max_dic: Dictionary
        A dictionary keyed by species name that has values representing the species importance to the system.  
    solution_object : obj
        Cantera Solution object
    threshold_value : int
        an edge weight threshold value
    keeper_species: list of strings
	a list of species for the mechinism to keep no matter what
    done: singleton
	a singleton boolean value that represents if their are more species to cut out or not.  

    Returns
    -------
    exclusion_list : list
        List of species to trim from mechanism
    """
        

    core_species = []
    species_objects = solution_object.species()
    
    #Take all species that are over the threshold value and add them to essentail species. 
    essential_species = []
    for sp in species_objects:
        if sp.name in max_dic:
            if max_dic[sp.name] > threshold_value and sp not in essential_species:
                essential_species.append(sp)
    done[0] = True
    for sp in species_objects: #If any more can be taken out, we are not done yet.  
        if sp.name in max_dic: 
            if max_dic[sp.name] > threshold_value:
                done[0] = False
    
    for sp in essential_species:
        if sp not in core_species:
            core_species.append(sp.name)

        #clear the graph for next individual set
        #graph.clear()
	#specific to nc7h_16_mech.cti
        #retained_species = ['n2', 'c5h11o2h-2','c5h11o-2','ic4ketit','ic5ketdb','o2c2h4o2h','ch2o2hcho','ch3choohcoch3','ic3h7coc2h5','ic3h6coc2h5','tc3h6coc2h5','ic3h7coc2h4p','ic3h7coc2h4s','ic3h5coc2h5','ac3h4coc2h5','ic3h5coc2h4p','ic3h5coc2h4s','nc6ket26','ar']
    #specific to gri30 N2 CO2 H2O
    #should also have targets of CH4 and O2
    
    retained_species = keeper_list #Specified by the user.  A list of species that also need to be kept.
    for sp in retained_species:
        if sp not in core_species:
            core_species.append(sp)

    exclusion_list = []

    for species in solution_object.species():
        if species.name not in core_species: #If its not one of our species we must keep, add it to the list of species to be trimmed. 
            exclusion_list.append(species.name)

    return exclusion_list
