import networkx as nx
import numpy as np
import h5py
from collections import Counter
import time as tm
from graph_search import graph_search

def trim_drg(total_edge_data, solution_object, threshold_value, keeper_list, done, target_species):
    """ Use the Direct Relation Graph (DRG) method to build a nodal graph of
        species and their edge-weights above a certain threshold value

    Parameters
    ----------
    total_edge_data: array
        A 3D array that contains information from the simulation used to calculate DRG interaction coefficients.  
    solution_object : obj
        Cantera Solution object
    threshold_value: int
        the threshold value at this point in the iteraton.  
    keeper_list: array
        an array of strings that represent the sames of the species to keep in the mechanism no matter what.  
    done: array
         A singleton boolean value that represents if the iteration should be finished or not.  
    target_species: list of strings
	a list containing the names of the target species.

    Returns
    -------
    exlusion_list: array
        An array of strings that represent the names of all of the species that should not be inlcuded in the final reduced model for the given threshold value.  
    """
    start_time = tm.time()
    
    #initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    graph = nx.DiGraph() #Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.

    
    safe = [] #A list of species that are to be retained for this threshold value
    
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
                            if weight > old_weight and weight <= 1 and weight > threshold_value: #Only include the weight if it is greater than the threshold value.   
                                graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
			    elif weight > 1:
			    	print "Error.  Edge weights should not be greater than one."
			    	exit()
                        elif weight <= 1 and weight > threshold_value:
                            graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
		        elif weight > 1:
		    		print "Error.  Edge weights should not be greater than one."
		     		exit()
                except IndexError:
                    print edge
                    continue
        
      	    dic = graph_search(graph, target_species) #Search graph for max values to each species based on targets
            for sp in dic: #Add to max dictionary if it is new or greater than the value already there. 
                if sp not in safe:
			safe.append(sp)
            graph.clear() #Reset graph

    core_species = []
    species_objects = solution_object.species()
    
    #Take all species that are over the threshold value and add them to essentail species. 
    essential_species = [] 
    for sp in species_objects:
        if sp.name in safe:
            if sp not in essential_species:
                essential_species.append(sp)
    done[0] = False
   
    #Add all species in essential species to core species 
    for sp in essential_species:
        if sp not in core_species:
            core_species.append(sp.name)

    #Add all of the must keep species to core species
    retained_species = keeper_list #Specified by the user.  A list of species that also need to be kept.
    for sp in retained_species:
        if sp not in core_species:
            core_species.append(sp)

    exclusion_list = []
    
    #Exclude everythong not in core species.
    for species in solution_object.species():
        if species.name not in core_species: #If its not one of our species we must keep, add it to the list of species to be trimmed. 
            exclusion_list.append(species.name)

    return exclusion_list
