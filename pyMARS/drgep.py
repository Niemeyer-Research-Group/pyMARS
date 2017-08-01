import networkx as nx
import numpy as np
import h5py
from collections import Counter
from graph_search_drgep import graph_search_drgep
import time as tm

def make_graph_drgep(solution_object, threshold_value, total_edge_data, target_species, done):
    """ Use the Direct Relation Graph (DRG) method to build a nodal graph of
        species and their edge-weights above a certain threshold value

    Parameters
    ----------
    solution_object : obj
        Cantera Solution object
    hdf5_file : str
        data file containing individual reaction production rates
    threshold_value : int
        an edge weight threshold value
    target_species: list of strings
	a list containing the names of the target species.
    done: singleton
	a singleton boolean value that represents if their are more species to cut out or not.  

    Returns
    -------
    graph: NetworkX graph
        A graph with edge weights representing species dependencies on each other.
    """
    start_time = tm.time()
    #initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    graph = nx.DiGraph() #Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.

    #add nodes to graph
    for species in species_objects:
        graph.add_node(species.name)
    ri_total = {}
    ri_partial = {}
    error_list = {}
    core_species = []

    #calculate edge weights based on list received from get_rate_data
    #initial condition
    for ic in total_edge_data.iterkeys():
        for species in species_objects:
            graph.add_node(species.name)
        #timestep
        for tstep in total_edge_data[ic].iterkeys():
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
    return graph
 

def run_drgep(graph,solution_object, threshold_value, target_species, keeper_list, done):
    """ Use the Direct Relation Graph (DRG) method to build a nodal graph of
        species and their edge-weights above a certain threshold value

    Parameters
    ----------
    graph: NetworkX graph
        A graph with edge weights representing species dependencies on each other.
    solution_object : obj
        Cantera Solution object
    threshold_value : int
        an edge weight threshold value
    target_species: list of strings
	a list containing the names of the target species.
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
    #search the compiled graph with a modified version of dikstras algorithm to determine what species need to stay in the model at this threshold value.
    essential_species = graph_search_drgep(graph, target_species, threshold_value, species_objects, done)
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
