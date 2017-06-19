import networkx as nx
import numpy as np
import h5py
from collections import Counter
from graph_search import graph_search
import time as tm

def make_graph(solution_object, threshold_value, total_edge_data, target_species):
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

    Returns
    -------
    exclusion_list : list
        List of species to trim from mechanism
    """
    start_time = tm.time()
    #initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    graph = nx.DiGraph()

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
                    weight = abs(float(numerator[edge])/float(denominator[species_a_name]))
                    #only add edge if greater than edge value from previous timesteps
                    if weight >= threshold_value:
                        #print weight
                        if graph.has_edge(species_a_name, species_b_name):
                            old_weight = graph[species_a_name][species_b_name]['weight']
                            if weight > old_weight:
                                graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
                        else:
                            graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
                except IndexError:
                    print edge
                    continue

        #search the compiled graph
        #temporarily input target species
        essential_species = graph_search(graph, target_species)
        for sp in essential_species:
            if sp not in core_species:
                core_species.append(sp)

        #clear the graph for next individual set
        graph.clear()
    #specific to nc7h_16_mech.cti
    #retained_species = ['n2', 'c5h11o2h-2','c5h11o-2','ic4ketit','ic5ketdb','o2c2h4o2h','ch2o2hcho','ch3choohcoch3','ic3h7coc2h5','ic3h6coc2h5','tc3h6coc2h5','ic3h7coc2h4p','ic3h7coc2h4s','ic3h5coc2h5','ac3h4coc2h5','ic3h5coc2h4p','ic3h5coc2h4s','nc6ket26','ar']

    #specific to gri30
    retained_species = ['n2', 'co2', 'h20']
    #should also have targets of CH4 and O2
    for sp in retained_species:
        if sp not in core_species:
            core_species.append(sp)

    #temporarily performs part of the graph search function, by finding
    #dicsonnected species
    exclusion_list = []

    for species in solution.species():
        if species.name not in core_species:
            exclusion_list.append(species.name)

    return exclusion_list
