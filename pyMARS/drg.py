import networkx as nx
import numpy as np
import h5py
from collections import Counter
from graph_search import graph_search

def make_graph(solution_object, hdf5_file, threshold_value):
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

    rate_file = h5py.File(hdf5_file, 'r')

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
    first_iteration = True
    core_species = []
    ri_total_test = {}
    #iterate through each initial condition set
    for initial_condition, data_group in rate_file.iteritems():
        #generate dict of sum production Rates, and make zero before evaluating
        #each initial condition set
        for species in species_objects:
            ri_total[species.name] = 0.0
            ri_total_test[species.name] = 0.0
        for pre_defined_edge in ri_partial:
            try:
                ri_partial[str(pre_defined_edge)] = 0.0
            except Exception:
                print pre_defined_edge
                continue
        #iterate through every timestep for the initial condition
        for timestep, data_group in rate_file[initial_condition].iteritems():
            rxn_prod_rates = np.array(data_group['Reaction Production Rates'])
            original_sp_npr = rate_file[initial_condition][timestep]['Species Net Production Rates Original']
            for species in species_objects:
                ri_total[species.name] = 0.0
                ri_total_test[species.name] = 0

            #build weights by iterating over every reaction
            for reaction_number, reaction in enumerate(reaction_objects):
                reaction_production_rate = float(rxn_prod_rates[reaction_number])
                reactants = reaction.reactants
                products = reaction.products
                all_species = reactants
                all_species.update(products)


                if reaction_production_rate != 0:

                    if reaction_production_rate > 0:
                        for species_a in products:
                            mcA = float(products[species_a])
                            ri_total[species_a] +=(reaction_production_rate*mcA)
                            for species_b in all_species:
                                if species_a != species_b:
                                    partial_name = species_a + '_' + species_b
                                    try:
                                        ri_partial[partial_name] += (reaction_production_rate*mcA)
                                    except KeyError:
                                        ri_partial[partial_name] = (reaction_production_rate*mcA)
                        for species_a in reactants:
                            mcA = float(reactants[species_a])
                            ri_total[species_a] +=(reaction_production_rate*-mcA)
                            for species_b in all_species:
                                if species_a != species_b:
                                    partial_name = species_a + '_' + species_b
                                    try:
                                        ri_partial[partial_name] += (reaction_production_rate*-mcA)
                                    except KeyError:
                                        ri_partial[partial_name] = (reaction_production_rate*-mcA)

                    if reaction_production_rate < 0:
                        for species_a in products:
                            mcA = float(products[species_a])
                            ri_total[species_a] +=(reaction_production_rate*mcA)
                            for species_b in all_species:
                                if species_a != species_b:
                                    partial_name = species_a + '_' + species_b
                                    try:
                                        ri_partial[partial_name] += (reaction_production_rate*mcA)
                                    except KeyError:
                                        ri_partial[partial_name] = (reaction_production_rate*mcA)
                        for species_a in reactants:
                            mcA = float(reactants[species_a])
                            ri_total[species_a] +=(reaction_production_rate*-mcA)
                            for species_b in all_species:
                                if species_a != species_b:
                                    partial_name = species_a + '_' + species_b
                                    try:
                                        ri_partial[partial_name] += (reaction_production_rate*-mcA)
                                    except KeyError:
                                        ri_partial[partial_name] = (reaction_production_rate*-mcA)
                    #alternate method for calculating total production rate
                    if reaction_production_rate > 0:
                        for species in products:
                            ri_total_test[species] += float(reaction_production_rate*products[species])
                        for species in reactants:
                            ri_total_test[species] += float(-reaction_production_rate*reactants[species])
                    if reaction_production_rate < 0:
                            for species in products:
                                ri_total_test[species] += float(-reaction_production_rate*products[species])
                            for species in reactants:
                                ri_total_test[species] += float(reaction_production_rate*reactants[species])



            #check to make sure calculated net production rate is correct
            # this is numerator
            names_list ={}
            total_error = 0

            #original_sp_npr is a dict written to the h5py file
            #ri_total_test is compiled above starting on line 118
            #this means there is an error in how reaction_production_rate is stored?
            #or how original_sp_npr is stored?
            for sp_name in ri_total_test:
                if abs(ri_total_test[sp_name] - original_sp_npr[sp_name]) > .01:
                    print '-----------------'
                    print 'species: %0.5s, error %0.5f' %(value, float(ri_total_test[value] - names_list[value]))
                    continue
                    total_error += float(ri_total_test[value] - names_list[value])
                ri_total[value] = original_sp_npr[value].value
            print '------------'
            print total_error
            print '------------'
                #print value
                #print ri_total[value]



        #divide progress related to species B by total progress
        #and make edge
        for edge in ri_partial:
            try:
                edge_name = edge.split('_', 1)
                species_a_name = edge_name[0]
                species_b_name = edge_name[1]
                try:
                    weight = abs(float(ri_partial[edge])/float(ri_total[species_a_name]))
                except ZeroDivisionError:
                    if ri_partial > 0.001:
                        continue
                    else:
                        continue
                #only add edge if > than edge value from previous timesteps
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
        essential_species = graph_search(solution, graph, 'nc7h16')
        for sp in essential_species:
            if sp not in core_species:
                core_species.append(sp)

        #clear the graph for next individual set
        graph.clear()


    #temporarily performs part of the graph search function, by finding
    #dicsonnected species
    exclusion_list = []
    for species in solution.species():
        if species.name not in core_species:
            exclusion_list.append(species.name)
    #return exclusion_list

    rate_file.close()

    if len(error_list) != 0:
        print 'error list'
        print error_list
    #return graph
    return exclusion_list
