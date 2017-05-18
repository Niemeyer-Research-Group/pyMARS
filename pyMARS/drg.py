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
    first_iteration = True
    core_species = []
    ri_total_test = {}

    #this code has been moved to the get rate data file
    #iterate through each initial condition set
    # for initial_condition, data_grp in rate_file.iteritems():
    #     #generate dict of sum production Rates, and make zero before evaluating
    #     #each initial condition set
    #     for species in species_objects:
    #         ri_total[species.name] = 0.0
    #         ri_total_test[species.name] = 0.0
    #     for pre_defined_edge in ri_partial:
    #         try:
    #             ri_partial[str(pre_defined_edge)] = 0.0
    #         except Exception:
    #             print pre_defined_edge
    #             continue
    #     #iterate through every timestep for the initial condition
    #     error_dict = {}
    #     for timestep, data_group in data_grp.iteritems():
    #         rxn_prod_rates = np.array(data_group['Reaction Production Rates'])
    #         #original_sp_npr = rate_file[initial_condition][timestep]['Species Net Production Rates Original']
    #         original_sp_npr = data_group['Species Net Production Rates Original']
    #
    #         #bring in reaction coefficient data from get_rate_data, and check that it matches with solution object here
    #         rxn_groups = data_group['Reactions']
    #         ij = 0
    #         for i, rxn_id_group in rxn_groups.iteritems():
    #             #print dir(rxn_id_group)
    #             #print rxn_id_group
    #             #print rxn_id_group['Reactants']
    #             for reac in solution.reaction(ij).reactants:
    #                 assert reac in rxn_id_group['Reactants']
    #             for prod in solution.reaction(ij).products:
    #                 assert prod in rxn_id_group['Products']
    #             ij += 1.0
    #         #the above block shows that coefficient data for species in reactions is appropriately matched
    #
    #
    #         #reset values for calculation of species production
    #         for species in species_objects:
    #             ri_total[species.name] = 0.0
    #             ri_total_test[species.name] = 0.0
    #         for pre_defined_edge in ri_partial:
    #             try:
    #                 ri_partial[str(pre_defined_edge)] = 0.0
    #             except Exception:
    #                 print pre_defined_edge
    #                 continue
    #
    #         #build weights by iterating over every reaction
    #         for reaction_number, reaction in enumerate(reaction_objects):
    #             reaction_production_rate = float(rxn_prod_rates[reaction_number])
    #             reactants = reaction.reactants
    #             products = reaction.products
    #             all_species = reactants
    #             all_species.update(products)
    #             if reaction_production_rate != 0:
    #                 if reaction_production_rate > 0:
    #                     for species_a in products:
    #                         mcA = float(products[species_a])
    #                         ri_total[species_a] += float((reaction_production_rate*mcA))
    #                         for species_b in all_species:
    #                             if species_a != species_b:
    #                                 partial_name = species_a + '_' + species_b
    #                                 try:
    #                                     ri_partial[partial_name] += float((reaction_production_rate*mcA))
    #                                 except KeyError:
    #                                     ri_partial[partial_name] = float((reaction_production_rate*mcA))
    #                                     continue
    #                     for species_a in reactants:
    #                         mcA = float(reactants[species_a])
    #                         ri_total[species_a] += float((reaction_production_rate*-mcA))
    #                         for species_b in all_species:
    #                             if species_a != species_b:
    #                                 partial_name = species_a + '_' + species_b
    #                                 try:
    #                                     ri_partial[partial_name] += float((-reaction_production_rate*mcA))
    #                                 except KeyError:
    #                                     ri_partial[partial_name] = float((-reaction_production_rate*mcA))
    #                                     continue
    #                     #alternate method for calculating total production rate
    #                     for species in products:
    #                         ri_total_test[species] += float(reaction_production_rate*products[species])
    #                     for species in reactants:
    #                         ri_total_test[species] += float(-reaction_production_rate*reactants[species])
    #
    #                 if reaction_production_rate < 0:
    #                     for species_a in products:
    #                         mcA = float(products[species_a])
    #                         ri_total[species_a] += float((reaction_production_rate*mcA))
    #                         for species_b in all_species:
    #                             if species_a != species_b:
    #                                 partial_name = species_a + '_' + species_b
    #                                 try:
    #                                     ri_partial[partial_name] += float((reaction_production_rate*mcA))
    #                                 except KeyError:
    #                                     ri_partial[partial_name] = float((reaction_production_rate*mcA))
    #                                     continue
    #                     for species_a in reactants:
    #                         mcA = float(reactants[species_a])
    #                         ri_total[species_a] +=float(-reaction_production_rate*mcA)
    #                         for species_b in all_species:
    #                             if species_a != species_b:
    #                                 partial_name = species_a + '_' + species_b
    #                                 try:
    #                                     ri_partial[partial_name] += float(-reaction_production_rate*mcA)
    #                                 except KeyError:
    #                                     ri_partial[partial_name] = float(-reaction_production_rate*mcA)
    #                                     continue
    #
    #                     #alternate method for calculating total production rate
    #                     for species in products:
    #                         ri_total_test[species] += float(reaction_production_rate*products[species])
    #                     for species in reactants:
    #                         ri_total_test[species] += float(-reaction_production_rate*reactants[species])
    #
    #
    #         #check to make sure calculated net production rate is correct
    #         # this is numerator
    #         names_list ={}
    #         total_error = 0
    #
    #
    #         #original_sp_npr is a dict written to the h5py file
    #         #ri_total_test is compiled above starting on line 118
    #         #this means there is an error in how reaction_production_rate is stored?
    #         #or how original_sp_npr is stored?
    #         #now know this is not correct, and the issue lies in either in how the value is calculated, or the molar coeff
    #         #right now, all of the computed values are below correct values
    #         error_list = []
    #         error_dict = {}
    #         for sp_name in ri_total_test:
    #             if abs(float(ri_total_test[sp_name]) - float(original_sp_npr[sp_name].value)) > .01:
    #                 if sp_name not in error_dict:
    #                     error_dict[sp_name] = ri_total_test[sp_name] - original_sp_npr[sp_name].value
    #                 else:
    #                     error_dict[sp_name] += ri_total_test[sp_name] - original_sp_npr[sp_name].value
    #                 #print '-----------------'
    #                 #print 'species: %0.5s, error %0.5f' %(sp_name, float(ri_total_test[sp_name] - original_sp_npr[sp_name].value))
    #                 total_error += float(ri_total_test[sp_name] - original_sp_npr[sp_name].value)
    #             ri_total[sp_name] = original_sp_npr[sp_name].value
    #             error_list.append(total_error)
    #         error_list.sort()
    #         print error_dict.keys()
    #
    #         #print '------------'
    #         #print 'Total Error %0.5f' %error_list[0]
    #         #print '------------'
    #             #print value
    #             #print ri_total[value]
    #
    #
    #
    #         #divide progress related to species B by total progress
    #         #and make edge
    #         for edge in ri_partial:
    #             try:
    #                 edge_name = edge.split('_', 1)
    #                 species_a_name = edge_name[0]
    #                 species_b_name = edge_name[1]
    #                 try:
    #                     weight = abs(float(ri_partial[edge])/float(ri_total[species_a_name]))
    #                 except ZeroDivisionError:
    #                     if ri_partial > 0.001:
    #                         continue
    #                     else:
    #                         continue
    #                 #only add edge if > than edge value from previous timesteps
    #                 if weight >= threshold_value:
    #                     #print weight
    #                     if graph.has_edge(species_a_name, species_b_name):
    #                         old_weight = graph[species_a_name][species_b_name]['weight']
    #                         if weight > old_weight:
    #                             graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
    #                     else:
    #                         graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
    #             except IndexError:
    #                 print edge
    #                 continue
    #         #search the compiled graph
    #         #temporarily input target species
    #         essential_species = graph_search(solution, graph, 'nc7h16')
    #         for sp in essential_species:
    #             if sp not in core_species:
    #                 core_species.append(sp)
    #
    #         #clear the graph for next individual set
    #         graph.clear()
    #

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
                    weight = abs(float(numerator[edge])/float(denominator[species_a_name]))
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
    # print '----core species----'
    # print len(core_species)
    # print '--------------------'
    # print '----solution species-----'
    # print len(solution.species())
    # print '--------------------------'
    for species in solution.species():
        if species.name not in core_species:
            exclusion_list.append(species.name)
    #return exclusion_list
    # assert 'nc7h16' not in exclusion_list
    # print '----------exclusion list-------'
    # print len(exclusion_list)
    # print '-------------------------------'
    # print 'drg time: %0.5f' %(tm.time() - start_time)
    return exclusion_list
