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
    #iterate through each initial condition set
    for initial_condition, data_group in rate_file.iteritems():
        #generate dict of sum production Rates, and make zero before evaluating
        #each initial condition set
        for species in species_objects:
            ri_total[species.name] = 0.0
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
            #print original_sp_npr['nc7h16'].value
            """
            if first_iteration:
                #generate dict of sum production Rates
                for species in species_objects:
                    ri_total[species.name] = 0
                for pre_defined_edge in ri_partial:
                    try:
                        ri_partial[str(pre_defined_edge)] = 0.0
                    except Exception:
                        continue
                first_iteration = False
            """
            #build weights by iterating over every reaction
            for reaction_number, reaction in enumerate(reaction_objects):
                reaction_production_rate = float(rxn_prod_rates[reaction_number])
                reactants = reaction.reactants
                products = reaction.products
                #A = Counter(reactants)
                #B = Counter(products)
                #for key in A:
                #    A[key] *= -1
                #all_species = A + B
                #generate list of all species
                all_species = reaction.reactants
                all_species.update(reaction.products)
                #for species_a in reactants:
                #    if species_a in products:
                #        if reactants[species_a] != products[species_a]:
                #            error_list[reaction] = [reaction_number, reactants[species_a], products[species_a]]
                if reaction_production_rate != 0:
                    #for species_a in all_species:
                    for species_a in products:
                        molar_coeff_A = float(products[species_a])
                        #molar_coeff_A = float(reactants[species_a])
                        #this is denominator
                        ri_total[species_a] += abs(reaction_production_rate* molar_coeff_A)
                        for species_b in all_species:
                        #for species_b in products:
                            partial_name = species_a + '_' + species_b
                            if species_a == species_b:
                                continue
                            #this is numerator
                            try:
                                ri_partial[partial_name] += abs((reaction_production_rate*molar_coeff_A))
                            except KeyError:
                                ri_partial[partial_name] = abs((reaction_production_rate*molar_coeff_A))
                    #alternate method for calculating total production rate
                    if reaction_production_rate > 0:
                        for species in products:
                            ri_total[species] += float(reaction_production_rate*products[species])
                        for species in reactants:
                            ri_total[species] += float(-reaction_production_rate*reactants[species])
                    if reaction_production_rate < 0:
                            for species in products:
                                ri_total[species] += float(reaction_production_rate*products[species])
                            for species in reactants:
                                ri_total[species] += float(reaction_production_rate*reactants[species])


            #check to make sure calculated net production rate is correct
            # this is numerator
            for value in ri_total:
                if abs(ri_total[value] - original_sp_npr[value].value) > .01:
                    #print ('species: %s and amount %0.5f')  %(value, abs(ri_total[value] - original_sp_npr[value].value))
                    continue
                ri_total[value] = original_sp_npr[value].value
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
