import networkx as nx
import numpy as np
import h5py
import matplotlib.pyplot as plt
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
    #iterate through each timestep
    for timestep, data_group in rate_file.iteritems():
        rxn_prod_rates = np.array(data_group['Reaction Production Rates'])
        #generate dict of sum production Rates
        for species in species_objects:
            ri_total[species.name] = 0
        for pre_defined_edge in ri_partial:
            try:
                ri_partial[str(pre_defined_edge)] = 0.0
            except Exception:
                continue
        #build weights
        for reaction_number, reaction in enumerate(reaction_objects):
            reaction_production_rate = float(rxn_prod_rates[reaction_number])
            products = reaction.products
            reactants = reaction.reactants
            #generate list of all species
            all_species = reaction.products
            all_species.update(reaction.reactants)
            for species_a in reactants:
                if species_a in products:
                    if reactants[species_a] != products[species_a]:
                        error_list[reaction] = [reaction_number, reactants[species_a], products[species_a]]
            if reaction_production_rate != 0:
                for species_a in all_species:
                    molar_coeff_A = float(all_species[species_a])
                    #this is denominator
                    ri_total[species_a] += abs(reaction_production_rate* molar_coeff_A)
                    for species_b in all_species:
                        partial_name = species_a + '_' + species_b
                        if species_a == species_b:
                            continue
                        #this is numerator
                        try:
                            ri_partial[partial_name] += abs((reaction_production_rate*molar_coeff_A))
                        except KeyError:
                            ri_partial[partial_name] = abs((reaction_production_rate*molar_coeff_A))

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
                    if graph.has_edge(species_a_name, species_b_name):
                        old_weight = graph[species_a_name][species_b_name]['weight']
                        if weight > old_weight:
                            graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
                    else:
                        graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
            except IndexError:
                print edge
                continue

    rate_file.close()

    if len(error_list) != 0:
        print 'error list'
        print error_list
    return graph
