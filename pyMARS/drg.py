import networkx as nx
import cantera as ct
import numpy as np
import h5py
import matplotlib.pyplot as plt

def make_graph(solution_object, hdf5_file, threshold_value, target_species):
    """ Use the Direct Relation Graph (DRG) method to choose species to
    eliminate from reaction mechanism.

    :param data_file:
        Mechanism file in the form of .cti
    :param hdf5_file:
        data file containing individual reaction production rates
    :param threshold_value:
        an edge weight threshold value
    :param target_species:
        species to start graph search from (usually the fuel)

    :return exclusion_list:
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
    error_list ={}
    #iterate through each timestep
    for timestep, data_group in rate_file.iteritems():
        rxn_prod_rates = np.array(data_group['Reaction Production Rates'])
        #generate dict of sum production Rates
        for species in species_objects:
            ri_total[species.name] = 0
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
                weight = abs(float(ri_partial[edge])/float(ri_total[species_a_name]))
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

    #get connected species
    target = target_species
    essential_nodes = list(nx.dfs_preorder_nodes(graph, target))
    nx.draw(graph, with_labels=True, width=.25)
    #get list of species to eliminate
    exclusion_list = list()
    ex_list = []
    for species in solution.species():
        if species.name not in essential_nodes:
            exclusion_list.append(species.name)
            ex_list.append(species.name)
    exclusion_list_string = '\''
    for species in exclusion_list:
        exclusion_list_string += str(species) + ', '
    exclusion_list_string = exclusion_list_string.rstrip(',')
    edges = nx.get_edge_attributes(graph, 'weight')
    for edge in edges.items():
        species_names = str(edge[0])
        weight = edge[1]
        edge_string = species_names + '   ' + str(weight) + '\n'
        if weight < .1:
            print edge_string
    #plt.show()
    rate_file.close()
    if len(error_list) != 0:
        print 'error list'
        print error_list
    return ex_list

#make_graph('gri301.cti', 'production_rates.hdf5')
