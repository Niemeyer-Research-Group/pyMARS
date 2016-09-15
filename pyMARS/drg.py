import networkx as nx
import cantera as ct
import numpy as np
import h5py
import matplotlib.pyplot as plt
import time

def make_graph(solution_object, hdf5_file, threshold_value, target_species):
    """ Use the Direct Relation Graph (DRG) method to choose species to
    eliminate from reaction mechanism.

    Parameters:
    -------------
    data_file: Mechanism file in the form of .cti
    hdf5_file: data file containing individual reaction production rates

    Returns:
    -------------
    exclusion_list: List of species to trim from mechanism
    """

    f = h5py.File(hdf5_file, 'r')

    #initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    #ask for threshold value
    threshold = threshold_value

    G = nx.MultiGraph()
    node_labels = {}
    count = 0
    start_time = time.time()

    #add nodes to graph
    for ind, sp in enumerate(species_objects):
        G.add_node(sp.name)
    #iterate through each timestep
    for nm, grp in f.iteritems():
        count += 1
        rxn_prod_rates = np.array(grp['Reaction Production Rates'])
        #generate dict of sum production Rates
        ri_total = {}
        ri_partial = {}
        for k in species_objects:
            ri_total[k.name] = 0
        #iterate through reactions and build weights
        for i, reac in enumerate(reaction_objects):
            products = reac.products
            product_names = reac.products.keys()
            reactants = reac.reactants
            reactant_names = reac.reactants.keys()
            #generate list of all species
            all_species = reac.products
            all_species.update(reac.reactants)
            reaction_production_rate = float(rxn_prod_rates[i])
            if reaction_production_rate != 0:
                for spx in all_species:
                    species_A = spx
                    molar_coeff_A = float(all_species[species_A])
                    ri_total[species_A] += abs(reaction_production_rate* molar_coeff_A)
                    for spy in all_species:
                        species_B = spy
                        partial_name = species_A + '_' + species_B
                        if spy == spx:
                            continue
                        try:
                            ri_partial[partial_name] += abs((reaction_production_rate*molar_coeff_A))
                        except KeyError:
                            ri_partial[partial_name] = abs((reaction_production_rate*molar_coeff_A))
        #clean up duplicate edges
        ri_partial_temp = dict(ri_partial)
        for ik in ri_partial:
            value_one = ri_partial[ik]
            split_name = ik.split('_', 1)
            switched_name = str(split_name[1]) + '_' + str(split_name[0])
            for il in ri_partial:
                value_two = ri_partial[il]
                if il == switched_name:
                    ri_partial_temp[ik] = value_one + value_two
                    del ri_partial_temp[il]
                else:
                    continue
        ri_partial = dict(ri_partial_temp)

        #divide progress related to species B by total progress
        #and make edge
        for ind in ri_partial:
            try:
                both = ind.split('_', 1)
                sp_A = both[0]
                sp_B = both[1]
                weight = abs(float(ri_partial[ind])/float(ri_total[sp_A]))
                #only add edge if > than edge value from previous timesteps
                if weight >= threshold:
                    if G.has_edge(sp_A, sp_B):
                        old_weight = G[sp_A][sp_B][0]['weight']
                        if weight > old_weight:
                            G.add_edge(sp_A, sp_B, weight=weight)
                    else:
                        G.add_edge(sp_A, sp_B, weight=weight)
            except IndexError:
                print ind
                continue
        progress =  str(count) + '/40 timesteps'
        #print progress
        #print("--- %s seconds ---" % (time.time() - start_time))
        #print (nm, G.number_of_edges())

    #get connected species
    target = target_species
    essential_nodes = list(nx.node_connected_component(G, target))

    #nx.draw(G, with_labels=True, width=.25)

    #get list of species to eliminate
    exclusion_list = list()
    ex_list = []
    for spec in solution.species():
        ind_name = spec.name
        if ind_name not in essential_nodes:
            exclusion_list.append(spec.name)
            ex_list.append(spec.name)
    exclusion_list_string = '\''
    for spc in exclusion_list:
        exclusion_list_string += spc + ', '
    exclusion_list_string = exclusion_list_string.rstrip(',')
    plt.show(block=False)
    print G.number_of_edges()
    print G.number_of_nodes()
    f.close()
    return ex_list

#make_graph('gri301.cti', 'production_rates.hdf5')
