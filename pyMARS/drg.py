import networkx as nx
import cantera as ct
import numpy as np
import h5py
import matplotlib.pyplot as plt

def make_graph(data_file, hdf5_file):
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

    #open reaction rate file
    f=h5py.File(hdf5_file, 'r')
    group_1=f['422']

    """
    #net rates of progress for each reaction
    rxn_prod_rates=np.array(group_1['Reaction Production Rates'])
    #net production rates for each species
    sp_prod_rates=np.array(group_1['Species Net Production Rates'])
    #orignal net production rates for each species
    sp_prod_rates_original=np.array(group_1['Species Net Production Rates Original'])
    #print sp_prod_rates_original[0]
    #print sp_prod_rates[6]
    """
    #initalize solution and components
    solution= ct.Solution(data_file)
    species_objects=solution.species()
    reaction_objects=solution.reactions()


    G=nx.MultiDiGraph()
    node_labels={}
    for nm, grp in f.iteritems():
        #net rates of progress for each reaction
        rxn_prod_rates=np.array(grp['Reaction Production Rates'])
        #net production rates for each species
        sp_prod_rates=np.array(grp['Species Net Production Rates'])
        #orignal net production rates for each species
        sp_prod_rates_original=np.array(grp['Species Net Production Rates Original'])

        for l, m in enumerate(species_objects):
            sp_A=m.name
            if G.has_node(sp_A) is False:
                G.add_node(sp_A)
                node_labels[sp_A] =l
            ri_total={}
            ri_partial={}
            #iterate through every species
            for k, n in enumerate(species_objects):
                name=n.name
                #if name==sp_A:
                    #continue

                ri_total[name] = 0
                ri_partial[name] = 0
                #for every species iterate through every reaction
                for i, reac in enumerate(reaction_objects):
                    reac_type=reac.__class__.__name__
                    #if reac.duplicate is False:
                    if reac.__contains__(name):
                        rxn_prod_rate=float(rxn_prod_rates[i])
                        product_species=reac.products
                        reactant_species=reac.reactants
                        if name in product_species and reactant_species:
                            if rxn_prod_rate > 0:
                                total_species = reactant_species
                            if rxn_prod_rate ==0:
                                total_species=reac.products
                                total_species.update(reac.reactants)
                            else:
                                total_species = product_species
                        else:
                            total_species=reac.products
                            total_species.update(reac.reactants)
                        #if name in total_species.keys():
                        molar_coeff= float(total_species[name])
                        ri_total[name] += (rxn_prod_rate*molar_coeff)
                        if sp_A in total_species.keys():
                            ri_partial[name] +=  (rxn_prod_rate*molar_coeff)
                        try:
                            weight_total= float(ri_partial[name])/float(ri_total[name])
                        except ZeroDivisionError:
                            weight_total = 0.0

                        edge_name=sp_A + '_'+ name
                        if G.has_edge(name, sp_A) is False:
                            if weight_total != 0.0 or -0.0:
                                G.add_edge(name, sp_A, weight=weight_total, color='b')
                            else:
                                pass
                                #G.add_edge(name, sp_A, weight=weight_total, color='r')
        print (nm, G.number_of_edges())


    nx.draw(G, with_labels =True, width=.25)

    plt.show()
    print G.number_of_edges()
    print G.number_of_nodes()



make_graph('gri301.cti', 'production_rates.hdf5')
