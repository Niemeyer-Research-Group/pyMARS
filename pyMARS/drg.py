import networkx as nx
import cantera as ct
import numpy as np
import h5py
import matplotlib.pyplot as plt
import time

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

    #initalize solution and components
    solution= ct.Solution(data_file)
    species_objects=solution.species()
    reaction_objects=solution.reactions()
    #ask for threshold value
    threshold=float(raw_input('Enter threshold value: '))


    G=nx.MultiGraph()
    node_labels={}
    count=0
    start_time=time.time()

    #add nodes to graph
    for ind, sp in enumerate(species_objects):
        G.add_node(sp.name)

    #iterate through each timestep
    for nm, grp in f.iteritems():
        count +=1
        rxn_prod_rates=np.array(grp['Reaction Production Rates'])
        #generate dict of sum production Rates
        ri_total={}
        ri_partial={}

        for k in species_objects:
            ri_total[k.name] = 0

        #iterate through reactions and build weights
        for i, reac in enumerate(reaction_objects):
            products=reac.products
            product_names=reac.products.keys()

            reactants=reac.reactants
            reactant_names=reac.reactants.keys()

            reaction_production_rate=float(rxn_prod_rates[i])

            if reaction_production_rate != 0:
                for reactant in reactant_names:
                    molar_coeff=float(reactants[reactant])
                    ri_total[reactant] += (reaction_production_rate*molar_coeff)

                    for product in product_names:
                        partial_name= reactant + '_' + product

                        if product == reactant:
                            continue
                        try:
                            ri_partial[partial_name] += (reaction_production_rate*molar_coeff)
                        except KeyError:
                            ri_partial[partial_name] = (reaction_production_rate*molar_coeff)

        #divide progress related to species B by total progress
        #and make edge
        for ind in ri_partial:
            try:
                both=ind.split('_', 1)
                sp_A=both[0]
                sp_B=both[1]
                weight = float(ri_partial[ind])/float(ri_total[sp_A])
                if weight >= threshold:
                    G.add_edge(sp_A, sp_B, weight=weight)
            except IndexError:
                print ind
                continue

        progress= str(count) + '/40 timesteps'
        print progress
        #print("--- %s seconds ---" % (time.time() - start_time))
        #print (nm, G.number_of_edges())

    #get connected species
    target=str(raw_input('Enter target starting species:'))
    essential_nodes=list(nx.node_connected_component(G, target))

    nx.draw(G, with_labels=True, width=.25)

    #get list of species to eliminate
    exclusion_list=list()
    ex_list=[]
    for spec in solution.species():
        ind_name=spec.name
        if ind_name not in essential_nodes:
            exclusion_list.append(spec.name)
            ex_list.append(spec.name)
    exclusion_list_string='\''
    for spc in exclusion_list:
        exclusion_list_string += spc + ', '
    exclusion_list_string= exclusion_list_string.rstrip(',')
    print 'Species to Exclude'
    print exclusion_list_string
    plt.show(block=False)
    print G.number_of_edges()
    print G.number_of_nodes()
    return ex_list


#make_graph('gri301.cti', 'production_rates.hdf5')



#original DRG code. takes too long
"""
for nm, grp in f.iteritems():
    count +=1
    #net rates of progress for each reaction
    rxn_prod_rates=np.array(grp['Reaction Production Rates'])
    #net production rates for each species
    #sp_prod_rates=np.array(grp['Species Net Production Rates'])
    #orignal net production rates for each species
    #sp_prod_rates_original=np.array(grp['Species Net Production Rates Original'])

    for l, m in enumerate(species_objects):
        sp_A=m.name
        #if G.has_node(sp_A) is False:
        if count == 1:
            G.add_node(sp_A)
            node_labels[sp_A] =l
        ri_total={}
        ri_partial={}
        #iterate through every species
        for k, n in enumerate(species_objects):
            name=n.name
            if name==sp_A:
                continue
            ri_total[name] = 0
            ri_partial[name] = 0
            #for every species iterate through every reaction
            for i, reac in enumerate(reaction_objects):
                #reac_type=reac.__class__.__name__
                #if reac.duplicate is False:
                if reac.__contains__(name):
                    rxn_prod_rate=float(rxn_prod_rates[i])
                    product_species=reac.products
                    reactant_species=reac.reactants
                    if name in product_species and reactant_species:
                        if rxn_prod_rate != 0:
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
                            if weight_total > threshold:
                                G.add_edge(name, sp_A, weight=weight_total, color='b')
    """
