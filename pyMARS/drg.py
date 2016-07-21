import networkx as nx
import cantera as ct
import numpy as np
import h5py

def make_graph(data_file, sim_data):
    """ Use the Direct Relation Graph (DRG) method to choose species to
    eliminate from reaction mechanism.

    Parameters:
    -------------
    data_file: Mechanism file in the form of .cti

    Returns:
    -------------
    exclusion_list: List of species to trim from mechanism
    """

    fname= sim_data
    f= h5py.File(fname, "r")
    class data(object):
            times=np.array(f.get('Times'))
            temps=np.array(f.get('Temps'))
            species=f.get('Species_Data')



    """
    for n in data.species:
        if n == 'OH':
            print n
            sp =np.array(data.species.get(n))
            print sp.shape

    sp=sp.astype(float)
    o=sp[1]
    z=np.diff(sp)
    """

    #initalize solution and components
    solution= ct.Solution(data_file)
    species_objects=solution.species()
    reaction_objects=solution.reactions()

    #only want elementary reactions
    for n, rxn in enumerate(reaction_objects):
        if rxn.reaction_type is not 1: #1 is for ElementaryReaction
            reaction_objects.remove(rxn)


    #Make graph and nodes (each species is a node)
    G=nx.MultiDiGraph()
    for i, species in enumerate(species_objects):
        G.add_node(i, name=species.name)
        mass_values = np.array(data.species.get(species.name))
        mass_values = mass_values.astype(float) #convert to float
        G.node[i]['masses'] = mass_values       #assign mass values from simulation to each node
        G.node[i]['rates'] = np.diff(mass_values)   #get rate values for each timestep

    print len(G.node[0] ['rates'])

    dt=np.diff(data.times)
    print type(reaction_objects[0]).__name__
    for i, t in enumerate(dt):          #enumerate over every time step
        if i == 0:                      #temporary measure to only construct DRG for first timestep
            for n, rx in enumerate(reaction_objects):
                species_contained= rx.products.keys() + rx.reactants.keys()




make_graph('gri301.cti', 'pym_trimmed_gri301_species_data.hdf5')
