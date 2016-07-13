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
    print len(reaction_objects)

    #only want elementary reactions
    for n, rxn in enumerate(reaction_objects):
        if rxn.reaction_type is not 1: #1 is for ElementaryReaction
            reaction_objects.remove(rxn)


    #Make graph and nodes
    G=nx.MultiDiGraph()
    for i, species in enumerate(species_objects):
        G.add_node(i, name=species.name)
    #   G.node[i]['masses'] =

    dt=np.diff(data.times)
    print type(reaction_objects[0]).__name__
    for i, t in enumerate(dt):          #enumerate over every time step
        for n, rx in enumerate(reaction_objects):
            pass



make_graph('gri301.cti', 'pym_trimmed_gri301_species_data.hdf5')
