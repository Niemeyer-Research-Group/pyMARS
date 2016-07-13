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

    times=np.array(f.get('Times'))
    temps=np.array(f.get('Temps'))
    species_data=f.get('Species_Data')
    

    solution= ct.Solution(data_file)
    species_objects=solution.species()
    reaction_objects=solution.reactions()

    #only want elementary reactions

    for j, rxn in enumerate(reaction_objects):
        equation_type = type(reaction_objects[j]).__name__
        if equation_type != 'ElementaryReaction':
            reaction_objects.remove(reaction_objects[j])

    print (str(len(reaction_objects)) + '  Elementary Reactions\n')







    #Make graph and nodes
    graph=nx.MultiDiGraph()
    for i, species in enumerate(species_objects):
        graph.add_node(species_objects[i])




make_graph('gri301.cti', 'pym_trimmed_gri301_species_data.hdf5')
