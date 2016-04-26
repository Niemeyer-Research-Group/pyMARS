# Local Imports
import cantera as ct


def create_trimmed_model(data_file, exclusion_list):

    """ Function to reduce list of species and corresponding reactions.

    Parameters
    ----------
    data_file:
        local .cti or .xml data file containing mechanism information
    exclusion_list:
        List of species that will be trimmed

    Returns
    -------
    new_solution :
        Original Cantera Solution Object
        Simplified Cantera Solution Object
    """

    # define initial solution objects
    initial_solution = ct.Solution(data_file)

    initial_species_objects = initial_solution.species
    initial_species_names = initial_solution.species_names

    initial_reaction_list = initial_solution.reaction_equations()
    initial_reaction_objects = initial_solution.reactions

    # Remove reactions
    initial_index=[]
    index_exclude=[]
    for i, reaction_string in enumerate(initial_reaction_list):
        initial_index.append(i)
        Reaction = initial_reaction_list[i]
        for n in exclusion_list:
            if n in Reaction:
                index_exclude.append(i)
    final_index=initial_index

    for value in index_exclude:
        if value in final_index:
            final_index.remove(value)

    final_reaction_objects= [initial_solution.reaction(r) for r in final_index]
    final_T=initial_solution.T
    final_P=initial_solution.P

    # Remove Species
    final_species_names=initial_species_names
    for n in exclusion_list:
        if n in initial_species_names:
            final_species_names.remove(n)

    final_species_objects =   [initial_solution.species(name)
                                for name in final_species_names]

    # New solution definition
    new_solution= ct.Solution(  species=final_species_objects,
                        reactions=final_reaction_objects,
                        thermo='IdealGas', kinetics='GasKinetics')
    new_solution.TP = initial_solution.TP

    return (initial_solution, new_solution)
