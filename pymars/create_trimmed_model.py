# Local Imports
import cantera as ct
import os

def trim(initial_solution, exclusion_list, file_name):
    """Reduces model using list of species to be removed.

    Parameters
    ----------
    initial_solution : obj
        Cantera solution object of initial model
    exclusion_list : list
        List of species that will be trimmed
    file_name : str 
        Name of original model file to be reduced

    Returns
    -------
        Trimmed Cantera solution object

    """

    # Remove reactions that use trimmed species
    final_reactions = []
    for reaction in initial_solution.reactions():
        reaction_species = list(reaction.products.keys()) + list(reaction.reactants.keys())
        if all([sp not in reaction_species for sp in exclusion_list]):
            final_reactions.append(reaction)
        

    # Remove species if in list to be removed
    final_species = [initial_solution.species(sp) 
                     for sp in initial_solution.species_names if sp not in exclusion_list
                     ]
                             
    # New solution definition
    new_solution= ct.Solution(species=final_species,
                              reactions=final_reactions,
                              thermo='IdealGas',
                              kinetics='GasKinetics'
                              )
    new_solution.TP = initial_solution.TP
    new_solution.name = 'trimmed_' + os.path.splitext(file_name)[0]

    return new_solution
