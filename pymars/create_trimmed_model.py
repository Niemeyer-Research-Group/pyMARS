import os

def trim(initial_solution, exclusion_list, file_name):
    """Reduces model using list of species to be removed.

    Parameters
    ----------
    initial_solution : cantera.Solution
        Cantera solution object of initial model
    exclusion_list : list of str
        List of species names that will be removed
    file_name : str
        Name of original model file to be reduced

    Returns
    -------
    cantera.Solution
        Cantera solution object with reduced model

    """
    # Remove species if in list to be removed
    final_species = [sp for sp in initial_solution.species() if sp.name not in exclusion_list]

    # Remove reactions that use eliminated species
    final_reactions = []
    for reaction in initial_solution.reactions():
        reaction_species = list(reaction.products.keys()) + list(reaction.reactants.keys())
        if all([sp in final_species for sp in reaction_species]):
            # remove any eliminated species from third-body efficiencies
            if hasattr(reaction, 'efficiencies'):
                reaction.efficiencies = {
                    sp:val for sp, val in reaction.efficiencies.items() 
                    if sp in final_species
                    }

            final_reactions.append(reaction)

    # Create new solution based on remaining species and reactions
    new_solution= ct.Solution(species=final_species,
                              reactions=final_reactions,
                              thermo='IdealGas',
                              kinetics='GasKinetics'
                              )
    new_solution.TP = initial_solution.TP
    new_solution.name = 'reduced_' + os.path.splitext(file_name)[0]

    return new_solution
