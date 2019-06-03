"""Module for general model reduction class and function."""
import os
from typing import NamedTuple

import cantera as ct

class ReducedModel(NamedTuple):
    """Represents reduced model and associated metadata
    """
    model: ct.Solution
    filename: str = ''
    error: float = 0.0
    limbo_species: list = []
    

def trim(initial_model_file, exclusion_list, new_model_file):
    """Function to eliminate species and corresponding reactions from model

    Parameters
    ----------
    initial_model_file : str
        Filename for initial model to be reduced
    exclusion_list : list of str
        List of species names that will be removed
    new_model_file : str
        Name of new reduced model file

    Returns
    -------
    new_solution : ct.Solution
        Model with species and associated reactions eliminated

    """
    solution = ct.Solution(initial_model_file)

    # Remove species if in list to be removed
    final_species = [sp for sp in solution.species() if sp.name not in exclusion_list]
    final_species_names = [sp.name for sp in final_species]

    # Remove reactions that use eliminated species
    final_reactions = []
    for reaction in solution.reactions():
        reaction_species = list(reaction.products.keys()) + list(reaction.reactants.keys())
        if all([sp in final_species_names for sp in reaction_species]):
            # remove any eliminated species from third-body efficiencies
            if hasattr(reaction, 'efficiencies'):
                reaction.efficiencies = {
                    sp:val for sp, val in reaction.efficiencies.items() 
                    if sp in final_species_names
                    }
            final_reactions.append(reaction)

    # Create new solution based on remaining species and reactions
    new_solution = ct.Solution(
        species=final_species, reactions=final_reactions,
        thermo='IdealGas', kinetics='GasKinetics'
        )
    new_solution.TP = solution.TP
    new_solution.name = os.path.splitext(new_model_file)[0]

    return new_solution
