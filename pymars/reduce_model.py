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
    

def trim(initial_model_file, exclusion_list, new_model_file, phase_name=''):
    """Function to eliminate species and corresponding reactions from model

    Parameters
    ----------
    initial_model_file : str
        Filename for initial model to be reduced
    exclusion_list : list of str
        List of species names that will be removed
    new_model_file : str
        Name of new reduced model file
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 

    Returns
    -------
    new_solution : ct.Solution
        Model with species and associated reactions eliminated

    """
    solution = ct.Solution(initial_model_file, phase_name)

    # Remove species if in list to be removed
    final_species = [sp for sp in solution.species() if sp.name not in exclusion_list]
    final_species_names = [sp.name for sp in final_species]

    # Remove reactions that use eliminated species
    final_reactions = []
    for reaction in solution.reactions():
        # remove reactions with an explicit third body that has been removed
        if hasattr(reaction, 'efficiencies') and not getattr(reaction, 'default_efficiency', 1.0):
            if (len(reaction.efficiencies) == 1 and 
                list(reaction.efficiencies.keys())[0] in exclusion_list
                ):
                continue

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
    if phase_name:
        new_solution.name = phase_name
    else:
        new_solution.name = os.path.splitext(new_model_file)[0]

    return new_solution
