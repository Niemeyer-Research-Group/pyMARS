

#Local Imports
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
    ReactionObjects :
        Short list of solution reaction objects
    New :
        Simplified Cantera Solution
    """

    initial_solution = ct.Solution(data_file)
    initial_species  = initial_solution.species_names

    for n in exclusion_list:
            if n in initial_species:
                initial_species.remove(n)

    Species_Objects =   [initial_solution.species(name) for name in initial_species]

    ReactionList    = 	initial_solution.reaction_equations()
    ReactionObjects =	initial_solution.reaction
    list	=	[]
    for i, Rxn in enumerate(ReactionList):
        Reactants 	= 	ReactionObjects(i).reactants
        Products	= 	ReactionObjects(i).products
        for k in exclusion_list:
            if k not in Reactants and k not in Products:
                list.append(initial_solution.reaction(i))

    ReactionObjects = list
    new_solution= ct.Solution(species=Species_Objects, reactions=ReactionObjects,
                    thermo='IdealGas', kinetics='GasKinetics')
    print(initial_solution.n_reactions)
    print(new_solution.n_reactions)


#calling the function
#list to exclude
SPexc=['O2'];
create_trimmed_model("gri30.cti", SPexc)
