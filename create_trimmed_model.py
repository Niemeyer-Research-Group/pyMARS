

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


    #define initial reaction lists
    ReactionList    = 	initial_solution.reaction_equations()
    ReactionObjects =	initial_solution.reaction
    list	=	[]
    #for every reaction in the list
    for i, Rxn in enumerate(ReactionList):

        Reactants 	= 	ReactionObjects(i).reactants
        Products	= 	ReactionObjects(i).products

        #for every species in the exclusion list
        for k in exclusion_list:
            #if the species is not in the reactants or products list, add it to the blank list
            if k not in Reactants and k not in Products:
                list.append(ReactionObjects(i))



    #need to find a way to do this after eliminating reactions
    """for n in exclusion_list:
        if n in initial_species:
            initial_species.remove(n)
    """
    #define initial species list
    Species_Objects =   [initial_solution.species(name) for name in initial_species]

    ReactionObjects = list
    new_solution= ct.Solution(species=Species_Objects, reactions=ReactionObjects,
                    thermo='IdealGas', kinetics='GasKinetics')

    print('start %s initial species') %initial_solution.n_species
    print('end %s final species') % new_solution.n_species
    print('start %s initial reactions') %initial_solution.n_reactions
    print('end %s final reactions') %len(ReactionObjects)

    x=initial_solution.reaction(0)
    print(x)
    a=x.reactant_string

    print(a)

#calling the function
#list to exclude
SPexc=['O2', 'CO2'];
create_trimmed_model("gri30.cti", SPexc)
