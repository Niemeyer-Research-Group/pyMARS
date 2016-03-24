

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
    #<---------------------------define initial parameters----------------->
    initial_solution = ct.Solution(data_file)

    initial_species_objects = initial_solution.species   #this is the index of species objects (referenced with parenthesis)
    initial_species_names  = initial_solution.species_names #this is a index of species name strings (referenced with brackets)

    initial_reaction_list    = 	initial_solution.reaction_equations() #index of reactions strings (referenced with brackets)
    ReactionObjects =	 initial_solution.reactions   #index of reaction objects (not callable)

    #<----------------------------Reaction Removal ------------------------->
    initial_indices=[]
    for v, m in enumerate(ReactionObjects()):
            initial_indices.append(v)

    indices_exclude=[] #these will be indices of reactions containing excluded species
    #for every reaction in the list
    for i, k in enumerate(initial_reaction_list):
        Reaction = initial_reaction_list[i]
        for z in exclusion_list:
            if z in Reaction:
                #print(Reaction)
                indices_exclude.append(i)       #adds to index

    for value in indices_exclude:
        if value in initial_indices:
            initial_indices.remove(value)
    final_indices=initial_indices


    final_reaction_objects= [initial_solution.reaction(r) for r in final_indices]


    #<---------------------------Species Removal ------------------------>

    #uses the name index to iterate over and remove.
    final_species_names=initial_species_names
    for n in exclusion_list:
        if n in initial_species_names:
            final_species_names.remove(n)
    #gets the species objects from the name index
    final_species_objects =   [initial_solution.species(name) for name in final_species_names]
    #final_species_objects into the new definition

    #<----------------------New Solution Definition ----------------------->

    #this new definition requires lists of both the species and reaciton objects
    #works
    #new_solution= ct.Solution(species=initial_solution.species() , reactions=initial_solution.reactions(), thermo='IdealGas', kinetics='GasKinetics')
    #in progress
    new_solution= ct.Solution(species=final_species_objects, reactions=final_reaction_objects, thermo='IdealGas', kinetics='GasKinetics')


    #output for user
    print('start %s initial species, %s initial reactions') %(initial_solution.n_species, initial_solution.n_reactions)
    print('end %s final species, %s final reactions') % (len(final_species_objects), len(final_reaction_objects))



#calling the function
#list to exclude
SPexc=['H2', 'O2'];
create_trimmed_model("gri30.cti", SPexc)
