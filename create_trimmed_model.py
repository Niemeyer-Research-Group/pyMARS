#function to read in an xml/cti file and remove species and reactions

import cantera as ct

def datareadin(data_file, exclusion_list):
    Solution = ct.Solution(data_file)
    Species = Solution.species_names

    for n in exclusion_list:
            if n in Species:
                Species.remove(n)
    SpeciesObjects	=	[Solution.species(name) for name in Species]

    ReactionList 	= 	Solution.reaction_equations()
    ReactionObjects	=	Solution.reaction
    list	=	[]
    for i, Rxn in enumerate(ReactionList):
        Reactants 	= 	ReactionObjects(i).reactants
        Products	= 	ReactionObjects(i).products
        for k in exclusion_list:
            if k not in Reactants and k not in Products:
                list.append(Solution.reaction(i))

    ReactionObjects = list
    New= ct.Solution(species=SpeciesObjects, reactions=ReactionObjects, thermo='IdealGas', kinetics='GasKinetics')
    print(Solution())
    print(New())


#calling the function
#list to exclude
SPexc=['O2'];
datareadin("gri30.cti", SPexc)
