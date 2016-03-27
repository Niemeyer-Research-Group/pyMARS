#test
"""
file:
    create_trimmed_model

test parameter:
    number and name of species eliminated

"""
#import cantera
import cantera as ct

#import function
from .. import create_trimmed_model as func

#define function inputs
data_file=('gri30.xml')
exclusion_list = ['02', 'CO2']



#call Function
result=func.create_trimmed_model(data_file, exclusion_list)

print(result)

initial_solution = result[0]
new_solution = result[1]

n_initial_species= len(initial_solution.species())
n_new_species= len(new_solution.species())

def test_n_species():
    assert n_initial_species > n_new_species

n_initial_reactions = len(initial_solution.reactions())
n_new_reactions = len(new_solution.reactions())

def test_n_reactions():
    assert n_initial_reactions > n_new_reactions


#output for user
#print('start %s initial species, %s initial reactions') %(initial_solution.n_species, initial_solution.n_reactions)
#print('end %s final species, %s final reactions')   %(len(final_species_objects), len(final_reaction_objects))
