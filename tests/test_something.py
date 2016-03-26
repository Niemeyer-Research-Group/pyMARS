#test
"""
file:
    create_trimmed_model

test parameter:
    number and name of species eliminated

"""
#import cantera
import Cantera as ct

#import function
from create_trimmed_model import create_trimmed_model as test_func

#define function inputs
data_file=('gri30.xml')
exclusion_list = ['02', 'CO2']



#call Function
result=test_func(data_file, exclusion_list)

initial_solution = result[0]
new_solution = result[1]

print(initial_solution.nSpecies())

def test_answer():
    assert func(3) == 5

'''
#output for user
print('start %s initial species, %s initial reactions') %(initial_solution.n_species, initial_solution.n_reactions)
print('end %s final species, %s final reactions')   %(len(final_species_objects), len(final_reaction_objects))

'''
