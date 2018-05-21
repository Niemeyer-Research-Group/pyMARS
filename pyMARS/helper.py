##############
# This function takes the conditions array and converts it into an array of simulation species_objects
#
# conditions array: array of initial conditions read in from the ic file
#  model: model to base the simulation off of
#
# simulation array: an array of unran simulation objects from simulation.py
##############

from simulation import Simulation

def setup_simulations(conditions_array, model):
    sim_array = []
    i = 0
    while (i < len(conditions_array)): #For all conditions
        properties = {} #create properties dictionary
        properties['temperature'] = float(conditions_array[i].temperature)
        properties['pressure'] = float(conditions_array[i].pressure)
        properties['equivalence_ratio'] = float(conditions_array[i].equi)
        properties['fuel'] = conditions_array[i].fuel
        properties['oxidizer'] = conditions_array[i].oxid

        sim_array.append(Simulation(i,properties,model,'./h5files')) #create simulation object and add it to the array
        i = i + 1

    return sim_array
