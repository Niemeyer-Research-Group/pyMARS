import numpy as np

from .simulation import Simulation

def setup_simulations(conditions_array, model):
    """Sets up simulations based on list of initial conditions

    Parameters
    ----------
    conditions_array : list of Condition
        List of initial conditions read from input file
    model : cantera.Solution
        Model to use for simulations
    
    Returns
    -------
    sim_array : list of Simulation
        List of simulations to be performed

    """
    sim_array = []
    i = 0
    while (i < len(conditions_array)): #For all conditions
        properties = {} #create properties dictionary
        properties['temperature'] = float(conditions_array[i].temperature)
        properties['pressure'] = float(conditions_array[i].pressure)
        properties['equivalence_ratio'] = float(conditions_array[i].equi)
        properties['fuel'] = conditions_array[i].fuel
        properties['oxidizer'] = conditions_array[i].oxid

        # create simulation object and add it to the list
        sim_array.append(Simulation(i, properties, model))
        i = i + 1

    return sim_array


def simulate(sim_array):
    """Performs simulations and processes results based on input list.

    Parameters
    ----------
    sim_array : list of Condition

    Returns
    -------
    ignition_delay : numpy.array
        Array of simulated ignition delays

    """
    tau = [] #Ignition delays
    sample_points = [] #Sample information
    for case in sim_array: #Run simulations for original model and process results
        tau.append(case.run_case())
        sample_points.append(case.process_results())

    ignition_delay = np.array(tau) #Turn tau array into a numpy array
    return ignition_delay
