from autoignition_module import run_sim
from readin_initial_conditions import readin_conditions


def autoignition_loop_control(solution_object, args, plot=False):
    """Controls autoignition module

    Parameters
    ----------
    solution_object : class
        Cantera solution object
    args : class
        Arguments from terminal
    plot : boolean
        A boolean value that represents if this run will use additional features such as making plots or csv files.  

    Returns
    -------
    sim_result : obj
        Object containing autoignition results
            .time
            .temp
            .initial_temperature_array
            .sp_data
            .test (h5py object)
            .tau
            .Temp
            .frac
    """
    class condition(object):
        def __init__(self, name):
            self.name = str(name)
            self.pressure = float(raw_input('Initial Pressure (atm): '))
            self.temperature = float(raw_input('Initial Temperature (K): '))
            self.reactants = str(raw_input('List of reactants and mole fractions (ex.Ch4:1,o2:11): '))
            self.reactant_list = self.reactants.split(',')
            self.species={}
            for reactant in self.reactant_list:
                self.species[reactant.split(':')[0]] = reactant.split(':')[1]
    conditions_array = []

    #get initial conditions from readin, or prompt user if none given
    if args.conditions_file:
        conditions_array = readin_conditions(str(args.conditions_file))
    elif not args.conditions_file:
        print('No initial conditions file found. Please enter below')
        number_conditions = int(float(raw_input('Enter # of initial conditions: ')))
        if number_conditions > 1.0:
            for i in range(number_conditions):
                conditions_array.append(condition(i))
        else:
            conditions_array.append(condition(1.0))
            args.multiple_conditions = False

    tau_array = []
    initial_temperature_array = []

    #send initial conditions to autoignition script
    i = 1 #iterate which condition number it is for file naming.
    for condition in conditions_array:
        sim_result = run_sim(i,solution_object, condition, args, plot)
        if (sim_result == 0):
            return 0
        initial_temperature_array.append(condition.temperature)
        i = i + 1
        try:
            tau_array.append(sim_result.tau)
        except AttributeError:
            tau_array.append(0.0)
    
    #Store the information from the simulation in an object and return it.
    sim_result.tau_array = tau_array
    sim_result.initial_temperature_array = initial_temperature_array
    return sim_result
