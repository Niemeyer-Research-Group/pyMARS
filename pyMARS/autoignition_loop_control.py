from autoignition_module import run_sim
from readin_initial_conditions import readin_conditions


def autoignition_loop_control(solution_object, args):
    """Controls autoignition module

    Parameters
    ----------
    solution_object : class
        Cantera solution object
    args : class
        Arguments from terminal

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
        print 'No initial conditions file found. Please enter below'
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
    if args.multiple_conditions is True:
        for condition in conditions_array:
            sim_result = run_sim(solution_object, condition, args)
	    if (sim_result == 0):
	        return 0
            initial_temperature_array.append(condition.temperature)
            try:
                tau_array.append(sim_result.tau)
            except AttributeError:
                tau_array.append(0.0)
    else:
        sim_result = run_sim(solution_object, conditions_array[0], args)
        initial_temperature_array.append(conditions_array[0].temperature)
        try:
            tau_array.append(sim_result.tau)
        except AttributeError:
            tau_array.append(0.0)

    sim_result.tau_array = tau_array
    sim_result.initial_temperature_array = initial_temperature_array
    return sim_result
