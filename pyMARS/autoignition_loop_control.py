import cantera as ct
import h5py
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
    sim_result
        Object containing autoignition results
    """
    try:
        conditions_array = readin_conditions(str(args.conditions_file))
    except Exception:
        class condition():
            name = str(raw_input('Enter condition #: '))
            pressure = float(raw_input('Initial Pressure (kPa): '))
            temperature = float(raw_input('Initial Temperature (K): '))
            reactants = str(raw_input('List of reactants and mole fractions (ex.Ch4:1,o2:11): '))
            reactant_list = reactants.split(',')
            species={}
            for reactant in reactant_list:
                species[reactant.split(':')[0]] = reactant.split(':')[1]
        conditions_array = condition()
        args.multiple_conditions = False

    tau_array = []
    initial_temperature_array = []

    if args.multiple_conditions is True:
        for condition in conditions_array:
            sim_result = run_sim(solution_object, condition, args)
            initial_temperature_array.append(condition.temperature)
            try:
                tau_array.append(sim_result.tau)
            except AttributeError:
                tau_array.append(0.0)
    else:
        sim_result = run_sim(solution_object, condition, args)
        initial_temperature_array.append(condition.temperature)
        try:
            tau_array.append(sim_result.tau)
        except AttributeError:
            tau_array.append(0.0)

    sim_result.tau_array = tau_array
    sim_result.initial_temperature_array = initial_temperature_array
    return sim_result
