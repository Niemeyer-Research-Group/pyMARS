import cantera as ct
import h5py
from autoignition_module import run_sim
from readin_initial_conditions import readin_conditions


def autoignition_loop_control(solution_object, args):
    """Controls autoignition module

    :param solution_object:
        Cantera solution object
    :param args:
        Arguments from terminal
    """
    conditions_array = readin_conditions(str(args.conditions_file))
    print 'read in conditions'
    if args.initial_sim is True:
        tau_array = []
        initial_temperature_array = []
        for condition in conditions_array:
            sim_result = run_sim(solution_object, condition, args)

            initial_temperature_array.append(condition.temperature)
            try:
                tau_array.append(sim_result.tau)
            except AttributeError:
                tau_array.append(0.0)
        sim_result.tau_array = tau_array
        sim_result.initial_temperature_array = initial_temperature_array
    else:
        frac = args.frac
        initial_temperature = args.Temp
        sim_result = run_sim(solution_object, args)
    return sim_result
