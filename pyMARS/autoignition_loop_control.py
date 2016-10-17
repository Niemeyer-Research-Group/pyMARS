import cantera as ct
import h5py
from autoignition_module import run_sim


def autoignition_loop_control(solution_object, args):
    """Controls autoignition module

    :param solution_object:
        Cantera solution object
    :param args:
        Arguments from terminal
    """
    if args.initial_sim is True:
        args.frac = raw_input('Enter mole fractions (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich) :  ')
        args.Temp = str(raw_input('Enter Solution Temperature Range in (K) (ex. 800-1000):'))
        args.tau_array=[]
        if '-' in args.Temp:
            t_low = float(args.Temp.split('-')[0])
            t_high = float(args.Temp.split('-')[1])
            t_array = [t_low]
            while t_low <= (t_high-100.0):
                t_low += 200.0
                t_array.append(t_low)
            for condition in t_array:
                args.Temp = float(condition)
                sim_result = run_sim(solution_object, args)
                try:
                    args.tau_array.append(sim_result.tau)
                except AttributeError:
                    args.tau_array.append(0.0)
                #sim_result.test.keys() results in [u'1700.0', u'1800.0']
            sim_result.tau_array = args.tau_array
            sim_result.initial_temperature_array = t_array
        else:
            args.Temp = float(args.Temp)
            sim_result = run_sim(solution_object, args)
            print sim_result.tau_array

    else:
        frac = args.frac
        initial_temperature = args.Temp
        sim_result = run_sim(solution_object, args)
    return sim_result
