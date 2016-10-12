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
        t_low = float(args.Temp.split('-')[0])
        t_high = float(args.Temp.split('-')[1])
        t_array = [t_low]
        while t_low <= (t_high-100.0):
            t_low += 100.0
            t_array.append(t_low)
        print t_array
        for condition in t_array:
            args.Temp = float(condition)
            print args.Temp
            sim_result = run_sim(solution_object, args)
            #sim_result.test.keys() results in [u'1700.0', u'1800.0']
    else:
        frac = args.frac
        initial_temperature = args.Temp
        args.Temp = float(raw_input(' else case Enter Solution Temperature to test error in (K) (ex. 800-1000):'))
        sim_result = run_sim(solution_object, args)
    return sim_result
