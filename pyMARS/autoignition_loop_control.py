import cantera as ct
from autoignition_module import run_sim

def autoignition_loop_control(solution_object, args='none'):
    """Controls autoignition module

    :param solution_object:
        Cantera solution object
    :param args:
        Arguments from terminal
    """
    print dir(args)
    if args is 'none':
        if 'initial_sim' in args:
            initial_sim = True
        else:
            initial_sim = False
    else:

        if args.initial_sim is True:
            initial_sim = True
        else:
            initial_sim = False
    if initial_sim is True:
        args.frac = raw_input('Enter mole fractions (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich) :  ')
        args.Temp = float(raw_input('Enter Solution Temperature in (K):'))
    else:
        frac = args.frac
        initial_temperature = args.Temp

    sim_result = run_sim(solution_object, args)

    return sim_result
