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
        args.Temp = str(raw_input('Enter Solution Temperature Range in (K) (ex. 800-1000):'))
    else:
        frac = args.frac
        initial_temperature = args.Temp

    t_low = float(args.Temp.split('-')[0])
    t_high = float(args.Temp.split('-')[1])
    t_array = [t_low]
    combined_rates = h5py.File('combined_rate_file.hdf5', 'w')
    while t_low <= (t_high-100.0):
        t_low += 100.0
        t_array.append(t_low)
    print t_array

    for condition in t_array:
        args.Temp = float(condition)
        sim_result = run_sim(solution_object, args)

    return sim_result
