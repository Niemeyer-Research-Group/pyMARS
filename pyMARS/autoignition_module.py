import sys
import numpy as np
import cantera as ct
import os
from get_sample_range import get_range
import h5py


def run_sim(solution_object, sys_args='none', **usr_args ):
    """
    Function to run Cantera reactor simulation

    Arguments
        Cantera mechanism file
        User arguments (plot='y',
                        points='y',
                        writehdf5='y',
                        writecsv='y',
                        initial_sim='y')
        *User argument initial_sim needed if running module independently
    ----------
    Output
        Plot of Temp vs Time
        Points of interest
        CSV file
        Hdf5 file
    ----------
    Example
        run_sim(gas_solution, points='y', plot='y', initial_sim='y')

    """
    solution = solution_object
    # temporary fix, allow initial conditions to be carried in if already set from previous
    #sim.
    if sys_args is 'none':
        if 'initial_sim' in usr_args:
            initial_sim = True
    else:
        if sys_args.initial_sim is True:
            initial_sim = True
        else:
            initial_sim = False
    if initial_sim is True:
        frac = raw_input('Enter mole fractions (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich) :  ')
        initial_temperature = float(raw_input('Enter Solution Temperature in (K):'))
        initial_temperature = initial_temperature
    else:
        frac = sys_args.frac
        initial_temperature = sys_args.Temp
    solution.TPX = initial_temperature, ct.one_atm, frac #1001.325 kPa
    species = solution.species()
    reactions = solution.reactions()

    #run sim to find ignition delay from dT/dt max
    reactor = ct.Reactor(solution)
    simulation = ct.ReactorNet([reactor])
    current_time = 0.0
    stop_time = 5.0e-3
    group_index = 0
    times1 = []
    temps = [] #first column is time, second is temperature
    mass = reactor.mass
    sdata = np.zeros([0, len(reactor.Y)])
    production_data = np.zeros([0, len(solution.net_production_rates)])
    state_list = list()

    f1 = h5py.File('mass_fractions.hdf5', 'w')
    while current_time < stop_time:
        group_index += 1
        current_time = simulation.step(stop_time)
        times1.append(current_time)
        temps.append(reactor.T)
        species_data = reactor.Y

        grp = f1.create_group(str(group_index))
        grp['Temp'] = reactor.T
        grp['Time'] = current_time
        grp['Pressure'] = reactor.thermo.P
        species_production_rates = reactor.thermo.net_production_rates
        grp.create_dataset('Species Mass Fractions', data=species_data)
        grp.create_dataset('Species Net Production Rates Original', data=species_production_rates)

        species_data = species_data[:,np.newaxis].T #translate from [n, 1] to [1,n]
        sdata = np.vstack((sdata, species_data))

        production_rates = np.array(solution.net_production_rates)
        production_rates = production_rates[:,np.newaxis].T
        production_data = np.vstack((production_data, production_rates))
    print('\n')

    function_data = get_range(times1,temps,sdata, production_data)

    tau = function_data.tau
    i = function_data.index
    times1 = function_data.times
    temps = function_data.temps
    sdata = function_data.species_data
    production_data = function_data.production_data
    times_total = function_data.times_total
    temps_total = function_data.temps_total
    sdata_total = function_data.species_data_total
    production_data_total = function_data.production_data_total
    deriv_max =function_data.derivative_max
    initial_point = function_data.initial_point
    final_point = function_data.final_point

    for grp in f1.keys():
        if int(grp) not in range((i-20), (i+20)):
            f1.__delitem__(str(grp))

    f1.close()

    #utility functions
    def plot():
        import matplotlib.pyplot as plt
        plt.clf()
        #plot combustion point
        plt.plot(deriv_max[0], deriv_max[1], 'ro', ms=7, label= 'ignition point')
        #plot initial and final sample points
        plt.plot(initial_point[0], initial_point[1], 'rx', ms=5, mew=2)
        plt.plot(final_point[0], final_point[1], 'rx', ms=5, mew=2)
        #plot temp vs time
        plt.plot(times_total, temps_total)
        plt.xlabel('Time (s)')
        plt.title('Mixture Temperature vs Time')
        plt.legend()
        plt.ylabel('Temperature (C)')
        #plt.axis([0, 1.2, 900, 2800])
        plt.show()

    def writecsv(sdata):
        #format matrix for csv
        names = str(solution.species_names)
        tt = ['Time (s)', 'Temp (K)']
        names = solution.species_names
        name_array = np.append(tt, names)
        sdata = sdata.astype('|S10')
        file_data = np.vstack((name_array, sdata))
        #open and write to file
        input_file_name_stripped = os.path.splitext(data_file)[0]
        output_file_name = os.path.join(input_file_name_stripped + '_species_data.csv')
        print output_file_name
        with open(output_file_name, 'wb') as f:
            np.savetxt(f, file_data, fmt=('%+12s'),  delimiter=',')
        #os.system('atom '+ output_file_name)

    def writehdf5(sdata):
        #format matrix for hdf5
        names = str(solution.species_names)
        tt = ['Time (s)', 'Temp (K)']
        names = solution.species_names
        name_array = np.append(tt, names)
        sdata = sdata.astype('|S10')
        file_data = np.vstack((name_array, sdata))
        #open and write to file
        input_file_name_stripped = os.path.splitext(data_file)[0]
        output_file_name = os.path.join(input_file_name_stripped + '_species_data.hdf5')
        with h5py.File(output_file_name, 'w') as f:
            Times = f.create_dataset("Times", data=times1)
            Temps = f.create_dataset("Temps", data=temps)
            sgroup= f.create_group('Species_Data')
            for i, sp in enumerate(solution.species_names):
                    sgroup.create_dataset(sp, data=sdata[:,i+2])

    def points():
        print("\nTime[s]            Temp[K]        Index        Point")
        print( str(initial_point[0]) +  "       " + str("{0:.2f}".format(initial_point[1]))\
         + "       " + str(initial_point[2]) + "     " + "Initial sample point")
        print(str(tau) + "        " + str("{0:.2f}".format(deriv_max[1])) + "       " + str(deriv_max[2])\
                + "     " + "Ignition point")
        print( str(final_point[0]) +  "       " + str("{0:.2f}".format(final_point[1]))\
         + "       " + str(final_point[2]) + "     " + "Final sample point")

    #terminal use case
    if sys_args is not 'none':
        if sys_args.plot:
            plot()
        if sys_args.writecsv:
            writecsv(sdata)
        if sys_args.writehdf5:
            writehdf5(sdata)
        if sys_args.points:
            points()
    #individual use case
    if sys_args is 'none':
        if 'plot' in usr_args:
            plot()
        if 'writecsv' in usr_args:
            writecsv(sdata)
        if 'writehdf5' in usr_args:
            writehdf5(sdata)
        if 'points' in usr_args:
            points()


    class return_obj:
        def __init__(self, time, temp, sp_data, f1, tau, Temp, frac):
            self.time = time
            self.temp = temp
            self.sp_data = sp_data
            self.test = f1
            self.tau = tau
            self.Temp = initial_temperature
            self.frac = frac

    return return_obj(times1, temps, sdata, f1, tau, initial_temperature, frac)

    "sdata is an array of 40 timesteps, with each instance containing an array of species"
    "mass fractions at that instant"
