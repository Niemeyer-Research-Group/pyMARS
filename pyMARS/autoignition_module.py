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
                        writecsv='y',)
        *User argument initial_sim needed if running module independently
    ----------
    Output
        Plot of Temp vs Time
        Points of interest
        CSV file
        Hdf5 file
        mass_fractions.hdf5 : [initial_temp]
                                    [index]
                                        [Temp]
                                        [Time]
                                        [Pressure]
                                        [Species Mass Fractions]
                                        [Species Net Production Rates Original]

    ---------
    returns:
        sim result
            .time
            .temp
            .initial_temperature_array
            .sp_data
            .test (h5py object)
            .tau
            .Temp
            .frac
    ----------
    Example
        run_sim(gas_solution, points='y', plot='y', initial_sim='y')

    """
    solution = solution_object
    # temporary fix, allow initial conditions to be carried in if already set from previous
    #sim.
    print 'running sim from initial_temp %s' %sys_args.Temp
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

    f1 = h5py.File('mass_fractions.hdf5', 'a')
    individual = f1.create_group(str(initial_temperature))
    while current_time < stop_time:
        group_index += 1
        try:
            current_time = simulation.step(stop_time)
        except Exception:
            error_string = 'Cantera autoignition_error @ %sK initial temperature' %initial_temperature
            print error_string
            return
        times1.append(current_time)
        temps.append(reactor.T)
        species_data = reactor.Y
        grp = individual.create_group(str(group_index))
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

    sample = get_range(times1,temps,sdata, production_data)

    for grp in f1[str(initial_temperature)].keys():
        if int(grp) not in range((sample.index-20), (sample.index+20)):
            f1[str(initial_temperature)].__delitem__(str(grp))

    #f1.close()

    #utility functions
    def plot():
        import matplotlib.pyplot as plt
        plt.clf()
        #plot combustion point
        plt.plot(sample.derivative_max[0], sample.derivative_max[1], 'ro', ms=7, label= 'ignition point')
        #plot initial and final sample points
        plt.plot(sample.initial_point[0], sample.initial_point[1], 'rx', ms=5, mew=2)
        plt.plot(sample.final_point[0], sample.final_point[1], 'rx', ms=5, mew=2)
        #plot temp vs time
        plt.plot(sample.times_total, sample.temps_total)
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
            Times = f.create_dataset("Times", data=sample.times)
            Temps = f.create_dataset("Temps", data=sample.temps)
            sgroup= f.create_group('Species_Data')
            for i, sp in enumerate(solution.species_names):
                    sgroup.create_dataset(sp, data=sdata[:,i+2])

    def points():
        print("\nTime[s]            Temp[K]        Index        Point")
        print( str(sample.initial_point[0]) +  "       " + str("{0:.2f}".format(sample.initial_point[1]))\
         + "       " + str(sample.initial_point[2]) + "     " + "Initial sample point")
        print(str(sample.tau) + "        " + str("{0:.2f}".format(sample.derivative_max[1])) + "       " + str(sample.derivative_max[2])\
                + "     " + "Ignition point")
        print( str(sample.final_point[0]) +  "       " + str("{0:.2f}".format(sample.final_point[1]))\
         + "       " + str(sample.final_point[2]) + "     " + "Final sample point")

    #terminal use case
    if sys_args is not 'none':
        if sys_args.plot:
            plot()
        if sys_args.writecsv:
            writecsv(sample.species_data)
        if sys_args.writehdf5:
            writehdf5(sample.species_data)
        if sys_args.points:
            points()
    #individual use case
    if sys_args is 'none':
        if 'plot' in usr_args:
            plot()
        if 'writecsv' in usr_args:
            writecsv(sample.species_data)
        if 'writehdf5' in usr_args:
            writehdf5(sample.species_data)
        if 'points' in usr_args:
            points()


    class return_obj:
        def __init__(self, time, temp, sp_data, f1, tau, Temp, frac):
            self.time = time
            self.temp = temp
            self.initial_temperature_array = []
            self.sp_data = sp_data
            self.test = f1
            self.tau = tau
            self.Temp = initial_temperature
            self.frac = frac
            sim_result.tau_array = []

    return return_obj(sample.times, sample.temps, sample.species_data, f1, sample.tau, initial_temperature, frac)

    "sdata is an array of 40 timesteps, with each instance containing an array of species"
    "mass fractions at that instant"
