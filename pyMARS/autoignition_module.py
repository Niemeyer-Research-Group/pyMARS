import sys
import numpy as np
import cantera as ct
import os
from get_sample_range import get_range
import h5py
import time as tm
import matplotlib
matplotlib.use("Agg")


def run_sim(i,solution_object, condition, sys_args='none', info=False,**usr_args):
    """
    Function to run Cantera reactor simulation for autoigntion conditions

    Parameters
    ----------
    i : int
        An integer value that represents which number initial condition the program is on
    solution_object : obj
        Cantera solution object
    condition
        An object contining initial conditions (temperature, pressure, mole fractions)
    info : boolean
        A boolean value that is true if additional features such as plotting or writing csv files will be used.  

    Returns
    ----------
    Output
        Plot of Temp vs Time
        Points of interest
        CSV file
        Hdf5 file
        mass_fractions.hdf5 : [initial_condition]
                                [index]
                                    [Pressure]
                                    [Reaction Rates of Progress]
                                    [Species Mass Fractions]
                                    [Species Net Production Rates Original]
                                    [Temp]
                                    [Time]

        return_obj : obj
            sim result
                .time
                .temp
                .initial_temperature_array
                .sp_data
                .test (h5py object)
                .tau
                .Temp
                .frac

    Example
    -------
        run_sim(gas_solution, points='y', plot='y', initial_sim='y')

    """
    #Set up variables
    func_start_time = tm.time()
    solution = solution_object
    initial_temperature = float(condition.temperature)
    pressure = float(condition.pressure)*float(ct.one_atm)
    
    #Set up species fractions
    frac = ''
    for reactant in condition.species.iteritems():
        if reactant[0] in solution.species_names:
            frac += str(reactant[0]) + ':' + str(reactant[1]) + ','
    frac = frac[:-1]
    
    #Set up variables
    solution.TPX = initial_temperature, pressure, frac #101325 Pa
    species = solution.species()
    reactions = solution.reactions()

    #run sim to find ignition delay from when the tempurature first reaches 400 Kelvin above its original value.
   
    #Set up variables for the simulation.
    reactor = ct.IdealGasReactor(solution) #Set up reactor
    simulation = ct.ReactorNet([reactor])
    current_time = 0.0
    stop_time = 25
    group_index = 0
    times1 = []
    temps = [] #first column is time, second is temperature
    mass = reactor.mass
    sdata = np.zeros([0, len(reactor.Y)])
    production_data = np.zeros([0, len(solution.net_production_rates)])
    state_list = list()
    
    #Prepare mass_fractions.hdf5 file.
    f1 = h5py.File('mass_fractions.hdf5', 'a')
    try:
        group_name = str(initial_temperature) + '_' + str(pressure) + '_' + str(frac)
    except ValueError:
        print "Duplicate initial conditions, or check to make sure mass fractions file isn't in directory. If it is, delete it"""
    individual = f1.create_group(group_name)
    
    #Run simulation 
    timer_start =tm.time()
    while current_time < stop_time:
        group_index += 1
        try:
            current_time = simulation.step()
        except Exception:
            error_string = 'Cantera autoignition_error @ %sK initial temperature' %initial_temperature
            print error_string
            return
        #Store information at this timestep
        times1.append(current_time)
        temps.append(reactor.T)
        species_data = reactor.Y
        grp = individual.create_group(str(group_index))
        grp['Temp'] = reactor.T
        grp['Time'] = current_time
        grp['Pressure'] = reactor.thermo.P
        species_production_rates = reactor.kinetics.net_production_rates
	try:
            net_rates_of_progress = reactor.kinetics.net_rates_of_progress
	except Exception:
	    return 0
        net_rates_of_progress = reactor.kinetics.net_rates_of_progress
        grp.create_dataset('Species Mass Fractions', data=species_data)
        grp.create_dataset('Reaction Rates of Progress', data=net_rates_of_progress)
        grp.create_dataset('Species Net Production Rates Original', data=species_production_rates)

        species_data = species_data[:, np.newaxis].T #translate from [n, 1] to [1,n]
        sdata = np.vstack((sdata, species_data))

        production_rates = np.array(solution.net_production_rates)
        production_rates = production_rates[:, np.newaxis].T
        production_data = np.vstack((production_data, production_rates))
    
    #Organize information collected from the simulation
    sample = get_range(times1, temps, sdata, production_data)
    timer_stop = tm.time()


    #strips all data except that within a 40 point sample range around ignition
    for grp in f1[group_name].keys():
        if int(grp) not in sample.index:
            f1[group_name].__delitem__(str(grp))


    #utility functions. Currently not active.  Work to be done here.  
    def plot(i):
        import matplotlib.pyplot as plt
        plt.clf()
        #plot combustion point
        #plt.plot(sample.derivative_max[0], sample.derivative_max[1], 'ro', ms=7, label='ignition point')
        #plot initial and final sample points
        #plt.plot(sample.initial_point[0], sample.initial_point[1], 'rx', ms=5, mew=2)
        #plt.plot(sample.final_point[0], sample.final_point[1], 'rx', ms=5, mew=2)
        #plot temp vs time
        plt.plot(sample.times_total, sample.temps_total)
        plt.xlabel('Time (s)')
        plt.title('Mixture Temperature vs Time')
        plt.ylabel('Temperature (K)')
        plt.axis([sample.times_total[0], sample.tau * 2, sample.temps_total[0] - 200, sample.temps_total[len(sample.temps_total) - 1] + 200])
        #plt.legend()
        plt.savefig("fig" + "_ic" + str(i) + ".png", bbox_inches='tight')
	plt.close()

    def writecsv(sdata, i):
        names = str(solution.species_names)
        tt = ['Time (s)', 'Temp (K)']
        names = solution.species_names
        name_array = np.append(tt, names)
        #sdata = sdata.astype('|S10')
        file_data = np.vstack((name_array, sdata))
        #open and write to file
        input_file_name_stripped = os.path.splitext(sys_args.data_file)[0]
        output_file_name = os.path.join(input_file_name_stripped + '_species_data' + 'ic_' + str(i) + '.csv')
        with open(output_file_name, 'wb') as f:
            np.savetxt(f, file_data, fmt=('%+12s'), delimiter=',')
        #os.system('atom '+ output_file_name)

    def writehdf5(sdata,i):
        
        input_file_name_stripped = os.path.splitext(sys_args.data_file)[0]
        output_file_name = os.path.join(input_file_name_stripped + '_ic_' + str(i) + '.csv')
        if not (os.path.exists("./hdf5_files")):
	    os.system("mkdir hdf5_files")
	os.system("cp mass_fractions.hdf5 mass_fractions_" + output_file_name)
        os.system("mv mass_fractions_" + output_file_name + " ./hdf5_files")
	os.system("cp production_rates.hdf5 production_rates_" + output_file_name)
        os.system("mv production_rates_" + output_file_name + " ./hdf5_files")
        #format matrix for hdf5
        #names = str(solution.species_names)
        #tt = ['Time (s)', 'Temp (K)']
        #names = solution.species_names
        #name_array = np.append(tt, names)
        #sdata = sdata.astype('|S10')
        #file_data = np.vstack((name_array, sdata))
        #open and write to file
        #input_file_name_stripped = os.path.splitext(sys_args.data_file)[0]
        #output_file_name = os.path.join(input_file_name_stripped + '_species_data.hdf5')
        #with h5py.File(output_file_name, 'w') as f:
            #Times = f.create_dataset("Times", data=sample.times_total)
            #Temps = f.create_dataset("Temps", data=sample.temps_total)
            #sgroup = f.create_group('Species_Data')
            #for i, sp in enumerate(solution.species_names):
                #sgroup.create_dataset(sp, data=sdata[i])

    def write_ai_times():
        f = open("autoignition_times.txt", "a")
        f.write(str(sample.temps_total[0]) + ", " + str(sample.tau) + "\n")
        f.close()

    def points():
        print("\nTime[s]            Temp[K]        Index        Point")
        print(str(sample.initial_point[0]) +  "       " + str("{0:.2f}".format(sample.initial_point[1]))\
         + "       " + str(sample.initial_point[2]) + "     " + "Initial sample point")
        print(str(sample.tau) + "        " + str("{0:.2f}".format(sample.derivative_max[1])) + "       " + str(sample.derivative_max[2])\
                + "     " + "Ignition point")
        print(str(sample.final_point[0]) +  "       " + str("{0:.2f}".format(sample.final_point[1]))\
         + "       " + str(sample.final_point[2]) + "     " + "Final sample point")

    #terminal use case
    if sys_args is not 'none' and info:
        if sys_args.plot:
            plot(i)
        if sys_args.writecsv:
            writecsv(sample.species_data,i)
        if sys_args.writehdf5:
            writehdf5(sample.species_data,i)
        if sys_args.points:
            points()
        if sys_args.write_ai_times:
            write_ai_times()
    
    #Create and return an object that contains criticial information about the simulation.
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
            self.tau_array = []
    return return_obj(sample.times, sample.temps, sample.species_data, f1, sample.tau, initial_temperature, frac)

