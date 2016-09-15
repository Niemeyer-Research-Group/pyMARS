import sys
import numpy as np
import cantera as ct
import os
#import progressbar as pb
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

    solution1 = solution_object
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
        frac = raw_input('Enter mole fractions (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich) :  ') # (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich)and 100 for TPY where Y is mass fractions
        initial_temperature = float(raw_input('Enter Solution Temperature in (K):'))
        initial_temperature = initial_temperature#+273.0 #convert to kelvin
    else:
        frac = sys_args.frac
        initial_temperature = sys_args.Temp
    solution1.TPX = initial_temperature, ct.one_atm, frac #1001.0
    species = solution1.species()
    reactions = solution1.reactions()

    class state:
        def __init__(self, time, species_list, reactions):
            self.time = time
            for sp in species_list:
                coeff={}
                for i, rxn in enumerate(reactions):
                    name = 'Reaction ' + str(i)
                    if sp.name in rxn.reactants.keys():
                        coeff[name] = (rxn.reactants.get(sp.name))
                    if sp.name in rxn.products.keys():
                        coeff[name] = (rxn.products.get(sp.name))

                setattr(self, 'sp_'+ sp.name, coeff)



    """-------------------------------------------------------------------------
    run sim to find ignition delay from dT/dt max
    -------------------------------------------------------------------------"""

    r1 = ct.Reactor(solution1)
    sim1 = ct.ReactorNet([r1])
    tnow = 0.0
    tfinal = 5.0e-3
    index1 = 0
    times1 = []
    temps = [] #first column is time, second is temperature
    mass=r1.mass
    sdata=np.zeros([0, len(r1.Y)])
    production_data=np.zeros([0, len(solution1.net_production_rates)])
    state_list=list()

    f1=h5py.File('mass_fractions.hdf5', 'w')
    while tnow < tfinal:
        index1 += 1
        tnow = sim1.step(tfinal)
        times1.append(tnow)
        temps.append(r1.T)
        species_data=r1.Y

        grp=f1.create_group(str(index1))
        grp['Temp'] = r1.T
        grp['Time'] = tnow
        grp['Pressure'] =r1.thermo.P
        species_production_rates=r1.thermo.net_production_rates
        grp.create_dataset('Species Mass Fractions', data=species_data)
        grp.create_dataset('Species Net Production Rates Original', data=species_production_rates)

        species_data=species_data[:,np.newaxis].T #translate from [n, 1] to [1,n]
        sdata= np.vstack((sdata, species_data))

        production_rates=np.array(solution1.net_production_rates)
        production_rates=production_rates[:,np.newaxis].T
        production_data=np.vstack((production_data, production_rates))


    print('\n')
    """-------------------------------------------------------------------------
    concatenate time, temp and coeff values
    -------------------------------------------------------------------------"""
    times1=np.array(times1)
    temps=np.array(temps)
    timetemp=np.vstack((times1, temps)).T
    sdata= np.hstack((timetemp, sdata))
    production_data= np.hstack((timetemp, production_data))


    #get ignition point from dT/dt
    T=np.array(temps)
    dt= np.ones(len(times1)-1)*(times1[1]-times1[0])
    dT= np.diff(T)
    deriv= dT/dt
    i=deriv.argmax()
    deriv_max=[times1[i], T[i], i]
    tau= times1[i]

    """-------------------------------------------------------------------------
    find initial and final sample points
    -------------------------------------------------------------------------"""
    try:
        initial_point=[times1[i-20], T[i-20], i-20]
    except IndexError:
        initial_point=[times1[0], T[0], 0]
        print 'not enough timesteps before ignition'
        print 'timesteps before ignition: %s' %i
        print 'total timesteps: %s' % len(T)
    try:
        final_point=[times1[i+20], T[i+20], i+20]
    except IndexError:
        final_point=[times1[len(times1)-1], T[len(T)-1], (len(T)-1)]
        print 'not enough timesteps after ignition'
        print 'timesteps after ignition: %s' %(len(T)-i)
        print 'total timesteps: %s' % len(T)


    """-------------------------------------------------------------------------
    remove unnecessary data points (slice)
    -------------------------------------------------------------------------"""
    times_total=times1
    temps_total=temps
    sdata_total=sdata
    production_data_total=production_data

    for grp in f1.keys():
        if int(grp) not in range((i-20), (i+20)):
            f1.__delitem__(str(grp))

    f1.close()

    times1=times1[ i-20:i+20 ]
    temps=temps[i-20:i+20 ]
    sdata=sdata[i-20:i+20, :]
    production_data=production_data[i-20:i+20, :]


    """-------------------------------------------------------------------------
    utility functions
    -------------------------------------------------------------------------"""


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
        names=str(solution1.species_names)
        tt=['Time (s)', 'Temp (K)']
        names = solution1.species_names
        name_array = np.append(tt, names)
        sdata=sdata.astype('|S10')
        file_data= np.vstack((name_array, sdata))
        #open and write to file
        input_file_name_stripped=os.path.splitext(data_file)[0]
        output_file_name=os.path.join(input_file_name_stripped + '_species_data.csv')
        print output_file_name
        with open(output_file_name, 'wb') as f:
            np.savetxt(f, file_data, fmt=('%+12s'),  delimiter=',')
        #os.system('atom '+ output_file_name)

    def writehdf5(sdata):
        #format matrix for hdf5
        names=str(solution1.species_names)
        tt=['Time (s)', 'Temp (K)']
        names = solution1.species_names
        name_array = np.append(tt, names)
        sdata=sdata.astype('|S10')
        file_data= np.vstack((name_array, sdata))

        #open and write to file
        input_file_name_stripped=os.path.splitext(data_file)[0]
        output_file_name=os.path.join(input_file_name_stripped + '_species_data.hdf5')
        with h5py.File(output_file_name, 'w') as f:
            Times = f.create_dataset("Times", data=times1)
            Temps = f.create_dataset("Temps", data=temps)
            sgroup= f.create_group('Species_Data')

            for i, sp in enumerate(solution1.species_names):
                    sgroup.create_dataset(sp, data=sdata[:,i+2])

    def points():
        print("\nTime[s]            Temp[K]        Index        Point")
        print( str(initial_point[0]) +  "       " + str("{0:.2f}".format(initial_point[1]))\
         + "       " + str(initial_point[2]) + "     " + "Initial sample point")
        print(str(tau) + "        " + str("{0:.2f}".format(deriv_max[1])) + "       " + str(deriv_max[2])\
                + "     " + "Ignition point")
        print( str(final_point[0]) +  "       " + str("{0:.2f}".format(final_point[1]))\
         + "       " + str(final_point[2]) + "     " + "Final sample point")





    if sys_args is not 'none':
        if sys_args.plot:
            plot()
        if sys_args.writecsv:
            writecsv(sdata)
        if sys_args.writehdf5:
            writehdf5(sdata)
        if sys_args.points:
            points()
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
            self.time=time
            self.temp=temp
            self.sp_data=sp_data
            self.test=f1
            self.tau=tau
            self.Temp = initial_temperature
            self.frac = frac

    return return_obj(times1, temps, sdata, f1, tau, initial_temperature, frac)

    "sdata is an array of 40 timesteps, with each instance containing an array of species"
    "mass fractions at that instant"
