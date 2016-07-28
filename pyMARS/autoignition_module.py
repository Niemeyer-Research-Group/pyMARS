import sys
import numpy as np
import cantera as ct
import os
#import progressbar as pb
import h5py




def run_sim(mech_file, sys_args='none', **usr_args ):
    """
    Function to run Cantera reactor simulation

    Arguments
        Cantera mechanism file
        User arguments (plot='y', points='y', writehdf5='y', writecsv='y' )
    ----------
    Output
        Plot of Temp vs Time
        Points of interest
        CSV file
        Hdf5 file
    ----------
    Example
        run_sim('gri30.cti', points='y', plot='y')

    """

    data_file=mech_file
    solution1 = ct.Solution(mech_file)
    frac = raw_input('Enter mass fractions (ex.H2:2,O2:1,N2:4) :  ') # (ex.CH4:1, O2:2, N2:7.52 for Gri30 Stoich)and 100 for TPY where Y is mass fractions
    Temp= float(raw_input('Enter Solution Temperature in (K):'))
    solution1.TPX = Temp, ct.one_atm, frac #1001.0
    species=solution1.species()
    reactions=solution1.reactions()

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
    tnow=0.0
    tfinal=5.0e-3
    index1 = 0
    times1 = []
    temps = [] #first column is time, second is temperature
    mass=r1.mass
    sdata=np.zeros([0, len(r1.Y)])
    production_data=np.zeros([0, len(solution1.net_production_rates)])
    state_list=list()


    while tnow < tfinal:

        index1 += 1
        tnow = sim1.step(tfinal)
        times1.append(tnow)
        temps.append(r1.T)
        species_data=np.array(r1.Y) #*mass (optional)
        species_data=species_data[:,np.newaxis].T #translate from [n, 1] to [1,n]
        sdata= np.vstack((sdata, species_data))

        production_rates=np.array(solution1.net_production_rates)
        production_rates=production_rates[:,np.newaxis].T
        production_data=np.vstack((production_data, production_rates))

        """
        st=state(tnow, species, reactions)
        st.index = index1
        state_list.append(st)
        """


    print('\n')
    """-------------------------------------------------------------------------
    concatenate time, temp and coeff values
    -------------------------------------------------------------------------"""
    times1=np.array(times1)
    temps=np.array(temps)
    timetemp=np.vstack((times1, temps)).T
    sdata= np.hstack((timetemp, sdata))
    production_data= np.hstack((timetemp, production_data))

    """
    base=state_list.__getitem__(0)

    #concatenate coefficient values over time
    for idex, instance in enumerate(state_list):
        if idex > 0:                                                                #skip the base instance
            for spx in dir(instance):                                               #iterate through attributes (spx is a string)
                if spx.startswith('sp_'):                                           #select species attributes
                    base_spx = getattr(base, spx)
                    new_spx = getattr(instance, spx)
                    for ky in new_spx.keys():                                       #iterate through reaction keys
                        if ky in base_spx:                                          #if current reaction is in base reaction
                            cs=np.array(base_spx[ky])                               #convert to numpy array for easy appending
                            cs= np.append(cs, new_spx[ky])                          #append new coeff value to aray
                            base_spx[ky] = cs                                       #update reaction with new array of coefficient values

    #test length of concatenated matrices
    for spx in dir(base):
        if spx.startswith('sp_'):
            obj= getattr(base, spx)
            for att in obj.values():
                assert len(att) == len(times1), 'Number of coefficients does not match number of timesteps'
    """

    #get ignition point from dT/dt
    T=np.array(temps)
    dt= np.ones(len(times1)-1)*(times1[1]-times1[0])
    dT= np.diff(T)
    deriv= dT/dt
    i=deriv.argmax()
    deriv_max=[times1[i], T[i], i]
    tau= times1[i]
    print tau

    """-------------------------------------------------------------------------
    find initial and final sample points
    -------------------------------------------------------------------------"""

    initial_point=[times1[i-20], T[i-20], i-20]
    final_point=[times1[i+20], T[i+20], i+20]


    """-------------------------------------------------------------------------
    remove unnecessary data points (slice)
    -------------------------------------------------------------------------"""
    times_total=times1
    temps_total=temps
    sdata_total=sdata
    production_data_total=production_data


    times1=times1[ i-20:i+20 ]
    temps=temps[i-20:i+20 ]
    sdata=sdata[i-20:i+20, :]
    production_data=production_data[i-20:i+20, :]

    """
    for spx in dir(base):
        if spx.startswith('sp_'):
            obj= getattr(base, spx)
            for att in obj.values():
                att=att[i-20:i+20]
    """


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
        plt.ylabel('Temperature (K)')
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
        os.system('atom '+ output_file_name)

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
        def __init__(self, time, temp, sp_data):
            self.time=time
            self.temp=temp
            self.sp_data=sp_data

    return return_obj(times1, temps, sdata)

    "sdata is an array of 40 timesteps, with each instance containing an array of species"
    "mass fractions at that instant"









    """
    #find initial sample point
    for j, dTi in enumerate(dT):
        if dTi > .2:     #when dT > 5 degrees kelvin
            initial_point=[times1[j], T[j], j] # initial point, Temp, index
            break
    try:
        initial_point
    except NameError:
        try:
            #sys.exit("initial sample point not found") #alternative sys exit option
            initial_point=[times1[i-20], T[i-20], j]
            print("\nInitial Sample Point not found based on dTmax. Alternative method"+\
            " of 5 steps before tau is used.")
        except IndexError:
            try:
                initial_point=[times1[i], T[i], j]
            except IndexError:
                print 'Error: Initial sample point cannot be located'
                return

    #find final sample point
    for k, dti in enumerate(dT):
        if k > i:
            final_point=[times1[k+20], T[k+20], k+20]    #final point, Temp, index
            break
    if dti < .1:
        final_point=[times1[k], T[k], k]    #final point, Temp, index
        break
    initial_point[2] = i-20
    print initial_point[2]
    #final_point[2] = i+20
    print final_point[2]

    """
