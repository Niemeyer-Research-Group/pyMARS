import cantera as ct
import h5py
import numpy as np
import os

def get_rates(hdf5_file, solution_object, initial_temperature):
    """ Takes mass_fractions hdf5 file of species mole fractions at a point in
        time, and initalizes a cantera solution object to get species production
        rates at that point

    :param hdf5_file:
        A hdf5 file containing time, temp and species mole fractions used to
        set solution state

    :param solution_object:
        A Cantera solution object used to get net production rates

    :param initial_temperature:
        Initial temperature used in autoignition

    """

    #read in data file
    f=h5py.File(hdf5_file, 'r')

    #create file for production rates
    g=h5py.File('production_rates.hdf5', 'w')
    #initialize solution
    solution=solution_object
    ind = str(initial_temperature)
    #iterate through all initial conditions
    for grp in f.iterkeys():
        #get solution data at individual timestep
        ic_group = f[grp]
        print ic_group
        print 'rates converted'
        #iterate through all timesteps
        for tstep in ic_group.iterkeys():
            group=f[str(grp)][str(tstep)]
            time=group['Time'].value
            temp=group['Temp'].value
            pressure=group['Pressure'].value
            mass_fractions=np.array(group['Species Mass Fractions'])
            #set solution state
            solution.TPY= temp, pressure, mass_fractions

            #create groups and write production data to new file
            species_production_rates=solution.net_production_rates
            species_production_rates_original=group['Species Net Production Rates Original']
            reaction_production_rates=solution.net_rates_of_progress
            new_grp=g.create_group(str(grp)+'_'+str(tstep))
            new_grp['Temp'] = solution.T
            new_grp['Time'] = time
            #net production rates for each species (kmol/m^3/s for bulk phases, or kmol/m^2/s for surface)
            new_grp.create_dataset('Species Net Production Rates', data=species_production_rates)
            #net rates of progress for each reaction  (kmol/m^3/s for bulk phases, or kmol/m^2/s for surface)
            new_grp.create_dataset('Reaction Production Rates', data=reaction_production_rates)
            new_grp.create_dataset('Species Net Production Rates Original', data=species_production_rates_original)
    g.close()
    f.close()
    os.system('rm mass_fractions.hdf5')
