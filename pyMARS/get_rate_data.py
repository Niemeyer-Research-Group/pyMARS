import cantera as ct
import h5py
import numpy as np
import os

def get_rates(hdf5_file, solution_object):
    """ Takes mass_fractions hdf5 file of species mole fractions at a point in
        time, and initalizes a cantera solution object to get species production
        rates at that point

    Parameters
    ----------
    hdf5_file : str
        A hdf5 file containing time, temp and species mole fractions used to
        set solution state
    solution_object : obj
        A Cantera solution object used to get net production rates

    Returns
    -------
        production_rates.hdf5
            [initial condition]
                [timesteps]
                    [temperature array]
                    [time array]
                    [Reaction Production rates dataset]
    """
    #read in data file
    f = h5py.File(hdf5_file, 'r')
    #create file for production rates
    g = h5py.File('production_rates.hdf5', 'w')
    #initialize solution
    solution = solution_object
    #iterate through all initial conditions
    for grp in f.iterkeys():
        #get solution data at individual timestep
        ic_group = g.create_group(grp.title())
        """
        #iterate through all timesteps
        for tstep in ic_group.iterkeys():
            group = f[str(grp)][str(tstep)]
            time = group['Time'].value
            temp = group['Temp'].value
            pressure = group['Pressure'].value
            mass_fractions = np.array(group['Species Mass Fractions'])
            #set solution state
            solution.TPY = temp, pressure, mass_fractions

            reaction_production_rates = solution.net_rates_of_progress
            new_grp = g.create_group(str(grp)+'_'+str(tstep))
            new_grp['Temp'] = solution.T
            new_grp['Time'] = time
            new_grp.create_dataset('Reaction Production Rates', data=reaction_production_rates)
        """
        for tstep in f[grp].iterkeys():
            group = f[grp][tstep]
            time = group['Time'].value
            temp = group['Temp'].value
            pressure = group['Pressure'].value
            mass_fractions = np.array(group['Species Mass Fractions'])
            #set solution state and get rates of progress of reactions
            solution.TPY = temp, pressure, mass_fractions
            reaction_production_rates = solution.net_rates_of_progress

            #create new groups and datasets in production rates file
            new_grp = g.create_group(str(tstep))
            new_grp['Temp'] = solution.T
            new_grp['Time'] = time
            new_grp.create_dataset('Reaction Production Rates', data=reaction_production_rates)

    g.close()
    f.close()
