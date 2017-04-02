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
                    [Species Net Production Rates Original]
                        ['H2'] = -3.47
                        ['CO2'] = 0
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
            #reading from mass fractions file
            group = f[grp][tstep]
            time = group['Time'].value
            temp = group['Temp'].value
            pressure = group['Pressure'].value
            mass_fractions = np.array(group['Species Mass Fractions'])
            original_species_prod_rates=np.array(group['Species Net Production Rates Original'])
            #set solution state and get rates of progress of reactions
            solution.TPY = temp, pressure, mass_fractions
            reaction_production_rates = solution.net_rates_of_progress

            #check that species net production rates are the same
            new_species_prod_rates=solution.net_production_rates
            for i, n in enumerate(original_species_prod_rates):
                #if new_species_prod_rates[i] != n:
                #    print ('old %0.7f, new %0.7f') %(n, new_species_prod_rates[i])
                assert abs(n - new_species_prod_rates[i]) <= .001
            #assert abs(original_species_prod_rates == new_species_prod_rates


            #create new groups and datasets in production rates file
            new_grp = ic_group.create_group(str(tstep))
            new_grp['Temp'] = solution.T
            new_grp['Time'] = time

            species_production_list = {}
            for i, n in enumerate(solution.species()):
                species_production_list[solution.species(i).name] = solution.net_production_rates[i]
            sp_data = new_grp.create_group('Species Net Production Rates Original')
            for j in species_production_list:
                sp_data[str(j)] = species_production_list[str(j)]
            #new_grp['Species Net Production Rates Original'] = species_production_list
            new_grp.create_dataset('Reaction Production Rates', data=reaction_production_rates)

    g.close()
    f.close()
