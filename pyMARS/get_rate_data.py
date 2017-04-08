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
                    [Reaction Production Rates Original]
                        ['H2'] = -3.47
                        ['CO2'] = 0
    """
    #read in data file
    f = h5py.File(hdf5_file, 'r')
    #create file for production rates
    g = h5py.File('production_rates.hdf5', 'w')
    #initialize solution
    old_solution = solution_object
    #iterate through all initial conditions
    for grp in f.iterkeys():
        #get solution data at individual timestep
        ic_group = g.create_group(grp.title())
        for tstep in f[grp].iterkeys():
            #--------------------------------
            #reading from mass fractions file
            #--------------------------------
            group = f[grp][tstep]
            time = group['Time'].value
            temp = group['Temp'].value
            pressure = group['Pressure'].value
            mass_fractions = np.array(group['Species Mass Fractions'])
            old_species_prod_rates=np.array(group['Species Net Production Rates Original'])
            old_reaction_rates_of_progress=np.array(group['Reaction Rates of Progress'])
            #---------------------------------------------------------
            #set solution state and get rates of progress of reactions
            #---------------------------------------------------------
            new_solution = old_solution
            new_solution.TPY = temp, pressure, mass_fractions
            new_reaction_production_rates = new_solution.net_rates_of_progress
            new_species_prod_rates=new_solution.net_production_rates

            #----------------------------------------------------
            #check that reaction_rates_of_progress are same
            #----------------------------------------------------
            for i, n in enumerate(old_reaction_rates_of_progress):
                if abs(old_reaction_rates_of_progress[i] - new_reaction_production_rates[i]) > .0001:
                    print '-----------------------'
                    print 'reaction rate of progress is off'
                    print solution.reaction(i)

            #----------------------------------------------------
            #check that species net production rates are the same
            #----------------------------------------------------

            for i, n in enumerate(old_species_prod_rates):
                if abs(n - new_species_prod_rates[i]) > .0001:
                    print '------------------------------------------'
                    print 'species production rates did not match in: '
                    print ic_group


            #----------------------------------------------------
            #check that species mass fractions are the same
            #----------------------------------------------------
            for i, n in enumerate(mass_fractions):
                if abs(mass_fractions[i]-new_solution.Y[i]) > .0001:
                    print '------------------------------------------'
                    print 'species mass fractions did not match for: '
                    print solution.species(i).name

            #-------------------------------------------------------
            #create new groups and datasets in production rates file
            #-------------------------------------------------------
            new_grp = ic_group.create_group(str(tstep))
            new_grp['Temp'] = new_solution.T
            new_grp['Time'] = time

            species_production_list = {}
            for i, n in enumerate(new_solution.species()):
                species_production_list[new_solution.species(i).name] = new_solution.net_production_rates[i]
            sp_data = new_grp.create_group('Species Net Production Rates Original')
            for j in species_production_list:
                sp_data[str(j)] = species_production_list[str(j)]
            #new_grp['Species Net Production Rates Original'] = species_production_list
            new_grp.create_dataset('Reaction Production Rates', data=new_reaction_production_rates)

    g.close()
    f.close()
