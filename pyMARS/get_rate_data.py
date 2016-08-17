import cantera as ct
import h5py
import numpy as np

def get_rates(hdf5_file, mechanism_file):
    #read in data file
    f=h5py.File(hdf5_file, 'r')

    #create file for production rates
    g=h5py.File('production_rates.hdf5', 'w')
    #iterate through all 40 timesteps
    for grp in f:
        #get solution data at individual timestep
        group=f[str(grp)]
        time=group['Time'].value
        temp=group['Temp'].value
        pressure=group['Pressure'].value
        mass_fractions=np.array(group['Species Mass Fractions'])
        #initalize solution, and set state
        solution=ct.Solution(mechanism_file)
        solution.TPY= temp, pressure, mass_fractions

        #create groups and write production data to new file
        species_production_rates=solution.net_production_rates
        species_production_rates_original=group['Species Net Production Rates Original']
        reaction_production_rates=solution.net_rates_of_progress
        new_grp=g.create_group(str(grp))
        new_grp['Temp'] = solution.T
        new_grp['Time'] = time
        #net production rates for each species (kmol/m^3/s for bulk phases, or kmol/m^2/s for surface)
        new_grp.create_dataset('Species Net Production Rates', data=species_production_rates)
        #net rates of progress for each reaction  (kmol/m^3/s for bulk phases, or kmol/m^2/s for surface)
        new_grp.create_dataset('Reaction Production Rates', data=reaction_production_rates)
        new_grp.create_dataset('Species Net Production Rates Original', data=species_production_rates_original)
    g.close()

#get_rates('mass_fractions.hdf5', 'gri301.cti')
