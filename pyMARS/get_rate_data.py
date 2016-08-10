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
        mass_fractions=np.array(group['Species Mass Fractions'])
        #initalize solution, and set state
        solution=ct.Solution(mechanism_file)
        solution.Y=mass_fractions

        #create groups and write production data to new file
        production_rates=solution.net_production_rates
        new_grp=g.create_group(str(grp))
        new_grp['Temp'] = solution.T
        new_grp['Time'] = time
        new_grp.create_dataset('Reaction Production Rates', data=production_rates)
    g.close()

get_rates('mass_fractions.hdf5', 'gri30.cti')
