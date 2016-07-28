import cantera as ct
import h5py
import numpy as np

f=h5py.File('test_file.hdf5', 'w')

solution=ct.Solution('gri301.cti')
reactions=solution.reactions()
species=solution.species()
time = [3, 4, 5]



#this builds a state constaining all coefficients at one instant in time
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


#generate list of every instance in time, and append to it
state_list=list()
for i, t in enumerate(time):
    st=state(t, species,  reactions)
    st.index = i
    state_list.append(st)

base = state_list.__getitem__(0)    #base instance to be appended to


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



#this is how to get to an array of coefficient values
z=state_list.__getitem__(0).sp_O2['Reaction 30']

for spx in dir(base):
    if spx.startswith('sp_'):
        obj= getattr(base, spx)                                                 #species object
        for att in obj.keys():                                                 #iterate through reaction keys
            obj1=obj[att]

print base.sp_AR['Reaction 142']
print dir(base.sp_AR.keys())
