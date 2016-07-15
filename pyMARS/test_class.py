import cantera as ct
import h5py
import numpy as np

f=h5py.File('test_file.hdf5', 'w')

solution=ct.Solution('gri301.cti')
reactions=solution.reactions()
species=solution.species()
time = [3, 4, 5]


#sp=reactions[0].products.keys()
#.extend(reactions[0].reactants.keys())


"""
for spc in species:
    sp=f.create_group(spc.name)
    for i, rxn in enumerate(reactions):
        name= 'Reaction ' + str(i)
        local_products=rxn.products.keys()
        local_reactants=rxn.reactants.keys()
        local_total=rxn.products.keys()
        local_total.extend(rxn.reactants.keys())
        if spc.name in local_total:
            coeff=[]
            for j, t in enumerate(time):
                try:
                    coeff.append(rxn.reactants.get(spc.name))
                    continue
                except TypeError:
                    coeff.append(rxn.products.get(spc.name))
            sp.create_dataset(name, data=coeff)

print np.array(f.get('AR').get("Reaction 142"))
"""

"""
if spc.name in local_products
    coeff.append(rxn.products.get(spc.name))
if spc.name in local_reactants:
    coeff.append(rxn.reactants.get(spc.name))
"""

#this builds a state constaining all coefficients at one instant in time
class state:
    def __init__(self, time, species_list, solution, reactions):
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



state_list=list()

for i, t in enumerate(time):
    st=state(t, species, solution, reactions)
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



#this is how to append to a key value
z=state_list.__getitem__(0).sp_AR

print z['Reaction 142']
#print dir(state_list.__getitem__(0))
#print dir(state_list.__getitem__(0).sp_O.keys())
#print dir(state_list.__getitem__(0))

"""
l=[]
for t in time:
    rxn=reactions[0]
    sp=rxn.products.keys()
    sp.extend(rxn.reactants.keys())
    print type(sp)
    for name in sp:
        pass
        name =[]
        #name.append(reactions[0])
        #l.append(reactions[0].products.values())
"""

"""
rxn={}

for i, rxns in enumerate(reactions):

    sp=rxns.products
    sp.update(rxns.reactants)

    rxn[i] = sp
    for spc in rxn[i]:
        rxn.setdefault(spc, rxn.get(spc)).append(1)
    if i is 0:
        #print spc
        print rxn[i]



#print rxn.get(0)

#print dir(rxn.get(0).get('O2'))



rxn={}
for i, rxns in enumerate(reactions):
    index= 'Reaction ' + str(i)
    sp=rxns.products
    sp.update(rxns.reactants)

    rxn[index] = sp



print rxn.get('Reaction 1')
"""


"""
rxn_list=f.create_group('Reactions')
for t in time:
    n= 'tstep ' + str(t)
    tstep = rxn_list.create_group(n)
    for i, k in enumerate(reactions):
        index= 'Reaction ' + str(i)
        sp=k.products
        sp.update(k.reactants)

        rxn = tstep.create_group(index)  #makes a group for each reaction
        for ind in sp:
            rxn.create_dataset(ind, data=np.array(sp.get(ind)))   #makes a dataset for each species in each reaction
        #f.create_dataset(name, data=)

#print type(rxn_list.get('tstep 0').get('Reaction 1'))
"""
