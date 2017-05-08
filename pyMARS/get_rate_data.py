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
    total_edge_data = {}
    for grp in f.iterkeys():
        #get solution data at individual timestep
        ic_group = g.create_group(grp.title())
        ic_edge_data = {}
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
            #
            #for i, n in enumerate(old_species_prod_rates):
            #    if abs(n - new_species_prod_rates[i]) > .0001:
            #        print '------------------------------------------'
            #        print 'species production rates did not match in: '
            #        print ic_group


            #----------------------------------------------------
            #check that species mass fractions are the same
            #----------------------------------------------------
            #for i, n in enumerate(mass_fractions):
            #    if abs(mass_fractions[i]-new_solution.Y[i]) > .0001:
            #        print '------------------------------------------'
            #        print 'species mass fractions did not match for: '
            #        print solution.species(i).name

            #-------------------------------------------------------
            #create new groups and datasets in production rates file
            #-------------------------------------------------------
            new_grp = ic_group.create_group(str(tstep))
            new_grp['Temp'] = new_solution.T
            new_grp['Time'] = time
            sp_data = new_grp.create_group('Species Net Production Rates Original')

            #generate list of net species production rates
            species_production_list = {}
            for i, n in enumerate(new_solution.species()):
                species_production_list[new_solution.species(i).name] = new_species_prod_rates[i]

            #write list of net species production rates to hdf5 object
            for j in species_production_list:
                sp_data[str(j)] = species_production_list[str(j)]
            #new_grp['Species Net Production Rates Original'] = species_production_list

            #write reaction production rates to hdf5 object
            new_grp.create_dataset('Reaction Production Rates', data=new_reaction_production_rates)


            #writes reaction ID and products/reactants with coeffs to hdf5 object
            rxn_groups = new_grp.create_group('Reactions') #= str([rx.ID for rx in new_solution.reactions()])
            for rx in new_solution.reactions():
                rxi = rxn_groups.create_group(rx.ID)
                rxiReactants = rxi.create_group('Reactants')
                rxiProducts = rxi.create_group('Products')
                for spi in rx.products:
                    rxiProducts[spi] = rx.products[spi]
                for spi in rx.reactants:
                    rxiReactants[spi] = rx.reactants[spi]

            #some extra stuff that is going to slow everything down, but I need right now for troubleshooting
            #checking method for calculating species production rates
            #list_A = species_production_list
            list_A = sp_data
            list_B = {}
            list_C = {}
            for spc in new_solution.species():
                list_B[spc.name] = 0
            #calculate species production rates as in the DRG. I've proven this method
            #to work in a few other test functions
            for i, reac in enumerate(new_solution.reactions()):
                reac_prod_rate = float(new_reaction_production_rates[i])
                reactants = reac.reactants
                products = reac.products
                all_species = reac.reactants
                all_species.update(reac.products)
                if reac_prod_rate != 0:
                    if reac_prod_rate > 0:
                        for species in products:
                            list_B[species] += abs(float(reac_prod_rate*products[species]))

                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in list_C:
                                        list_C[partial_name] += abs(float(reac_prod_rate*products[species]))
                                    else:
                                        list_C[partial_name] = abs(float(reac_prod_rate*products[species]))
                        for species in reactants:
                            list_B[species] += abs(float(-reac_prod_rate*reactants[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in list_C:
                                        list_C[partial_name] += abs(float(-reac_prod_rate*reactants[species]))
                                    else:
                                        list_C[partial_name] = abs(float(-reac_prod_rate*reactants[species]))
                    if reac_prod_rate < 0:
                        for species in products:
                            list_B[species] += float(reac_prod_rate*products[species])
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in list_C:
                                        list_C[partial_name] += abs(float(reac_prod_rate*products[species]))
                                    else:
                                        list_C[partial_name] = abs(float(reac_prod_rate*products[species]))
                        for species in reactants:
                            list_B[species] += float(-reac_prod_rate*reactants[species])

                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in list_C:
                                        list_C[partial_name] += abs(float(-reac_prod_rate*reactants[species]))
                                    else:
                                        list_C[partial_name] = abs(float(-reac_prod_rate*reactants[species]))
            # for blah in list_A:
            #     if abs(float(list_A[blah].value) - np.float64(list_B[blah])) > .0001:
            #         print '--------//------'
            #         print blah
            #         print list_A[blah] - list_B[blah]
            #         print '--------//-------'
            #         print '---------------'

            ic_edge_data[tstep] = [list_A, list_B, list_C]
        total_edge_data[grp.title()] = ic_edge_data
    #print total_edge_data['1400.0_101325.0_N2:41.36,O2:11.0,Nc7H16:1.0'].keys()
    # print total_edge_data['1400.0_101325.0_N2:41.36,O2:11.0,Nc7H16:1.0']['1423']
    # print len(total_edge_data['1400.0_101325.0_N2:41.36,O2:11.0,Nc7H16:1.0']['1423'])
    return total_edge_data



    g.close()
    f.close()
