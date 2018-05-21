import networkx as nx
import numpy as np
import h5py
from collections import Counter
import time as tm
import os, sys, argparse
import cantera as ct
import soln2ck
import soln2cti
from readin_initial_conditions import readin_conditions
from simulation import Simulation
import helper
from create_trimmed_model import trim
from numpy import genfromtxt
import math
from dijkstra import ss_dijkstra_path_length_modified

os.environ['Cantera_Data'] =os.getcwd()

def make_dic_drgep(solution_object, total_edge_data, target_species):
    """ Use the Direct Relation Graph with Error Propegation (DRGEP) method to build a nodal graph of
        species and use the graph to determine each species overall interaction coefficents.

    Parameters
    ----------
    solution_object : obj
        Cantera Solution object
    total_edge_data : array
        A 3-D array containing information for calculating direct interaction coeffiecnets which will serve as edge weights on the graphs.
    target_species: list of strings
	a list containing the names of the target species.

    Returns
    -------
    max_dic: Dictionary
        A Dictionary keyed by species name with values that represent that species importance to the system.
    """
    #initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    graph = nx.DiGraph() #Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.

    max_dic = {} #Dictionary holding the maximum values for the iteration

    #calculate edge weights based on list received from get_rate_data and use them to create a graph
    for ic in total_edge_data.keys(): #For each initial condition
        for species in species_objects: #Make graph
            graph.add_node(species.name)
        #timestep
        for tstep in total_edge_data[ic].keys(): #Make a graph at each timestep
            numerator = total_edge_data[ic][tstep][1] #DRGEP calculations of direct interaction coeffients are done in the total_edge_data function.
            denominator = total_edge_data[ic][tstep][0]
            for edge in numerator: #For each edge, determine its weight amnd add it to the graph
                try:
                    edge_name = edge.split('_', 1)
                    species_a_name = edge_name[0]
                    species_b_name = edge_name[1]
                    #dgep weight between two species
                    if denominator[species_a_name] != 0:
                        weight = abs(float(numerator[edge])/float(denominator[species_a_name]))
                        if graph.has_edge(species_a_name, species_b_name):
                            old_weight = graph[species_a_name][species_b_name]['weight']
                            if weight > old_weight and weight <= 1:
                                graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
                            elif weight > 1:
                                print("Error.  Edge weights should not be greater than one.")
                                exit()
                        elif weight <= 1:
                            graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
                        elif weight > 1:
                            print("Error.  Edge weights should not be greater than one.")
                            exit()
                except IndexError:
                    print(edge)
                    continue

            #Search the graph for overall interaction coefficents and add them to max_dic if they belong
            dic = graph_search_drgep(graph, target_species) #Search graph for max values to each species based on targets
            for sp in dic: #Add to max dictionary if it is new or greater than the value already there.
                if sp not in max_dic:
                    max_dic[sp] = dic[sp]
                elif dic[sp] > max_dic[sp]:
                    max_dic[sp] = dic[sp]
            graph.clear() #Reset graph
    return max_dic


def trim_drgep(max_dic, solution_object, threshold_value, keeper_list, done):
    """ Use the dictionary created by the drgep method to determine what should
        be cut out of the model at a specific threshold value.

    Parameters
    ----------
    max_dic: Dictionary
        A dictionary keyed by species name that has values representing the species importance to the system.
    solution_object : obj
        Cantera Solution object
    threshold_value : int
        an edge weight threshold value
    keeper_species: list of strings
	a list of species for the mechinism to keep no matter what
    done: singleton
	a singleton boolean value that represents if their are more species to cut out or not.

    Returns
    -------
    exclusion_list : list
        List of species to trim from mechanism
    """


    core_species = []
    species_objects = solution_object.species()

    #Take all species that are over the threshold value and add them to essentail species.
    essential_species = []
    for sp in species_objects:
        if sp.name in max_dic:
            if max_dic[sp.name] > threshold_value and sp not in essential_species:
                essential_species.append(sp)
    done[0] = True
    for sp in species_objects: #If any more can be taken out, we are not done yet.
        if sp.name in max_dic:
            if max_dic[sp.name] > threshold_value:
                done[0] = False

    for sp in essential_species:
        if sp not in core_species:
            core_species.append(sp.name)

    retained_species = keeper_list #Specified by the user.  A list of species that also need to be kept.
    for sp in retained_species:
        if sp not in core_species:
            core_species.append(sp)

    exclusion_list = []

    for species in solution_object.species():
        if species.name not in core_species: #If its not one of our species we must keep, add it to the list of species to be trimmed.
            exclusion_list.append(species.name)

    return exclusion_list

def run_drgep(args, solution_object,error,past):
	"""
    Function to run the drgep method for model reduction

    Parameters
    ----------
    solution_object: ct._Solutuion object
	An object that represents the solution to be trimmed.
    args: args object
	An object that contains the following:
    file:
        Input mechanism file (ex. file='gri30.cti')
    species:
        Species to eliminate (ex. species='H, OH')
    thermo:
        Thermo data file if Chemkin format (ex. thermo= 'thermo.dat')
    transport:
        Transport data file if Chemkin format
    plot:
        plot ignition curve (ex. plot='y')
    points:
        print ignition point and sample range (ex. points='y')
    writecsv:
        write data to csv (ex. writecsv='y')
    writehdf5:
        write data to hdf5 (ex. writehdf5='y')
    run_drg:
        Run DRG model reduction
    run_drgep:
	Run drgep model.
    error:
	Maximum ammount of error allowed.
    keepers: list of strings
	The string names of the species that should be kept in the model no matter what.
    targets: list of strings
	The string names of the species that should be used as target species.

	"""

        #Set up variables
	if len(args.target) == 0: #If the target species are not specified, puke and die.
		print("Please specify a target species.")
		exit()
	done = [] #Singleton to hold wether or not any more species can be cut from the simulation.
	done.append(False)
	threshold = .1 #Starting threshold value
	threshold_i = .1
	n = 1
	error = [10.0] #Singleton to hold the error value of the previously ran simulation.

	#Set up simulations

	if args.conditions_file:
		conditions_array = readin_conditions(str(args.conditions_file))
	elif not args.conditions_file:
		print("Conditions file not found")
		exit()

	sim_array = helper.setup_simulations(conditions_array,solution_object) #Turn conditions array into unran simulation objects



    #Find overall interaction coefficients.
	original_tau = []
	sample_points = []
	for case in sim_array: #Run simulations for original model and process results
		original_tau.append(case.run_case())
		sample_points.append(case.process_results())

	ignition_delay_detailed = np.array(original_tau) #Turn tau array into a numpy array
	print(ignition_delay_detailed)


	rate_edge_data = get_rates(sim_array,solution_object) #Get edge weight calculation data.
	max_dic = make_dic_drgep(solution_object, rate_edge_data, args.target) #Make a dictionary of overall interaction coefficients.
	print(max_dic)

	print("Testing for starting threshold value")
	drgep_loop_control(solution_object, args, error, threshold, done, max_dic,ignition_delay_detailed,conditions_array) #Trim the solution at that threshold and find the error.
	while error[0] != 0: #While the error for trimming with that threshold value is greater than allowed.
		threshold = threshold / 10 #Reduce the starting threshold value and try again.
		threshold_i = threshold_i / 10
		n = n + 1
		drgep_loop_control(solution_object, args, error, threshold, done, max_dic,ignition_delay_detailed,conditions_array)
		if error[0] <= .02:
			error[0] = 0

	print("Starting with a threshold value of " + str(threshold))
	sol_new = solution_object
	past[0] = 0 #An integer representing the error introduced in the past simulation.
	done[0] = False

	while not done[0] and error[0] < args.error: #Run the simulation until nothing else can be cut.
		sol_new = drgep_loop_control( solution_object, args, error, threshold, done, max_dic,ignition_delay_detailed,conditions_array) #Trim at this threshold value and calculate error.
		if args.error > error[0]: #If a new max species cut without exceeding what is allowed is reached, save that threshold.
			max_t = threshold
                        #if (past == error[0]): #If error wasn't increased, increase the threshold at a higher rate.
		        #threshold = threshold + (threshold_i * 4)
			past[0] = error[0]
		        #if (threshold >= .01):
                        #threshold_i = .01
			threshold = threshold + threshold_i
			threshold = round(threshold, n)

	print("\nGreatest result: ")
	sol_new = drgep_loop_control( solution_object, args, error, max_t, done, max_dic,ignition_delay_detailed,conditions_array)
	#if os.path.exists("production_rates.hdf5"):
	#	os.system('rm production_rates.hdf5')
	#if os.path.exists('mass_fractions.hdf5'):
	#	os.system('rm mass_fractions.hdf5')
	drgep_trimmed_file = soln2cti.write(sol_new) #Write the solution object with the greatest error that isn't over the allowed ammount.
	return sol_new[1]
	#try:
	#	os.system('rm production_rates.hdf5')
	#except Exception:
	#	pass


def drgep_loop_control(solution_object, args, stored_error, threshold, done, max_dic,ignition_delay_detailed,conditions_array):
    """ Controls the trimming and error calculation of a solution with the graph already created using the DRGEP method.

        Parameters
        ----------
        solution_object : obj
            Cantera solution object
        args : obj
            function arguments object
	stored_error: float singleton
	    The error introduced by the last simulation (to be replaced with this simulation).
	done: singleton
	    a singleton boolean value that represnts wether or not more species can be excluded from the graph or not.
	max_dic: dictionary
	    a dictionary keyed by species name that represents the species importance to the model.

        Returns
        -------
        new_solution_objects : obj
            Cantera solution object That has been reduced.
        """

    target_species = args.target

    #run detailed mechanism and retain initial conditions
    species_retained = []
    printout = ''
    print('Threshold     Species in Mech      Error')

    #run DRGEP and create new reduced solution
    exclusion_list = trim_drgep(max_dic, solution_object, threshold, args.keepers, done) #Find out what to cut from the model
    new_solution_objects = trim(solution_object, exclusion_list, args.data_file) #Cut the exclusion list from the model.
    species_retained.append(len(new_solution_objects[1].species()))

    #simulated reduced solution
    new_sim = helper.setup_simulations(conditions_array,new_solution_objects[1])
    new_tau = []
    new_sample_points = []
    for case in new_sim: #Run simulations for original model and process results
        new_tau.append(case.run_case())
        new_sample_points.append(case.process_results())

    ignition_delay_reduced = np.array(new_tau) #Turn tau array into a numpy array
    if (ignition_delay_detailed.all() == 0):
        print("Original model did not ignite.  Check initial conditions.")
        exit()

    #Calculate and print error.
    error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100 #Calculate error
    printout += str(threshold) + '                 ' + str(len(new_solution_objects[1].species())) + '              '+  str(round(np.max(error), 2)) +'%' + '\n'
    print(printout)
    stored_error[0] = round(np.max(error), 2)

    #Return new model
    new_solution_objects = new_solution_objects[1]
    return new_solution_objects

def get_rates(sim_array, solution_object):
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
                [timesteps 1x40]
                    [temperature 1x1]
                    [time 1x1]
                    [Reaction Production Rates Original]
                        ['H2'] = -3.47
                        ['CO2'] = 0
    """


    #initialize solution
    old_solution = solution_object
    #iterate through all initial conditions
    total_edge_data = {}
    for ic in sim_array:
        ic_edge_data = {}
        for tstep in ic.sample_points:
            print(tstep)
            temp = tstep[0]
            pressure = tstep[1]
            mass_fractions = np.array(tstep[2])

            new_solution = old_solution
            new_solution.TPY = temp, pressure, mass_fractions
            new_reaction_production_rates = new_solution.net_rates_of_progress
            new_species_prod_rates=new_solution.net_production_rates

            denom = {}
            numerator = {}
            for spc in new_solution.species():
                denom[spc.name] = []
                denom[spc.name].append(0)
                denom[spc.name].append(0)
            #calculate species interaction coefficeints as in DRGEP. I've proven this method
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
                            denom[species][1] += abs(float(reac_prod_rate*products[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator:
                                        numerator[partial_name] += float(reac_prod_rate*products[species])
                                    else:
                                        numerator[partial_name] = float(reac_prod_rate*products[species])

                        for species in reactants:
                            denom[species][0] += abs(float(reac_prod_rate*reactants[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator:
                                        numerator[partial_name] += float(-reac_prod_rate*reactants[species])
                                    else:
                                        numerator[partial_name] = float(-reac_prod_rate*reactants[species])

                    if reac_prod_rate < 0:
                        for species in products:
                            denom[species][0] += abs(float(reac_prod_rate*products[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator:
                                         numerator[partial_name] += float(reac_prod_rate*products[species])
                                    else:
                                         numerator[partial_name] = float(reac_prod_rate*products[species])

                        for species in reactants:
                            denom[species][1] += abs(float(reac_prod_rate*reactants[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator:
                                         numerator[partial_name] += float(-reac_prod_rate*reactants[species])
                                    else:
                                         numerator[partial_name] = float(-reac_prod_rate*reactants[species])

            for species in new_solution.species():
                if abs(denom[species.name][0]) > abs(denom[species.name][1]):
                    denom[species.name] = abs(denom[species.name][0])
                else:
                    denom[species.name] = abs(denom[species.name][1])

            for name in numerator:
                numerator[name] = abs(numerator[name])
            ic_edge_data[temp] = [denom, numerator]
        total_edge_data[ic] = ic_edge_data

    return total_edge_data

def graph_search_drgep(nx_graph, target_species):
    """Search nodal graph and generate a dictionary of the greatest paths to all species from one of the targets.

    Parameters
    ----------
    nx_graph : obj
        networkx graph object of solution\
    target_species : list
        List of target species to search from

    Returns
    -------
    max_dic : dictionary
        Values of the greatest possible path to each species from one of the targets on this graph keyed by species name.
    """

    max_dic = {} #A dictionary holding the maximum path to each species.
    for target in target_species:
        dic = ss_dijkstra_path_length_modified(nx_graph, target) #Get dictionary of each values maximum path
        for sp in dic:    #If the species are not in the max dictionary or the new value for that species is greater than the one in the max dictionary, add it to the max dictionary.
            if sp not in max_dic:
                max_dic[sp] = dic[sp]
            elif max_dic[sp] < dic[sp]:
                max_dic[sp] = dic[sp]
    return max_dic
