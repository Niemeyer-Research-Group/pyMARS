import networkx as nx
import numpy as np
import h5py
from collections import Counter
import time as tm
from drg import graph_search
import os, sys, argparse
import cantera as ct
import soln2ck
import soln2cti
import math
from create_trimmed_model import trim
from numpy import genfromtxt
from readin_initial_conditions import readin_conditions
import helper

os.environ['Cantera_Data'] =os.getcwd()


#################
# This function determines what species should be excluded from the reduced model based on their DICs compared to the threshold value and a simple graph search.
#
# total_edge_data: information containing the DICs for the graph edge weights
# solution_object: The solution being reduced
# threshold_value: User specified threshold value
# keeper_list: Speicies that should always be kept
# done: Determines wether or not the reduction is complete
# target_species: The target species for the search in the array
#
# Returns an array of species that should be excluded from the original model at this threshold level
#################

def trim_pfa(total_edge_data, solution_object, threshold_value, keeper_list, done, target_species):

    start_time = tm.time()

    #initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    graph = nx.DiGraph() #Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.


    safe = [] #A list of species that are to be retained for this threshold value

    #calculate edge weights based on list received from get_rate_data
    #initial condition
    for ic in total_edge_data.keys(): #For each initial condition
        for species in species_objects: #Make graph
            graph.add_node(species.name)
        #timestep
        for tstep in total_edge_data[ic].keys(): #Set edge values for the graph
            number = total_edge_data[ic][tstep]
            #each species
            for edge in number:
                try:
                    edge_name = edge.split('_', 1)
                    species_a_name = edge_name[0]
                    species_b_name = edge_name[1]
                    #pfa weight between two species
                    weight = number[edge]
                    if graph.has_edge(species_a_name, species_b_name):
                        old_weight = graph[species_a_name][species_b_name]['weight']
                        if weight > old_weight and weight > threshold_value: #Only include the weight if it is greater than the threshold value.
                            graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
                        #elif weight > 1:
                        #    print("Error.  Edge weights should not be greater than one.")
                        #    exit()
                    elif weight > threshold_value:
                        graph.add_weighted_edges_from([(species_a_name, species_b_name, weight)])
                    #elif weight > 1:
                    #        print("Error.  Edge weights should not be greater than one.")
                    #        exit()
                except IndexError:
                    print(edge)
                    continue

            dic = graph_search(graph, target_species) #Search graph for max values to each species based on targets
            for sp in dic: #Add to max dictionary if it is new or greater than the value already there.
                if sp not in safe:
                    safe.append(sp)
            graph.clear() #Reset graph

    core_species = []
    species_objects = solution_object.species()

    #Take all species that are over the threshold value and add them to essentail species.
    essential_species = []
    for sp in species_objects:
        if sp.name in safe:
            if sp not in essential_species:
                essential_species.append(sp)
    done[0] = False

    #Add all species in essential species to core species
    for sp in essential_species:
        if sp not in core_species:
            core_species.append(sp.name)

    #Add all of the must keep species to core species
    retained_species = keeper_list #Specified by the user.  A list of species that also need to be kept.
    for sp in retained_species:
        if sp not in core_species:
            core_species.append(sp)

    exclusion_list = []

    #Exclude everythong not in core species.
    for species in solution_object.species():
        if species.name not in core_species: #If its not one of our species we must keep, add it to the list of species to be trimmed.
            exclusion_list.append(species.name)

    return exclusion_list

############
# This is the MAIN top level function for running PFA
#
# args: a list of user specified command line arguements
# solution_object: a Cantera object of the solution to be reduced
# error: To hold the error level of the simulation
# past: To hold error level of previous simulation
#
# Writes reduced Cantera file and returns reduced Catnera solution object
############

def run_pfa(args, solution_object,error,past):

	if len(args.target) == 0: #If the target species are not specified, puke and die.
		print("Please specify a target species.")
		exit()
	done = [] #Singleton to hold wether or not any more species can be cut from the simulation.
	done.append(False)
	threshold = .1 #Starting threshold value
	threshold_i = .1
	n = 1
	error = [10.0] #Singleton to hold the error value of the previously ran simulation.

	#Check to make sure that conditions exist
	if args.conditions_file:
		conditions_array = readin_conditions(str(args.conditions_file))
	elif not args.conditions_file:
		print("Conditions file not found")
		exit()

	sim_array = helper.setup_simulations(conditions_array,solution_object) #Turn conditions array into unran simulation objects for the original solution
	ignition_delay_detailed = helper.simulate(sim_array) #Run simulations and process results
	rate_edge_data = get_rates_pfa(sim_array, solution_object) #Get edge weight calculation data.

	print("Testing for starting threshold value")
	pfa_loop_control(solution_object, args, error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array) #Trim the solution at that threshold and find the error.
	while error[0] != 0: #While the error for trimming with that threshold value is greater than allowed.
		threshold = threshold / 10 #Reduce the starting threshold value and try again.
		threshold_i = threshold_i / 10
		n = n + 1
		pfa_loop_control(solution_object, args, error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array)
		if error[0] <= .02:
			error[0] = 0

	print("Starting with a threshold value of " + str(threshold))
	sol_new = solution_object
	past[0] = 0 #An integer representing the error introduced in the past simulation.
	done[0] = False

	while not done[0] and error[0] < args.error: #Run the simulation until nothing else can be cut.
		sol_new = pfa_loop_control( solution_object, args, error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array) #Trim at this threshold value and calculate error.
		if args.error > error[0]: #If a new max species cut without exceeding what is allowed is reached, save that threshold.
			max_t = threshold
		#if (past == error[0]): #If error wasn't increased, increase the threshold at a higher rate.
		#	threshold = threshold + (threshold_i * 4)
			past[0] = error[0]
		#if (threshold >= .01):
                #        threshold_i = .01
		threshold = threshold + threshold_i
		threshold = round(threshold, n)

	print("\nGreatest result: ")
	sol_new = pfa_loop_control( solution_object, args, error, max_t, done, rate_edge_data, ignition_delay_detailed, conditions_array)
	drgep_trimmed_file = soln2cti.write(sol_new) #Write the solution object with the greatest error that isn't over the allowed ammount.

	return sol_new[1]

#############
# This function handles the reduction, simulation, and comparision for a single threshold value
#
# solution_object: object being reduced
# args: arguments from the command line
# stored_error: past error
# threshold: current threshold value
# done: are we done reducing yet? Boolean
# rate_edge_data: the DICs for reduction
# ignition_delay_detailed: ignition delay of detailed model
# conditions_array: array holding information about initial conditions
#
# Returns the reduced solution object for this threshold and updates error value
#############

def pfa_loop_control(solution_object, args, stored_error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array):

    target_species = args.target

    #run detailed mechanism and retain initial conditions
    species_retained = []
    printout = ''
    print('Threshold     Species in Mech      Error')

    #run DRG and create new reduced solution
    pfa = trim_pfa(rate_edge_data, solution_object, threshold, args.keepers, done,target_species) #Find out what to cut from the model
    exclusion_list = pfa
    new_solution_objects = trim(solution_object, exclusion_list, args.data_file) #Cut the exclusion list from the model.
    species_retained.append(len(new_solution_objects[1].species()))

    #simulated reduced solution
    new_sim = helper.setup_simulations(conditions_array,new_solution_objects[1]) #Create simulation objects for reduced model for all conditions
    ignition_delay_reduced = helper.simulate(new_sim) #Run simulations and process results

    if (ignition_delay_detailed.all() == 0): #Ensure that ignition occured
        print("Original model did not ignite.  Check initial conditions.")
        exit()

    #Calculate error
    error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100 #Calculate error
    printout += str(threshold) + '                 ' + str(len(new_solution_objects[1].species())) + '              '+  str(round(np.max(error), 2)) +'%' + '\n'
    print(printout)
    stored_error[0] = round(np.max(error), 2)

    #Return new model.
    new_solution_objects = new_solution_objects[1]
    return new_solution_objects

############
# This function calculates values to be used in the calculation of Direct Interaction Coefficients
#
# sim_array: Array of simulated simulation objects
# solution_object: Cantera object of the solution being reduced
#
# Returns:
#   total_edge_data: a dictionary with keys of initial conditions and values of dictionarys that hold information for caculating DICs at each timestep.
#       *the subdictionaries have the timestep as their keys and their values hold an array of DICs
############

def get_rates_pfa(sim_array, solution_object):

    old_solution = solution_object
    #iterate through all initial conditions
    total_edge_data = {}
    for ic in sim_array:
        ic_edge_data = {}
        for tstep in ic.sample_points: #iterate through all timesteps
            temp = tstep[0] #Set up variables
            pressure = tstep[1]
            mass_fractions = np.array(tstep[2])

            #Set up solution at current timestep
            new_solution = old_solution
            new_solution.TPY = temp, pressure, mass_fractions
            new_reaction_production_rates = new_solution.net_rates_of_progress
            new_species_prod_rates=new_solution.net_production_rates

            DIC = {}

            single = get_PA(new_solution,new_reaction_production_rates) #Get PA and CA
            PA = single[0]
            CA = single[1]

            double = get_PAB(new_solution,new_reaction_production_rates) #Get PAB and CAB
            PAB = double[0]
            CAB = double[1]

            r1 = get_rAB_1(new_solution,PA,CA,PAB,CAB)
            rAB_p1 = r1[0]
            rAB_c1 = r1[1]

            r2 = get_rAB_2(new_solution,rAB_p1,rAB_c1)
            rAB_p2 = r2[0]
            rAB_c2 = r2[1]

            s_names = new_solution.species_names
            for species_a in s_names:
                for species_b in s_names:
                    if (species_a != species_b):
                        full_name = species_a + "_" + species_b
                        add = rAB_p1[full_name] + rAB_c1[full_name] + rAB_p2[full_name] + rAB_c2[full_name]
                        DIC[full_name] = add

            ic_edge_data[temp] = DIC
        total_edge_data[ic] = ic_edge_data
    return total_edge_data

######################
# Function: get_PA
# Description: Gets the PA (and CA) values of all species in a given solution.
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates.
# Input: A solution class from the Cantera library and an array of production rates.
# Output: Two dictionaries both keyed by species name.  One for PA, one for CA.
# Returns: PA and CA dictionaries.
#######################

def get_PA(new_solution, new_reaction_production_rates):
	PA = {} #Dictionary that will hold the PA values for each species.
	CA = {} #Dictionary that will hold the CA values for each species.

	s_names = new_solution.species_names
	for species in s_names: #For species in the solutuion
		PA[species] = 0
		CA[species] = 0

		for i, reac in enumerate(new_solution.reactions()): #For all reactions
			reac_prod_rate = float(new_reaction_production_rates[i]) #Set up values
			all_species = reac.reactants
			all_species.update(reac.products)

			if reac_prod_rate != 0:
				if species in all_species: #If the species is in a productive reaction, sum up omega times v in the appropriate dictionary.
					add = float(reac_prod_rate * all_species[species])
					if add > 0:
						PA[species] += abs(add)
					else:
						CA[species] += abs(add)

	return PA,CA


######################
# Function: get_PAB
# Description: Gets the PAB (and CAB) values of all species in a given solution.
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates.
# Input: A solution class from the Cantera library and an array of production rates.
# Output: Two dictionaries both keyed by species name.  One for PAB, one for CAB.
# Returns: PAB and CAB dictionaries.
#######################

def get_PAB(new_solution, new_reaction_production_rates):
	PAB = {} #Set up dictionaries
	CAB = {}

	s_names = new_solution.species_names
	for species_a in s_names: #For every pair of species A and B in the solution
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b
				PAB[full_name] = 0
				CAB[full_name] = 0

				for i, reac in enumerate(new_solution.reactions()): #For all reactions
					reac_prod_rate = float(new_reaction_production_rates[i]) #Set up values
					all_species = reac.reactants
					all_species.update(reac.products)

					if reac_prod_rate != 0: #If both species exsist in the reaction, add the calculated value to the correct dictionary.
						if species_a in all_species:
							if species_b in all_species:
								add = float(reac_prod_rate * all_species[species_a])
								if add > 0:
									PAB[full_name] += abs(add)
								else:
									CAB[full_name] += abs(add)

	return PAB, CAB

######################
# Function: get_rAB_1
# Description: Gets the rAB_p1 (and rAB_c1) values of all species in a given solution.
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates and 4 dictionaries.
# Input: A solution class from the Cantera library and all of the PA,CA,PAB, and CAB values.
# Output: Two dictionaries both keyed by species name.  One for rAB_p1, one for rAB_c1.
# Returns: rAB_p1 and rAB_c1 dictionaries.
#######################

def get_rAB_1(new_solution,PA,CA,PAB,CAB):
	rAB_p1 = {} #Set up dictionaries
	rAB_c1 = {}

	s_names = new_solution.species_names
	for species_a in s_names: #For all pairs of species
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b #Set up
				rAB_p1[full_name] = 0
				rAB_c1[full_name] = 0

				top_p = PAB[full_name] #Get top
				top_c = CAB[full_name]

				if (PA[species_a] > CA[species_a]): #Get bot
					bot = PA[species_a]
				else:
					bot = CA[species_a]

				if (bot != 0): #Calculate
					rAB_p1[full_name] = top_p/bot
					rAB_c1[full_name] = top_c/bot

	return rAB_p1, rAB_c1


######################
# Function: get_rAB_2
# Description: Gets the rAB_p2 (and rAB_c2) values of all species in a given solution.
# Parameters: Must be a valid model stored in the cantera format with corresponding production rates and 2 dictionaries.
# Input: A solution class from the Cantera library and all of the rAB_p1 and rAB_c1 values.
# Output: Two dictionaries both keyed by species name.  One for rAB_p2, one for rAB_c2.
# Returns: rAB_p2 and rAB_c2 dictionaries.
#######################

def get_rAB_2(new_solution,rAB_p1,rAB_c1):
	rAB_p2 = {} #Set up dictionaries
	rAB_c2 = {}

	s_names = new_solution.species_names
	for species_a in s_names: #For all pairs of species
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b #Set up
				rAB_p2[full_name] = 0
				rAB_c2[full_name] = 0

				for species_m in s_names: #Look through all possible middle step species
					if (species_m != species_a and species_m != species_b):
						am_name = species_a + "_" + species_m
						mb_name = species_m + "_" + species_b

						add_p = rAB_p1[am_name] * rAB_p1[mb_name] #Get what to add for species_m
						add_c = rAB_c1[am_name] * rAB_c1[mb_name]

						rAB_p2[full_name] += add_p #Add that value
						rAB_c2[full_name] += add_c

	return rAB_p2,rAB_c2
