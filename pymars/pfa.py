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

def trim_pfa(total_edge_data, solution_object, threshold_value, keeper_list, done, target_species, model_file):

    """
    This function determines what species should be excluded from the reduced model based on their DICs compared to the threshold value and a simple graph search.

    Parameters
    ----------

    total_edge_data: information containing the DICs for the graph edge weights
    solution_object: The solution being reduced
    threshold_value: User specified threshold value
    keeper_list: Speicies that should always be kept
    done: Determines wether or not the reduction is complete
    target_species: The target species for the search in the array
    model_file: File holding the model being reduced

    Returns 
    -------

    Returns an array of species that should be excluded from the original model at this threshold level
   
    """          

    start_time = tm.time()

    # Initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    
    # Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.
    graph = nx.DiGraph()


    safe = [] # A list of species that are to be retained for this threshold value

    # Calculate edge weights based on list received from get_rate_data
    # Initial condition
    for ic in total_edge_data.keys(): # For each initial condition
        for species in species_objects: # Make graph
            graph.add_node(species.name)
        # Timestep
        for tstep in total_edge_data[ic].keys(): # Set edge values for the graph
            number = total_edge_data[ic][tstep]
            # Each species
            for edge in number:
                try:
                    edge_name = edge.split('_', 1)
                    species_a_name = edge_name[0]
                    species_b_name = edge_name[1]
                    # PFA weight between two species
                    weight = number[edge]
                    if graph.has_edge(species_a_name, species_b_name):
                        old_weight = graph[species_a_name][species_b_name]['weight']
                        if weight > old_weight and weight > threshold_value: # Only include the weight if it is greater than the threshold value.
                            graph.add_weighted_edges_from(
                                [(species_a_name, species_b_name, weight)])
                        #elif weight > 1:
                        #    print("Error.  Edge weights should not be greater than one.")
                        #    exit()
                    elif weight > threshold_value:
                        graph.add_weighted_edges_from(
                            [(species_a_name, species_b_name, weight)])
                    #elif weight > 1:
                    #        print("Error.  Edge weights should not be greater than one.")
                    #        exit()
                except IndexError:
                    print(edge)
                    continue

            dic = graph_search(graph, target_species) # Search graph for max values to each species based on targets
            for sp in dic: # Add to max dictionary if it is new or greater than the value already there.
                if sp not in safe:
                    safe.append(sp)
            graph.clear() # Reset graph

    core_species = []
    species_objects = solution_object.species()

    # Take all species that are over the threshold value and add them to essentail species.
    essential_species = []
    for sp in species_objects:
        if sp.name in safe:
            if sp not in essential_species:
                essential_species.append(sp)
    done[0] = False

    # Add all species in essential species to core species
    for sp in essential_species:
        if sp not in core_species:
            core_species.append(sp.name)

    # Add all of the must keep species to core species
    retained_species = keeper_list # Specified by the user.  A list of species that also need to be kept.
    for sp in retained_species:
        if sp not in core_species:
            core_species.append(sp)

    exclusion_list = []

    # Exclude everything not in core species.
    for species in solution_object.species():
        # If its not one of our species we must keep, add it to the list of species to be trimmed.
        if species.name not in core_species:
            exclusion_list.append(species.name)

    return exclusion_list


def run_pfa(solution_object, conditions_file, error_limit, target_species, retained_species, model_file, final_error):

	""""
	This is the MAIN top level function for running PFA

	Parameters
	----------

	solution_object: a Cantera object of the solution to be reduced
	conditions_file: The file holding the initial conditions to simulate
	error_limit: The highest allowed error percentage
	target_species: The target species for reduction
	retained_species: A list of species to be retained even if they should be cut by the algorithm
	model_file: The path to the file where the solution object was generated from
	final_error: To hold the error level of the simulation

	Returns
	-------

	Writes reduced Cantera file and returns reduced Catnera solution object
	
	"""    

	if len(target_species) == 0: # If the target species are not specified, puke and die.
		print("Please specify a target species.")
		exit()
	done = [] # Singleton to hold wether or not any more species can be cut from the simulation.
	done.append(False)
	threshold = .1 # Starting threshold value
	threshold_i = .1
	n = 1
	error = [10.0] # Singleton to hold the error value of the previously ran simulation.

	# Check to make sure that conditions exist
	if conditions_file:
		conditions_array = readin_conditions(str(conditions_file))
	elif not conditions_file:
		print("Conditions file not found")
		exit()

	# Turn conditions array into unran simulation objects for the original solution
	sim_array = helper.setup_simulations(conditions_array,solution_object)
	ignition_delay_detailed = helper.simulate(sim_array) #Run simulations and process results
	rate_edge_data = get_rates_pfa(sim_array, solution_object) #Get edge weight calculation data.

	print("Testing for starting threshold value")
	
	# Trim the solution at that threshold and find the error.
	pfa_loop_control(
		solution_object, target_species, retained_species, model_file, error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array)
	while error[0] != 0 and threshold_i > .001: # While the error for trimming with that threshold value is greater than allowed.
		threshold = threshold / 10 # Reduce the starting threshold value and try again.
		threshold_i = threshold_i / 10
		n = n + 1
		pfa_loop_control(
			solution_object, target_species, retained_species, model_file, error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array)
		if error[0] <= .02:
			error[0] = 0

	print("Starting with a threshold value of " + str(threshold))
	
	sol_new = solution_object
	final_error[0] = 0 # An integer representing the error introduced in the final simulation.
	done[0] = False

	while not done[0] and error[0] < error_limit: # Run the simulation until nothing else can be cut.
		# Trim at this threshold value and calculate error.
		sol_new = pfa_loop_control(
			solution_object, target_species, retained_species, model_file, error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array)
		if error_limit >= error[0]: # If a new max species cut without exceeding what is allowed is reached, save that threshold.
			max_t = threshold
		#if (final_error == error[0]): #If error wasn't increased, increase the threshold at a higher rate.
		#	threshold = threshold + (threshold_i * 4)
			final_error[0] = error[0]
		#if (threshold >= .01):
                #        threshold_i = .01
		threshold = threshold + threshold_i
		threshold = round(threshold, n)

	print("\nGreatest result: ")
	sol_new = pfa_loop_control(
		solution_object, target_species, retained_species, model_file, error, max_t, done, rate_edge_data, ignition_delay_detailed, conditions_array)
	
	return sol_new


def pfa_loop_control(solution_object, target_species, retained_species, model_file, stored_error, threshold, done, rate_edge_data, ignition_delay_detailed, conditions_array):

    """
    This function handles the reduction, simulation, and comparision for a single threshold value

    Parameters
    ----------

    solution_object: object being reduced # target_species:
    target_species: The target species for reduction
    retained_species: A list of species to be retained even if they should be cut by the algorithm
    model_file: The path to the file where the solution object was generated from
    stored_error: Error from the previous reduction attempt
    threshold: current threshold value
    done: are we done reducing yet? Boolean
    rate_edge_data: the DICs for reduction
    ignition_delay_detailed: ignition delay of detailed model
    conditions_array: array holding information about initial conditions

    Returns
    -------

    Returns the reduced solution object for this threshold and updates error value
    
    """

    # Run detailed mechanism and retain initial conditions
    species_retained = []
    printout = ''
    print('Threshold     Species in Mech      Error')

    # Run DRG and create new reduced solution
    exclusion_list = trim_pfa(
        rate_edge_data, solution_object, threshold, retained_species, done,target_species,model_file) # Find out what to cut from the model
    new_solution_objects = trim(solution_object, exclusion_list, model_file) # Cut the exclusion list from the model.
    species_retained.append(len(new_solution_objects[1].species()))

    # Simulated reduced solution
    new_sim = helper.setup_simulations(conditions_array,new_solution_objects[1]) # Create simulation objects for reduced model for all conditions
    ignition_delay_reduced = helper.simulate(new_sim) # Run simulations and process results

    if (ignition_delay_detailed.all() == 0): # Ensure that ignition occured
        print("Original model did not ignite.  Check initial conditions.")
        exit()

    # Calculate error
    error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100 # Calculate error
    printout += str(threshold) + '                 ' + str(len(new_solution_objects[1].species())) + '              '+  str(round(np.max(error), 2)) +'%' + '\n'
    print(printout)
    stored_error[0] = round(np.max(error), 2)

    # Return new model
    new_solution_objects = new_solution_objects[1]
    return new_solution_objects


def get_rates_pfa(sim_array, solution_object):

    """     
    This function calculates values to be used in the calculation of Direct Interaction Coefficients

    Parameters
    ----------

    sim_array: Array of simulated simulation objects
    solution_object: Cantera object of the solution being reduced

    Returns 
    -------

      total_edge_data: a dictionary with keys of initial conditions and values of dictionarys that hold information for caculating DICs at each timestep.
          *the subdictionaries have the timestep as their keys and their values hold an array of DICs
    
    """
    
    old_solution = solution_object
    
    # Iterate through all initial conditions
    total_edge_data = {}
    for ic in sim_array:
        ic_edge_data = {}
        for tstep in ic.sample_points: # Iterate through all timesteps
            temp = tstep[0] # Set up variables
            pressure = tstep[1]
            mass_fractions = np.array(tstep[2])

            # Set up solution at current timestep
            new_solution = old_solution
            new_solution.TPY = temp, pressure, mass_fractions
            new_reaction_production_rates = new_solution.net_rates_of_progress
            new_species_prod_rates=new_solution.net_production_rates

            DIC = {}

            single = get_PA(new_solution,new_reaction_production_rates) # Get PA and CA
            PA = single[0]
            CA = single[1]

            double = get_PAB(new_solution,new_reaction_production_rates) # Get PAB and CAB
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


def get_PA(new_solution, new_reaction_production_rates):

	"""
	Gets the PA (and CA) values of all species in a given solution.
	
	Parameters
	----------

	new_solution: The object representing the cantera model
	new_reaction_production_rates: the production rates associated with the model

	Returns
	-------

	PA and CA dictionaries
	
	"""

	PA = {} # Dictionary that will hold the PA values for each species.
	CA = {} # Dictionary that will hold the CA values for each species.

	# Initalize all species
	s_names = new_solution.species_names
	for species in s_names: 
		PA[species] = 0
		CA[species] = 0

	for i, reac in enumerate(new_solution.reactions()): # For all reactions
		reac_prod_rate = float(new_reaction_production_rates[i]) # Set up values
				
		if reac_prod_rate != 0:
			if reac_prod_rate > 0: # For forward reactions

				# Add all products to PA
				for species in reac.products:
					add = float(reac_prod_rate * reac.products[species])
					PA[species] += abs(add)
				# Add all reactants to CA
				for species in reac.reactants:
					add = float(reac_prod_rate * reac.reactants[species])
					CA[species] += abs(add)
					
			if reac_prod_rate < 0: # For forward reactions
			
				# Add all products to CA	
				for species in reac.products:
					add = float(reac_prod_rate * reac.products[species])
					CA[species] += abs(add)

				# Add all reactants to PA
				for species in reac.reactants:
					add = float(reac_prod_rate * reac.reactants[species])
					PA[species] += abs(add)

	return PA,CA


def get_PAB(new_solution, new_reaction_production_rates):
	
	"""
	Gets the PAB (and CAB) values of all species in a given solution.
	
	Parameters
	----------

	new_solution: The object representing the cantera model
	new_reaction_production_rates: the production rates associated with the model

	Returns
	-------
	
	PAB and CAB dictionaries.
	
	"""
	
	PAB = {} # Set up dictionaries
	CAB = {}

	s_names = new_solution.species_names
	for species_a in s_names: # For every pair of species A and B in the solution
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b
				PAB[full_name] = 0
				CAB[full_name] = 0

				for i, reac in enumerate(new_solution.reactions()): # For all reactions
					reac_prod_rate = float(new_reaction_production_rates[i]) # Set up values
					all_species = reac.products
					all_species.update(reac.reactants)

					# If both species exsist in the reaction, add the calculated value to the correct dictionary.
					if reac_prod_rate != 0:
						if species_a in all_species:
							if species_b in all_species:
								
								# For forward reactions
								if reac_prod_rate > 0:
						
									# Add products to PAB
									if species_a in reac.products:
										add = float(reac_prod_rate * reac.products[species_a])
										PAB[full_name] += abs(add)

									# Add reactants to CAB
									if species_a in reac.reactants:
										add = float(reac_prod_rate * reac.reactants[species_a])
										CAB[full_name] += abs(add)
				
								# For backward reactions	
								if reac_prod_rate < 0:
									
									# Add products to CAB
									if species_a in reac.products:
										add = float(reac_prod_rate * reac.products[species_a])
										CAB[full_name] += abs(add)
									
									# Add reactants to PAB
									if species_a in reac.reactants:
										add = float(reac_prod_rate * reac.reactants[species_a])
										PAB[full_name] += abs(add)
	return PAB, CAB


def get_rAB_1(new_solution,PA,CA,PAB,CAB):
	
	"""
	Gets the rAB_p1 (and rAB_c1) values of all species in a given solution.
	
	Parameters
	----------
	
	new_solution: The object representing the cantera model
	PA: A dictionary containing the PA values for the reduction
	CA: A dictionary containing the CA values for the reduction
	PAB: A dictionary containing the PAB values for the reduction
	CAB: A dictionary containing the CAB values for the reduction
	
	Returns
	-------
	
	rAB_p1 and rAB_c1 dictionaries.

	"""
	
	rAB_p1 = {} # Set up dictionaries
	rAB_c1 = {}

	s_names = new_solution.species_names
	for species_a in s_names: # For all pairs of species
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b # Set up
				rAB_p1[full_name] = 0
				rAB_c1[full_name] = 0

				top_p = PAB[full_name] # Get numerator
				top_c = CAB[full_name]

				if (PA[species_a] > CA[species_a]): # Get denomonator
					bot = PA[species_a]
				else:
					bot = CA[species_a]

				if (bot != 0): # Calculate
					rAB_p1[full_name] = top_p/bot
					rAB_c1[full_name] = top_c/bot

	return rAB_p1, rAB_c1


def get_rAB_2(new_solution,rAB_p1,rAB_c1):
	
	"""
	Gets the rAB_p2 (and rAB_c2) values of all species in a given solution.
	
	Parameters
	----------
	
	new_solution: The object representing the cantera model
	rAB_p1: A dictionary containing the rAB_p1 values for the reduction
	rAB_c1: A dictionary containing the rAB_c1 values for the reduction
	
	Returns
	-------

	rAB_p2 and rAB_c2 dictionaries.
	
	"""

	rAB_p2 = {} # Set up dictionaries
	rAB_c2 = {}

	s_names = new_solution.species_names
	for species_a in s_names: # For all pairs of species
		for species_b in s_names:
			if species_a != species_b:
				full_name = species_a + "_" + species_b # Set up
				rAB_p2[full_name] = 0
				rAB_c2[full_name] = 0

				# Look through all possible middle step species
				for species_m in s_names:
					if (species_m != species_a and species_m != species_b):
						am_name = species_a + "_" + species_m
						mb_name = species_m + "_" + species_b

						# Get what to add for species_m
						add_p = rAB_p1[am_name] * rAB_p1[mb_name]
						add_c = rAB_c1[am_name] * rAB_c1[mb_name]

						# Add that value
						rAB_p2[full_name] += add_p
						rAB_c2[full_name] += add_c

	return rAB_p2,rAB_c2
