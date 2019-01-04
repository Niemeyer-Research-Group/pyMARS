import networkx as nx
import numpy as np
import h5py
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


def make_dic_drgep(solution_object, total_edge_data, target_species):

    """     
    Makes the dictionary of overall interaction coefficients for DRGEP by building a graph and searching it as explained in DRGEP

    Parameters
    ----------

    solution_object: the Cantera soluton object that is being reduced
    target_species: an array of target species for the reduction as specified by the user
    total_edge_data: a dictionary with keys of initial conditions and values of dictionarys that hold information for caculating DICs at each timestep.
      *the subdictionaries have the timestep as their keys and their values hold an array of numberator and denominator information for calculating DICs.

    Returns
    -------

    A dictionary with keys of species and values of that species OIC
    
    """

    # Initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    
    # Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.
    graph = nx.DiGraph()

    max_dic = {} # Dictionary holding the maximum values for the iteration

    # Calculate edge weights based on list received from get_rate_data and use them to create a graph
    for ic in total_edge_data.keys(): # For each initial condition
        for species in species_objects: # Make graph
            graph.add_node(species.name)
        # Timestep
        for tstep in total_edge_data[ic].keys(): # Make a graph at each timestep
            # DRGEP calculations of direct interaction coeffients are done in the total_edge_data function.
            numerator = total_edge_data[ic][tstep][1]
            denominator = total_edge_data[ic][tstep][0]
            for edge in numerator: # For each edge, determine its weight amnd add it to the graph
                try:
                    edge_name = edge.split('_', 1)
                    species_a_name = edge_name[0]
                    species_b_name = edge_name[1]
                    # DRGEP weight between two species
                    if denominator[species_a_name] != 0:
                        weight = abs(
                            float(numerator[edge])/float(denominator[species_a_name]))
                        if graph.has_edge(species_a_name, species_b_name):
                            old_weight = graph[species_a_name][species_b_name]['weight']
                            if weight > old_weight and weight <= 1:
                                graph.add_weighted_edges_from(
                                    [(species_a_name, species_b_name, weight)])
                            elif weight > 1:
                                print(
                                    "Error.  Edge weights should not be greater than one.")
                                exit()
                        elif weight <= 1:
                            graph.add_weighted_edges_from(
                                [(species_a_name, species_b_name, weight)])
                        elif weight > 1:
                            print(
                                "Error.  Edge weights should not be greater than one.")
                            exit()
                except IndexError:
                    print(edge)
                    continue

            # Search the graph for overall interaction coefficents and add them to max_dic if they belong
            dic = graph_search_drgep(graph, target_species) # Search graph for max values to each species based on targets
            for sp in dic: # Add to max dictionary if it is new or greater than the value already there.
                if sp not in max_dic:
                    max_dic[sp] = dic[sp]
                elif dic[sp] > max_dic[sp]:
                    max_dic[sp] = dic[sp]
            graph.clear() # Reset graph
    return max_dic


def trim_drgep(max_dic, solution_object, threshold_value, retained_species, done):
    
    """
    Determines what species should be excluded from the reduced model based on their OICs compared to the threshold value.

    Parameters
    ----------
    max_dic: Dictionary of OICs for all species
    solution_object: The solution being reduced
    threshold_value: User specified threshold value
    keeper_list: Speicies that should always be kept
    done: Determines wether or not the reduction is complete

    Returns
    -------
    An array of species that should be excluded from the original model at this threshold level
    
    """
    
    core_species = []
    species_objects = solution_object.species()

    for sp in retained_species:
        core_species.append(sp)
    
    # Take all species that are over the threshold value and add them to essential species.
    essential_species = []
    for sp in species_objects:
        if sp.name in max_dic:
            if max_dic[sp.name] > threshold_value and sp not in essential_species:
                essential_species.append(sp)
    done[0] = True
    for sp in species_objects: # If any more can be taken out, we are not done yet.
        if sp.name in max_dic:
            if max_dic[sp.name] > threshold_value:
                done[0] = False

    for sp in essential_species:
        if sp not in core_species:
            core_species.append(sp.name)

    exclusion_list = []

    for species in solution_object.species():
        # If its not one of our species we must keep, add it to the list of species to be trimmed.
        if species.name not in core_species: 
            exclusion_list.append(species.name)

    return exclusion_list


def run_drgep(solution_object, conditions_file, error_limit, target_species, retained_species, model_file, final_error):
    
	""""
	This is the MAIN top level function for running DRGEP

	Parameters
	----------

	solution_object: A Cantera object of the solution to be reduced
	conditions_file: A file holding the initial conditions for the simulation
	error_limit: The maximum error allowed for the reduced model
	target_species: An array of the target species for reduction
	retained_species: An array of species that should not be removed from the model
	model_file: The path to the file holding the original model
	final_error: A singleton holding the final error percentage

	Returns
	-------

	Writes reduced Cantera file and returns reduced Catnera solution object
    
	"""

	# Set up variables
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

	sim_array = helper.setup_simulations(conditions_array,solution_object) # Turn conditions array into unran simulation objects for the original solution
	ignition_delay_detailed = helper.simulate(sim_array) # Run simulations and process results

	rate_edge_data = get_rates(sim_array,solution_object) # Get edge weight calculation data.
	max_dic = make_dic_drgep(solution_object, rate_edge_data, target_species) # Make a dictionary of overall interaction coefficients.

	# Trim the solution at that treshold and find the error.
	print("Testing for starting threshold value")
	drgep_loop_control(
		solution_object, target_species, retained_species, model_file, error, threshold, done, max_dic, ignition_delay_detailed, conditions_array)
	
	while error[0] != 0 and threshold_i > .001: # While the error for trimming with that threshold value is greater than allowed.
		threshold = threshold / 10 # Reduce the starting threshold value and try again.
		threshold_i = threshold_i / 10
		n = n + 1
		drgep_loop_control(
			solution_object, target_species, retained_species, model_file, error, threshold, done, max_dic, ignition_delay_detailed, conditions_array)
		if error[0] <= .02:
			error[0] = 0

	print("Starting with a threshold value of " + str(threshold))
	sol_new = solution_object
	final_error[0] = 0 # An integer representing the error introduced in the past simulation.
	done[0] = False

	while not done[0] and error[0] < error_limit: # Run the simulation until nothing else can be cut.
		# Trim at this threshold value and calculate error.
		sol_new = drgep_loop_control(
			solution_object, target_species, retained_species, model_file, error, threshold, done, max_dic, ignition_delay_detailed, conditions_array)
		if error_limit > error[0]: # If a new max species cut without exceeding what is allowed is reached, save that threshold.
			max_t = threshold
            #if (final_error[0] == error[0]): #If error wasn't increased, increase the threshold at a higher rate.
		        #threshold = threshold + (threshold_i * 4)
			final_error[0] = error[0]
		    #if (threshold >= .01):
                #threshold_i = .01
			threshold = threshold + threshold_i
			threshold = round(threshold, n)

	print("\nGreatest result: ")
	sol_new = drgep_loop_control(
		solution_object, target_species, retained_species, model_file, error, max_t, done, max_dic, ignition_delay_detailed, conditions_array)
	
	return sol_new


def drgep_loop_control(solution_object, target_species, retained_species, model_file, stored_error, threshold, done, max_dic, ignition_delay_detailed, conditions_array):

    """      
    This function handles the reduction, simulation, and comparision for a single threshold value.

    Parameters
    ----------

    solution_object: object being reduced
    target_species: An array of the target species for reduction
    retained_species: An array of species that should not be removed from the model
    model_file: The path to the file holding the original model
    stored_error: past error
    threshold: current threshold value
    done: are we done reducing yet? Boolean
    max_dic: OIC dictionary for DRGEP
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

    # Run DRGEP and create new reduced solution
    exclusion_list = trim_drgep(max_dic, solution_object, threshold, retained_species, done) # Find out what to cut from the model
    new_solution_objects = trim(solution_object, exclusion_list, model_file) # Cut the exclusion list from the model.
    species_retained.append(len(new_solution_objects[1].species()))

    # Simulated reduced solution
    new_sim = helper.setup_simulations(conditions_array,new_solution_objects[1]) # Create simulation objects for reduced model for all conditions
    ignition_delay_reduced = helper.simulate(new_sim) # Run simulations and process results

    if ignition_delay_detailed.all() == 0: # Ensure that ignition occured
        print("Original model did not ignite.  Check initial conditions.")
        exit()

    # Calculate and print error.
    error = (abs(ignition_delay_reduced-ignition_delay_detailed)/ignition_delay_detailed)*100 # Calculate error
    printout += str(threshold) + '                 ' + str(len(new_solution_objects[1].species())) + '              '+  str(round(np.max(error), 2)) +'%' + '\n'
    print(printout)
    stored_error[0] = round(np.max(error), 2)

    # Return new model
    new_solution_objects = new_solution_objects[1]
    return new_solution_objects


def get_rates(sim_array, solution_object):

    """     
    This function calculates values to be used in the calculation of Direct Interaction Coefficients

    Parameters
    ----------
    
    sim_array: Array of simulated simulation objects
    solution_object: Cantera object of the solution being reduced

    Returns
    -------

    Returns a dictionary for calculating interaction coefficients as described in the make_dic_drgep function comments

    """

    # Initialize solution
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
            new_species_prod_rates = new_solution.net_production_rates

            denom = {}
            numerator = {}
            for spc in new_solution.species():
                denom[spc.name] = []
                denom[spc.name].append(0)
                denom[spc.name].append(0)

            # Calculate direct interaction coefficients as specified by the DRGEP method
            for i, reac in enumerate(new_solution.reactions()): # For all reactions
                reac_prod_rate = float(new_reaction_production_rates[i])
                reactants = reac.reactants
                products = reac.products
                all_species = reac.reactants
                all_species.update(reac.products)

                if reac_prod_rate != 0:
                    if reac_prod_rate > 0:
                        for species in products: # Add to denominator for all of the species in products
                            denom[species][1] += abs(float(reac_prod_rate * products[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator: # Add to numerator for all species pairs
                                        numerator[partial_name] += float(reac_prod_rate * products[species])
                                    else:
                                        numerator[partial_name] = float(reac_prod_rate * products[species])

                        for species in reactants: # For all reactants subtract instead of add
                            denom[species][0] += abs(float(reac_prod_rate * reactants[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator:
                                        numerator[partial_name] += float(-reac_prod_rate * reactants[species])
                                    else:
                                        numerator[partial_name] = float(-reac_prod_rate * reactants[species])

                    if reac_prod_rate < 0: # Same as above but for negative.
                        for species in products:
                            denom[species][0] += abs(float(reac_prod_rate * products[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator:
                                         numerator[partial_name] += float(reac_prod_rate * products[species])
                                    else:
                                         numerator[partial_name] = float(reac_prod_rate * products[species])

                        for species in reactants:
                            denom[species][1] += abs(float(reac_prod_rate * reactants[species]))
                            for species_b in all_species:
                                if species_b != species:
                                    partial_name = species + '_' + species_b
                                    if partial_name in numerator:
                                         numerator[partial_name] += float(-reac_prod_rate * reactants[species])
                                    else:
                                         numerator[partial_name] = float(-reac_prod_rate * reactants[species])

            for species in new_solution.species(): # Use greater value as denominator
                if abs(denom[species.name][0]) > abs(denom[species.name][1]):
                    denom[species.name] = abs(denom[species.name][0])
                else:
                    denom[species.name] = abs(denom[species.name][1])

            for name in numerator: # Use absolute value of numerator
                numerator[name] = abs(numerator[name])
            ic_edge_data[temp] = [denom, numerator] # Add to ic data
        total_edge_data[ic] = ic_edge_data # Add to information to be used to make the graph

    return total_edge_data


def graph_search_drgep(nx_graph, target_species):

    """
    This function searches the graph to generate a dictionary of the greatest paths to all species from one of the targets.
   
    Parameters
    ----------
    nx_graph : obj
        networkx graph object of solution
    target_species : list
        List of target species to search from

    Returns
    -------
    max_dic : dictionary
        Values of the greatest possible path to each species from one of the targets on this graph keyed by species name.
    
    """

    max_dic = {} # A dictionary holding the maximum path to each species.
    for target in target_species:
        dic = ss_dijkstra_path_length_modified(nx_graph, target) # Get dictionary of each values maximum path
        for sp in dic: # If the species are not in the max dictionary or the new value for that species is greater than the one in the max dictionary, add it to the max dictionary.
            if sp not in max_dic:
                max_dic[sp] = dic[sp]
            elif max_dic[sp] < dic[sp]:
                max_dic[sp] = dic[sp]
    return max_dic
