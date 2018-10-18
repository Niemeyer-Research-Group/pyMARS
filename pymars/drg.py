"""Module containing Directed Relation Graph (DRG) reduction method."""
from collections import Counter
import time
import os
import sys

import networkx
import numpy as np
import h5py
import cantera as ct

import soln2ck
import soln2cti
import helper
from simulation import Simulation
from create_trimmed_model import trim
from readin_initial_conditions import readin_conditions

os.environ['Cantera_Data'] = os.getcwd()

def trim_drg(total_edge_data, solution_object, threshold_value, keeper_list, done, target_species):
    """Determines which species to remove based on their DICs compared to the threshold value and a simple graph search.

    Parameters
    ----------
    total_edge_data :
        information for calculating the DICs for the graph edge weights
    solution_object :
        The solution being reduced
    threshold_value :
        User specified threshold value
    keeper_list :
        Species that should always be kept
    done :
        Determines wether or not the reduction is complete
    target_species :
        The target species for the search in the array

    Returns
    -------
    Array of species that should be reduced at this threshold level

    """

    start_time = time.time()

    # initalize solution and components
    solution = solution_object
    species_objects = solution.species()
    reaction_objects = solution.reactions()
    # Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.
    graph = networkx.DiGraph()

    safe = []  # A list of species that are to be retained for this threshold value

    # calculate edge weights based on list received from get_rate_data
    # initial condition
    for ic in total_edge_data.keys():  # For each initial condition
        for species in species_objects:  # Make graph
            graph.add_node(species.name)
        # timestep
        # Set edge values for the graph
        for tstep in total_edge_data[ic].keys():
            numerator = total_edge_data[ic][tstep][1]
            denominator = total_edge_data[ic][tstep][0]
            # each species
            for edge in numerator:
                try:
                    edge_name = edge.split('_', 1)
                    species_a_name = edge_name[0]
                    species_b_name = edge_name[1]
                    # dge weight between two species
                    if denominator[species_a_name] != 0:
                        weight = abs(
                            float(numerator[edge]) / float(denominator[species_a_name]))
                        if graph.has_edge(species_a_name, species_b_name):
                            old_weight = graph[species_a_name][species_b_name]['weight']
                            # Only include the weight if it is greater than the threshold value.
                            if weight > old_weight and weight <= 1 and weight > threshold_value:
                                graph.add_weighted_edges_from(
                                    [(species_a_name, species_b_name, weight)])
                            elif weight > 1:
                                print(
                                    "Error.  Edge weights should not be greater than one.")
                                exit()
                        elif threshold_value <= weight <= 1.0:
                            graph.add_weighted_edges_from(
                                [(species_a_name, species_b_name, weight)])
                        elif weight > 1:
                            print(
                                "Error.  Edge weights should not be greater than one.")
                            exit()
                except IndexError:
                    print(edge)
                    continue

            # Search graph for max values to each species based on targets
            dic = graph_search(graph, target_species)
            for sp in dic:  # Add to max dictionary if it is new or greater than the value already there.
                if sp not in safe:
                    safe.append(sp)
            graph.clear()  # Reset graph

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
    # Specified by the user.  A list of species that also need to be kept.
    retained_species = keeper_list
    for sp in retained_species:
        if sp not in core_species:
            core_species.append(sp)

    exclusion_list = []

    # Exclude everythong not in core species.
    for species in solution_object.species():
        # If its not one of our species we must keep, add it to the list of species to be trimmed.
        if species.name not in core_species:
            exclusion_list.append(species.name)

    return exclusion_list

def run_drg(solution_object, conditions, error_limit, target_species, retained_species):
    """Main function for running DRG reduction.

    Parameters
    ----------
    solution_object : ~cantera.Solution
        a Cantera object of the solution to be reduced.
    conditions : str
        Name of file with list of autoignition initial conditions.
    error_limit : float
        Maximum allowable error level for reduced model.
    target_species : list of str
        List of target species
    retained_species : list of str
        List of species to always be retained

    Returns
    -------

    Writes reduced Cantera file and returns reduced Cantera solution object

    """
    assert target_species, 'Need to specify at least one target species.'

    # Singleton to hold whether any more species can be cut from the simulation.
    done = []
    done.append(False)
    threshold = 0.1  # Starting threshold value
    threshold_increment = 0.1
    num_iterations = 1

    conditions_array = readin_conditions(conditions)

    # Turn conditions array into unrun simulation objects for the original solution
    sim_array = helper.setup_simulations(conditions_array, solution_object)
    # Run simulations and process results
    ignition_delay_detailed = helper.simulate(sim_array)
    # Get edge weight calculation data.
    rate_edge_data = get_rates_drg(sim_array, solution_object)

    print("Testing for starting threshold value")
    # Trim the solution at that threshold and find the error.
    drg_loop_control(
        solution_object, target_species, retained_species, threshold, done, rate_edge_data,
        ignition_delay_detailed, conditions_array
        )
    # While the error for trimming with that threshold value is greater than allowed.
    while error[0] != 0:
        # Reduce the starting threshold value and try again.
        threshold /= 10
        threshold_increment /= 10
        num_iterations += 1
        drg_loop_control(
            solution_object, target_species, retained_species, threshold, done, rate_edge_data,
            ignition_delay_detailed, conditions_array
            )
        if error[0] <= 0.02:
            error[0] = 0

    print("Starting with a threshold value of " + str(threshold))
    sol_new = solution_object

    # Run the simulation until nothing else can be cut.
    while not done[0] and error[0] < error:
        #Trim at this threshold value and calculate error.
        sol_new = drg_loop_control(
            solution_object, target_species, retained_species, threshold, done, rate_edge_data,
            ignition_delay_detailed, conditions_array
            )
        # If a new max species cut without exceeding what is allowed is reached, save that threshold.
        if error_limit > error[0]:
            max_t = threshold
        # if (past == error[0]): #If error wasn't increased, increase the threshold at a higher rate.
        #	threshold = threshold + (threshold_increment * 4)
        # if (threshold >= .01):
        #        threshold_increment = .01
        threshold += threshold_increment
        threshold = round(threshold, num_iterations)

    print("Greatest result: ")
    sol_new = drg_loop_control(
        solution_object, target_species, retained_species, max_t, done, rate_edge_data,
        ignition_delay_detailed, conditions_array
        )
    # Write the solution object with the greatest error that isn't over the allowed ammount.
    drgep_trimmed_file = soln2cti.write(sol_new)

    return sol_new[1]

def drg_loop_control(solution_object, target_species, retained_species, threshold, done, rate_edge_data,
                     ignition_delay_detailed, conditions_array
                     ):
    """Handles the reduction, simulation, and comparision for a single threshold value.

    Parameters
    ----------
    solution_object:
        object being reduced
    args:
        arguments from the command line
    threshold : float
        current threshold value
    done : bool
        are we done reducing yet?
    rate_edge_data :
        information for calculating the DICs for reduction
    ignition_delay_detailed :
        ignition delay of detailed model
    conditions_array :
        array holding information about initial conditions

    Returns
    -------
    Reduced solution object for this threshold and updates error value

    """

    species_retained = []
    printout = ''
    print('Threshold     Species in Mech      Error')

    # run DRG and create new reduced solution
    # Find out what to cut from the model
    drgep = trim_drg(rate_edge_data, solution_object, threshold, retained_species, done, target_species)
    exclusion_list = drgep
    # Cut the exclusion list from the model.
    new_solution_objects = trim(solution_object, exclusion_list)
    species_retained.append(len(new_solution_objects[1].species()))

    # simulated reduced solution
    # Create simulation objects for reduced model for all conditions
    new_sim = helper.setup_simulations(conditions_array, new_solution_objects[1])
    ignition_delay_reduced = helper.simulate(
        new_sim)  # Run simulations and process results

    if ignition_delay_detailed.all() == 0:  # Ensure that ignition occured
        print("Original model did not ignite.  Check initial conditions.")
        exit()

    # Calculate error
    error = (abs(ignition_delay_reduced -ignition_delay_detailed) /ignition_delay_detailed) *100  # Calculate error
    printout += str(threshold) + '                 ' + str(len(new_solution_objects[1].species())) + '              ' + str(round(np.max(error), 2)) + '%' + '\n'
    print(printout)

    # Return new model.
    new_solution_objects = new_solution_objects[1]
    return new_solution_objects

def get_rates_drg(sim_array, solution_object):
    """Calculates values to be used in the calculation of Direct Interaction Coefficients

    Parameters
    ----------
    sim_array :
        Array of simulated simulation objects
    solution_object :
        Cantera object of the solution being reduced

    Returns
    -------
    dict
        Initial conditions and information for calculating DICs at each timestep. The subdictionaries have the
        timestep as their keys and their values hold an array of numberator and denominator information for
        calculating DICs

    """

    old_solution = solution_object
    # iterate through all initial conditions
    total_edge_data = {}
    for ic in sim_array:
        ic_edge_data = {}
        for tstep in ic.sample_points:  # iterate through all timesteps
            temp = tstep[0]  # Set up variables
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
                for i, reac in enumerate(new_solution.reactions()):
                    reac_prod_rate = float(new_reaction_production_rates[i])
                    reactants = reac.reactants
                    products = reac.products
                    all_species = reac.reactants
                    all_species.update(reac.products)
                    if reac_prod_rate != 0:
                        if reac_prod_rate > 0:

                            for species in products:
                                if species in denom:
                                    denom[species] += abs(float(reac_prod_rate * products[species]))
                                else:
                                    denom[species] = abs(float(reac_prod_rate * products[species]))
                                for species_b in all_species:
                                    if species_b != species:
                                        partial_name = species + '_' + species_b
                                        if partial_name in numerator:
                                            numerator[partial_name] += abs(float(reac_prod_rate * products[species]))
                                        else:
                                            numerator[partial_name] = abs(float(reac_prod_rate * products[species]))

                            for species in reactants:
                                if species in denom:
                                    denom[species] += abs(float(reac_prod_rate * reactants[species]))
                                else:
                                    denom[species] = abs(float(reac_prod_rate * reactants[species]))
                                for species_b in all_species:
                                    if species_b != species:
                                        partial_name = species + '_' + species_b
                                        if partial_name in numerator:
                                            numerator[partial_name] += abs(float(reac_prod_rate * reactants[species]))
                                        else:
                                            numerator[partial_name] = abs(float(reac_prod_rate * reactants[species]))

                        if reac_prod_rate < 0:

                            for species in products:
                                if species in denom:
                                    denom[species] += abs(float(reac_prod_rate * products[species]))
                                else:
                                    denom[species] = abs(float(reac_prod_rate * products[species]))
                                for species_b in all_species:
                                    if species_b != species:
                                        partial_name = species + '_' + species_b
                                        if partial_name in numerator:
                                            numerator[partial_name] += abs(float(reac_prod_rate * products[species]))
                                        else:
                                            numerator[partial_name] = abs(float(reac_prod_rate * products[species]))

                            for species in reactants:
                                if species in denom:
                                    denom[species] += abs(float(reac_prod_rate * reactants[species]))
                                else:
                                    denom[species] = abs(float(reac_prod_rate * reactants[species]))
                                for species_b in all_species:
                                    if species_b != species:
                                        partial_name = species + '_' + species_b
                                        if partial_name in numerator:
                                            numerator[partial_name] += abs(float(reac_prod_rate * reactants[species]))
                                        else:
                                            numerator[partial_name] = abs(float(reac_prod_rate * reactants[species]))

            ic_edge_data[temp] = [denom, numerator]
        total_edge_data[ic] = ic_edge_data
    return total_edge_data

def graph_search(nx_graph, target_species):
    """Search nodal graph and generate list of species to remove

    Parameters
    ----------
    nx_graph : obj
        networkx graph object of solution
    target_species : list
        List of target species to search from

    Returns
    -------
    essential_nodes : str
        String containing names of essential species

    """

    if len(target_species) > 1:
        essential_nodes = list()
        for target in target_species:
            essential = list(networkx.dfs_preorder_nodes(nx_graph, target))
            for sp in essential:
                if sp not in essential_nodes:
                    essential_nodes.append(sp)
    else:
        essential_nodes = list(
            networkx.dfs_preorder_nodes(nx_graph, target_species[0]))

    return essential_nodes
