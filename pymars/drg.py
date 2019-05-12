"""Module containing Directed Relation Graph (DRG) reduction method."""
from collections import Counter
import time
import os
import sys

import networkx
import numpy as np
import cantera as ct

from . import soln2ck
from . import soln2cti
from . import helper
from .simulation import Simulation
from .create_trimmed_model import trim
from .readin_initial_conditions import readin_conditions


def trim_drg(total_edge_data, solution_object, threshold_value, keeper_list, done, target_species):
    """Determines which species to remove based on direct interaction coefficients

    Compares direct interaction coefficients (DICs) to a threshold value,
    and performs a simple graph search to find species reachable from targets.

    Parameters
    ----------
    total_edge_data : dict
        Information for calculating the DICs for the graph edge weights
    solution_object : cantera.Solution
        The solution being reduced
    threshold_value : float
        User specified threshold value
    keeper_list : list of str
        Species that should always be kept
    done : bool
        Determines whether or not the reduction is complete
    target_species : list of str
        The target species for the search in the array

    Returns
    -------
    Array of species that should be reduced at this threshold level

    """

    # Initalize solution and components
    solution = solution_object
    species_objects = solution.species()

    # Use the networkx library to create a weighted graph of all of the species and their dependencies on each other.
    graph = networkx.DiGraph()

    # A list of species that are to be retained for this threshold value
    safe = []

    # Calculate edge weights based on list received from get_rate_data
    # Initial condition
    for ic in total_edge_data.keys():  # For each initial condition
        # Timestep
        # Set edge values for the graph
        for tstep in total_edge_data[ic].keys():
            for species in species_objects:  # Make graph
                graph.add_node(species.name)
            numerator = total_edge_data[ic][tstep][1]
            denominator = total_edge_data[ic][tstep][0]
            # Each species
            for edge in numerator:
                try:
                    edge_name = edge.split('_', 1)
                    species_a_name = edge_name[0]
                    species_b_name = edge_name[1]

                    # DRG weight between two species
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
            # Add to the safe list if it is not already there.
            for sp in dic:
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


def run_drg(solution_object, conditions_file, error_limit, target_species,
            retained_species, model_file, final_error, epsilon_star=.1):
    """Main function for running DRG reduction.

    Parameters
    ----------
    solution_object : ~cantera.Solution
        A Cantera object of the solution to be reduced.
    conditions_file : str
        Name of file with list of autoignition initial conditions.
    error_limit : float
        Maximum allowable error level for reduced model.
    target_species : list of str
        List of target species
    retained_species : list of str
        List of species to always be retained
    model_file : string
        The path to the file where the solution object was generated from
    final_error: singleton float
        To hold the error level of simulation
    epsilon_star: float
        Epsilon star value for sensativity analysis

    Returns
    -------
    Tuple of reduced Cantera solution object [0] and list of limbo species for SA [1]

    """

    assert target_species, 'Need to specify at least one target species.'

    # Singleton to hold whether any more species can be cut from the simulation.
    done = []
    done.append(False)
    threshold = 0.1  # Starting threshold value
    threshold_increment = 0.1
    num_iterations = 1
    error = [10.0]

    conditions_array = readin_conditions(conditions_file)

    # Turn conditions array into unrun simulation objects for the original solution
    sim_array = helper.setup_simulations(conditions_array, solution_object)
    # Run simulations and process results
    ignition_delay_detailed = helper.simulate(sim_array)
    # Get edge weight calculation data.
    rate_edge_data = get_rates_drg(sim_array, solution_object)

    print("Testing for starting threshold value")
    # Trim the solution at that threshold and find the error.
    drg_loop_control(
        solution_object, target_species, retained_species, model_file, error, threshold, 
        done, rate_edge_data, ignition_delay_detailed, conditions_array
        )

    # While the error for trimming with that threshold value is greater than allowed.
    while error[0] != 0 and threshold_increment > .001:
        # Reduce the starting threshold value and try again.
        threshold /= 10
        threshold_increment /= 10
        num_iterations += 1
        drg_loop_control(
            solution_object, target_species, retained_species, model_file,
            error, threshold, done, rate_edge_data,
            ignition_delay_detailed, conditions_array
            )
        if error[0] <= 0.02:
            error[0] = 0

    print("Starting with a threshold value of " + str(threshold))
    sol_new = solution_object
    final_error[0] = 0
    done[0] = False

    # Run the simulation until nothing else can be cut.
    while not done[0] and error[0] < error_limit:
        # Trim at this threshold value and calculate error.
        sol_new = drg_loop_control(
            solution_object, target_species, retained_species, model_file,
            error, threshold, done, rate_edge_data,
            ignition_delay_detailed, conditions_array
            )
        # If a new max species cut without exceeding what is allowed is reached, save that threshold.
        if error_limit > error[0]:
            max_t = threshold
            final_error[0] = error[0]
        # if (final_error[0] == error[0]): #If error wasn't increased, increase the threshold at a higher rate.
        #	threshold = threshold + (threshold_increment * 4)
        # if (threshold >= .01):
        #        threshold_increment = .01
        threshold += threshold_increment
        threshold = round(threshold, num_iterations)

    limbo = []
    if epsilon_star:
        print("Calculating for DRGASA:")

        # Trim with ep star as threshold value and calculate error.
        epstar_sol = drg_loop_control(
            solution_object, target_species, retained_species, model_file,
            error, epsilon_star, done, rate_edge_data,
            ignition_delay_detailed, conditions_array
            )

        # Anything that was reduced at epstar threshold, add to limbo
        for sp in sol_new.species_names:
            if sp not in epstar_sol.species_names:
                limbo.append(sp)


    print("Greatest result: ")
    sol_new = drg_loop_control(
        solution_object, target_species, retained_species, model_file,
        error, max_t, done, rate_edge_data,
        ignition_delay_detailed, conditions_array
        )

    result = [sol_new, limbo]
    return result


def drg_loop_control(solution_object, target_species, retained_species, model_file, 
                     stored_error, threshold, done, rate_edge_data,
                     ignition_delay_detailed, conditions_array
                     ):
    """Handles the reduction, simulation, and comparision for a single threshold value.

    Parameters
    ----------
    solution_object : cantera.Solution
        object being reduced
    target_species : list of str
        List of target species
    retained_species : list of str
        List of species to always be retained
    model_file : string
        The path to the file where the solution object was generated from
    stored_error: signleton float
        Error of this reduced model simulation
    threshold : float
        current threshold value
    done : bool
        are we done reducing yet?
    rate_edge_data : dict
        information for calculating the DICs for reduction
    ignition_delay_detailed : numpy.array
        ignition delay of detailed model
    conditions_array : list of Condition
        list holding information about initial conditions

    Returns
    -------
    Reduced solution object for this threshold and updates error value

    """

    species_retained = []
    printout = ''
    print('Threshold     Species in Mech      Error')

    # Run DRG and create new reduced solution
    # Find out what to cut from the model
    exclusion_list = trim_drg(rate_edge_data, solution_object, threshold,
                              retained_species, done, target_species
                              )

    # Cut the exclusion list from the model.
    new_solution = trim(solution_object, exclusion_list, model_file)
    species_retained.append(len(new_solution.species()))

    # Simulated reduced solution
    # Create simulation objects for reduced model for all conditions
    new_sim = helper.setup_simulations(conditions_array, new_solution)
    # Run simulations and process results
    ignition_delay_reduced = helper.simulate(new_sim)

    if ignition_delay_detailed.all() == 0:  # Ensure that ignition occured
        print("Original model did not ignite.  Check initial conditions.")
        exit()

    # Calculate error
    error = 100 * abs(ignition_delay_reduced - ignition_delay_detailed) / ignition_delay_detailed
    printout += (str(threshold) + '                 ' + str(len(new_solution.species())) + 
                 '              ' + str(round(np.max(error), 2)) + '%' + '\n'
                 )
    print(printout)
    stored_error[0] = round(np.max(error), 2)

    # Return new model.
    return new_solution


def get_rates_drg(sim_array, solution_object):
    """Calculates values to be used in the calculation of Direct Interaction Coefficients

    Parameters
    ----------
    sim_array : 
        Array of simulated simulation objects
    solution_object : cantera.Solution
        Cantera object of the solution being reduced

    Returns
    -------
    dict
        Initial conditions and information for calculating DICs at each timestep. 
        The subdictionaries have the timestep as their keys and their values hold an 
        array of numberator and denominator information for calculating DICs.

    """

    old_solution = solution_object
    # Iterate through all initial conditions
    total_edge_data = {}
    for ic in sim_array:
        ic_edge_data = {}
        for tstep in ic.sample_points:  # Iterate through all timesteps
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
    """Search nodal graph and generate list of species to remove.

    Parameters
    ----------
    nx_graph : networkx.Graph
        Object of solution graph
    target_species : list of str
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
