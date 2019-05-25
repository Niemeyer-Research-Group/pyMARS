
import numpy as np
import cantera as ct
import networkx

from . import soln2cti
from .reduce_model import trim, ReducedModel
from .sampling import sample, sample_metrics, calculate_error, SamplingInputs
from .drg import graph_search

def trim_pfa(total_edge_data, solution_object, threshold_value, keeper_list,
             done, target_species, model_file
             ):
    """Determines what species to remove based on direct interaction coefficients.

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


    safe = [] # A list of species that are to be retained for this threshold value

    # Calculate edge weights based on list received from get_rate_data
    # Initial condition
    for ic in total_edge_data.keys(): # For each initial condition
        # Timestep
        for tstep in total_edge_data[ic].keys(): # Set edge values for the graph
            for species in species_objects: # Make graph
                graph.add_node(species.name)
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
            for sp in dic: # Add to safe list if it is not already there.
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

    """Top level function for running PFA, writes reduced file

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

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

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

    final_error[0] = 0 # An integer representing the error introduced in the final simulation.
    done[0] = False

    while not done[0] and error[0] < error_limit: # Run the simulation until nothing else can be cut.
        # Trim at this threshold value and calculate error.
        reduced_model = pfa_loop_control(
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
    reduced_model = pfa_loop_control(
        solution_object, target_species, retained_species, model_file, error, max_t, done, rate_edge_data, ignition_delay_detailed, conditions_array)

    return reduced_model


def pfa_loop_control(solution_object, target_species, retained_species, model_file, 
                     stored_error, threshold, done, rate_edge_data, 
                     ignition_delay_detailed, conditions_array
                     ):
    """Handles the reduction, simulation, and comparison for a single threshold value

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

    # Run detailed mechanism and retain initial conditions
    species_retained = []
    printout = ''
    print('Threshold     Species in Mech      Error')

    # Run DRG and create new reduced solution
    exclusion_list = trim_pfa(
        rate_edge_data, solution_object, threshold, retained_species, done,target_species,model_file) # Find out what to cut from the model
    new_solution = trim(solution_object, exclusion_list, model_file) # Cut the exclusion list from the model.
    species_retained.append(len(new_solution.species()))

    # Simulated reduced solution
    new_sim = helper.setup_simulations(conditions_array, new_solution) # Create simulation objects for reduced model for all conditions
    ignition_delay_reduced = helper.simulate(new_sim) # Run simulations and process results

    if (ignition_delay_detailed.all() == 0): # Ensure that ignition occured
        print("Original model did not ignite.  Check initial conditions.")
        exit()

    # Calculate error
    error = 100 * abs(ignition_delay_reduced - ignition_delay_detailed) / ignition_delay_detailed
    printout += (str(threshold) + '                 ' + str(len(new_solution.species())) + 
                 '              '+  str(round(np.max(error), 2)) +'%' + '\n'
                 )
    print(printout)
    stored_error[0] = round(np.max(error), 2)

    # Return new model
    return new_solution


def get_rates_pfa(sim_array, solution_object):
    """Calculates values for calculation of Direct Interaction Coefficients

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
    """Gets the PA (and CA) values of all species in a given solution.

    Parameters
    ----------
    new_solution : cantera.Solution
        The object representing the cantera model
    new_reaction_production_rates : 
        the production rates associated with the model

    Returns
    -------
    dict
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
    """Gets the PAB (and CAB) values of all species in a given solution.

    Parameters
    ----------
    new_solution : cantera.Solution
        The object representing the cantera model
    new_reaction_production_rates : 
        the production rates associated with the model

    Returns
    -------
    dict
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
            for species_a in all_species:
                for species_b in all_species:
                    if species_a != species_b:
                        full_name = species_a + "_" + species_b

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
    """Gets the rAB_p1 (and rAB_c1) values of all species in a given solution.

    Parameters
    ----------
    new_solution : cantera.Solution
        The object representing the cantera model
    PA : dict
        A dictionary containing the PA values for the reduction
    CA : dict
        A dictionary containing the CA values for the reduction
    PAB : dict
        A dictionary containing the PAB values for the reduction
    CAB : dict
        A dictionary containing the CAB values for the reduction

    Returns
    -------
    rAB_p1 : dict
        rAB_p1 dictionary
    rAB_c1 : dict
        rAB_c1 dictionaries.

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
    """Gets the rAB_p2 (and rAB_c2) values of all species in a given solution.

    Parameters
    ----------
    new_solution : cantera.Solution
        The object representing the cantera model
    rAB_p1 : dict
        A dictionary containing the rAB_p1 values for the reduction
    rAB_c1 : dict
        A dictionary containing the rAB_c1 values for the reduction

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
