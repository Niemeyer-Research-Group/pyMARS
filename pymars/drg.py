"""Module containing Directed Relation Graph (DRG) reduction method."""
import logging

import networkx
import numpy as np
import cantera as ct

from . import soln2cti
from .sampling import sample, sample_metrics, calculate_error, SamplingInputs
from .reduce_model import trim, ReducedModel


def create_drg_matrix(state, solution):
    """Creates DRG adjacency matrix based on direct interaction coefficients

    Parameters
    ----------
    state : tuple
        Tuple of state with temperature, pressure, and species mass fractions
    solution : cantera.Solution
        Cantera object of the solution being analyzed

    Returns
    -------
    adjacency_matrix : numpy.ndarray
        Adjacency matrix based on calculated direct interaction coefficients

    """
    temp, pressure, mass_fractions = state
    solution.TPY = temp, pressure, mass_fractions

    net_stoich = solution.product_stoich_coeffs() - solution.reactant_stoich_coeffs()
    flags = np.where(((solution.product_stoich_coeffs() != 0) |
                        (solution.reactant_stoich_coeffs() !=0 )
                        ), 1, 0)

    # only consider contributions from reactions with nonzero net rates of progress
    valid_reactions = np.where(solution.net_rates_of_progress != 0)[0]
    if valid_reactions.size:
        base_rates = np.abs(
            net_stoich[:, valid_reactions] *
            solution.net_rates_of_progress[valid_reactions]
            )
        denominator = np.sum(base_rates, axis=1)[:, np.newaxis]

        numerator = np.zeros((solution.n_species, solution.n_species))
        for sp_b in range(solution.n_species):
            numerator[:, sp_b] += np.sum(
                base_rates[:, np.where(flags[sp_b, valid_reactions])[0]], axis=1
                )
        #numerator = np.einsum('ij,kj->ik', base_rates, flags)

        # May get divide by zero if an inert species is present, and denominator
        # entry is zero.
        with np.errstate(divide='ignore', invalid='ignore'):
            adjacency_matrix = np.where(denominator != 0, numerator / denominator, 0)

    else:
        adjacency_matrix = np.zeros((solution.n_species, solution.n_species))

    # set diagonals to zero, to avoid self-directing graph edges
    np.fill_diagonal(adjacency_matrix, 0.0)

    return adjacency_matrix


def graph_search(graph, target_species):
    """Search nodal graph and generate list of species to remove

    Parameters
    ----------
    graph : networkx.DiGraph
        graph representing reaction system
    target_species : list
        List of target species to search from

    Returns
    -------
    reached_species : list of str
        List of species reachable from targets in graph

    """
    reached_species = []
    for target in target_species:
        reached_species += list(networkx.dfs_preorder_nodes(graph, target))
    reached_species = list(set(reached_species))

    return reached_species


def reduce_drg(solution, model_file, species_targets, species_safe, threshold,
               matrices, sample_inputs, sampled_metrics, threshold_upper=None
               ):
    """Given a threshold and DRG matrix, reduce the model and determine the error.

    Parameters
    ----------
    solution : cantera.Solution
        Model being reduced
    model_file : str
        Filename for model being reduced
    species_targets : list of str
        List of target species names
    species_safe : list of str
        List of species to always be retained
    threshold : float
        DRG threshold for trimming graph
    matrices : list of numpy.ndarray
        List of DRG adjacency matrices determined from thermochemical state data
    sample_inputs : SamplingInputs
        Filename information for sampling (e.g., autoignition inputs/outputs)
    sampled_metrics: numpy.ndarray
        Global metrics from original model used to evaluate error
    threshold_upper : float, optional
        Optional upper threshold (epsilon^star) used to identify species for
        further sensitivity analysis

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

    """
    name_mapping = {i: sp for i, sp in enumerate(solution.species_names)}
    species_retained = []
    for matrix in matrices:
        graph = networkx.DiGraph(np.where(matrix >= threshold, matrix, 0.0))
        networkx.relabel_nodes(graph, name_mapping, copy=False)
        species_retained += graph_search(graph, species_targets)
    # want to ensure retained species are the set of those reachable for each state
    species_retained = list(set(species_retained))

    species_removed = [sp for sp in solution.species_names
                       if sp not in (species_retained + species_safe)
                       ]

    # Cut the exclusion list from the model.
    new_solution = trim(solution, species_removed, model_file)
    new_model_file = soln2cti.write(new_solution)

    reduced_model_metrics = sample_metrics(sample_inputs, new_model_file)
    error = calculate_error(sampled_metrics, reduced_model_metrics)
    
    # If desired, now identify limbo species for future sensitivity analysis
    species_limbo = []
    if threshold_upper:
        graph = networkx.DiGraph(np.where(matrix >= threshold_upper, matrix, 0.0))
        networkx.relabel_nodes(graph, name_mapping, copy=False)
        species_retained += graph_search(graph, species_targets)
        species_limbo = [sp for sp in solution.species_names
                         if sp not in (species_retained + species_safe + species_removed)
                         ]

    return ReducedModel(
        model=new_solution, error=error, filename=new_model_file, limbo_species=species_limbo
        )


def run_drg(model_file, sample_inputs, error_limit, species_targets,
            species_safe, threshold_upper=None
            ):
    """Main function for running DRG reduction.

    Parameters
    ----------
    model_file : str
        Original model file
    sample_inputs : SamplingInputs
        Contains filenames for sampling
    error_limit : float
        Maximum allowable error level for reduced model
    species_targets : list of str
        List of target species names
    species_safe : list of str
        List of species names to always be retained
    threshold_upper: float, optional
        Upper threshold (epsilon^*) to identify limbo species for sensitivity analysis

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

    """
    solution = ct.Solution(model_file)

    assert species_targets, 'Need to specify at least one target species.'

    # first, sample thermochemical data and generate metrics for measuring error
    # (e.g, ignition delays). Also produce adjacency matrices for graphs, which
    # will be used to produce graphs for any threshold value.
    sampled_data, sampled_metrics = sample(sample_inputs, model_file)

    matrices = []
    for state in sampled_data:
        matrices.append(create_drg_matrix((state[0], state[1], state[2:]), solution))

    # begin reduction iterations
    logging.info('Beginning DRG reduction loop')
    logging.info(45 * '-')
    logging.info('Threshold | Number of species | Max error (%)')

    iterations = 0
    error_current = 0.0
    threshold = 0.01
    threshold_increment = 0.01
    while error_current <= error_limit:
        reduced_model = reduce_drg(
            solution, model_file, species_targets, species_safe, 
            threshold, matrices, sample_inputs, sampled_metrics,
            threshold_upper
            )
        error_current = reduced_model.error
        num_species = reduced_model.model.n_species

        # reduce threshold if past error limit on first iteration
        if not iterations and error_current > error_limit:
            error_current = 0.0
            threshold /= 10
            threshold_increment /= 10
            if threshold <= 1e-5:
                raise SystemExit(
                    'Threshold value dropped below 1e-5 without producing viable reduced model'
                    )
            logging.info('Threshold value too high, reducing by factor of 10')
            continue
        
        logging.info(f'{threshold:^9} | {num_species:^17} | {error_current:^.2f}')

        threshold += threshold_increment
        model_previous = reduced_model
    
    if error_current > error_limit:
        reduced_model = model_previous
        error_current = reduced_model.error
        num_species = reduced_model.model.n_species
    
    logging.info(45 * '-')
    logging.info('DRG reduction complete.')
    logging.info(f'Skeletal model: {num_species} species and '
                 f'{reduced_model.model.n_reactions} reactions.'
                 )
    logging.info(f'Maximum error: {error_current:.2f}%')
    return reduced_model
