"""Module containing Directed Relation Graph (DRG) reduction method."""
import logging
import os
import networkx
import numpy as np
import cantera as ct

from . import soln2cti
from .sampling import sample, sample_metrics, calculate_error
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

    net_stoich = solution.product_stoich_coeffs - solution.reactant_stoich_coeffs
    flags = np.where(((solution.product_stoich_coeffs != 0) |
                        (solution.reactant_stoich_coeffs !=0 )
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


def trim_drg(matrix, species_names, species_targets, threshold):
    """

    Parameters
    ----------
    matrix : np.ndarray
        Adjacency matrix representing graph
    species_names : list of str
        List of all species names
    species_targets : list of str
        List of target species names
    threshold : float
        DRG threshold for trimming graph

    Returns
    ------
    species_reached : list of str
        Names of species reached in graph search

    """
    name_mapping = {i: sp for i, sp in enumerate(species_names)}
    graph = networkx.DiGraph(np.where(matrix >= threshold, matrix, 0.0))
    networkx.relabel_nodes(graph, name_mapping, copy=False)
    species_reached = graph_search(graph, species_targets)

    return species_reached


def reduce_drg(model_file, species_targets, species_safe, threshold,
               matrices, ignition_conditions, sampled_metrics,
               phase_name='', previous_model=None, threshold_upper=None,
               num_threads=1, path=''
               ):
    """Given a threshold and DRG matrix, reduce the model and determine the error.

    Parameters
    ----------
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
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    sampled_metrics: numpy.ndarray
        Global metrics from original model used to evaluate error
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas').
    previous_model : ReducedModel, optional
        Model produced at previous threshold level; used to avoid repeated work.
    threshold_upper : float, optional
        Optional upper threshold (epsilon^star) used to identify species for
        further sensitivity analysis
    num_threads : int, optional
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    path : str, optional
        Optional path for writing files

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

    """
    solution = ct.Solution(model_file, phase_name)

    species_retained = []
    for matrix in matrices:
        species_retained += trim_drg(matrix, solution.species_names, species_targets, threshold)

    # want to ensure retained species are the set of those reachable for each state
    species_retained = list(set(species_retained))

    if previous_model and len(species_retained) == previous_model.model.n_species:
        return previous_model

    species_removed = [sp for sp in solution.species_names
                       if sp not in (species_retained + species_safe)
                       ]

    # Cut the exclusion list from the model.
    reduced_model = trim(
        model_file, species_removed, f'reduced_{model_file}', phase_name=phase_name
        )
    reduced_model_filename = soln2cti.write(
        reduced_model, f'reduced_{reduced_model.n_species}.cti', path=path
        )

    reduced_model_metrics = sample_metrics(
        reduced_model_filename, ignition_conditions, phase_name=phase_name,
        num_threads=num_threads, path=path
        )
    error = calculate_error(sampled_metrics, reduced_model_metrics)

    # If desired, now identify limbo species for future sensitivity analysis
    limbo_species = []
    if threshold_upper:
        species_retained = []
        for matrix in matrices:
            species_retained += trim_drg(
                matrix, solution.species_names, species_targets, threshold_upper
                )
        limbo_species = list(set(species_retained))

    return ReducedModel(
        model=reduced_model, filename=reduced_model_filename,
        error=error, limbo_species=limbo_species
        )


def run_drg(model_file, ignition_conditions, psr_conditions, flame_conditions,
            error_limit, species_targets, species_safe, phase_name='',
            threshold_upper=None, num_threads=1, path=''
            ):
    """Main function for running DRG reduction.

    Parameters
    ----------
    model_file : str
        Original model file
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    psr_conditions : list of InputPSR, optional
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame, optional
        List of laminar flame simulation conditions.
    error_limit : float
        Maximum allowable error level for reduced model
    species_targets : list of str
        List of target species names
    species_safe : list of str
        List of species names to always be retained
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas').
    threshold_upper: float, optional
        Upper threshold (epsilon^*) to identify limbo species for sensitivity analysis
    num_threads : int, optional
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    path : str, optional
        Optional path for writing files

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

    """
    solution = ct.Solution(model_file, phase_name)

    assert species_targets, 'Need to specify at least one target species.'

    # first, sample thermochemical data and generate metrics for measuring error
    # (e.g, ignition delays). Also produce adjacency matrices for graphs, which
    # will be used to produce graphs for any threshold value.
    sampled_metrics, sampled_data = sample(
        model_file, ignition_conditions, phase_name=phase_name, num_threads=num_threads, path=path
        )

    matrices = []
    for state in sampled_data:
        matrices.append(create_drg_matrix((state[0], state[1], state[2:]), solution))

    # begin reduction iterations
    logging.info('Beginning DRG reduction loop')
    logging.info(45 * '-')
    logging.info('Threshold | Number of species | Max error (%)')

    # start with detailed (starting) model
    previous_model = ReducedModel(model=solution, filename=model_file, error=0.0)

    first = True
    error_current = 0.0
    threshold = 0.01
    threshold_increment = 0.01
    while error_current <= error_limit:
        reduced_model = reduce_drg(
            model_file, species_targets, species_safe, threshold, matrices,
            ignition_conditions, sampled_metrics,
            phase_name=phase_name, previous_model=previous_model,
            threshold_upper=threshold_upper, num_threads=num_threads, path=path
            )
        error_current = reduced_model.error
        num_species = reduced_model.model.n_species

        # reduce threshold if past error limit on first iteration
        if first and error_current > error_limit:
            error_current = 0.0
            threshold /= 10
            threshold_increment /= 10
            if threshold <= 1e-5:
                raise SystemExit(
                    'Threshold value dropped below 1e-5 without producing viable reduced model'
                    )
            logging.info('Threshold value too high, reducing by factor of 10')
            continue

        logging.info(f'{threshold:^9.2e} | {num_species:^17} | {error_current:^.2f}')

        threshold += threshold_increment
        first = False

        # cleanup files
        if previous_model.model.n_species != reduced_model.model.n_species:
            os.remove(reduced_model.filename)

        previous_model = ReducedModel(
            model=reduced_model.model, filename=reduced_model.filename,
            error=reduced_model.error, limbo_species=reduced_model.limbo_species
            )

    if error_current > error_limit:
        threshold -= (2 * threshold_increment)
        reduced_model = reduce_drg(
            model_file, species_targets, species_safe, threshold, matrices,
            ignition_conditions, sampled_metrics, phase_name=phase_name,
            threshold_upper=threshold_upper, num_threads=num_threads, path=path
            )
    else:
        soln2cti.write(reduced_model, f'reduced_{reduced_model.model.n_species}.cti', path=path)

    logging.info(45 * '-')
    logging.info('DRG reduction complete.')
    logging.info(f'Skeletal model: {reduced_model.model.n_species} species and '
                 f'{reduced_model.model.n_reactions} reactions.'
                 )
    logging.info(f'Maximum error: {reduced_model.error:.2f}%')
    logging.info('Final reduced model saved as ' + reduced_model.filename)
    return reduced_model
