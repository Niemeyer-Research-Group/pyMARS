"""Module containing Directed Relation Graph with Error Propagation (DRGEP) reduction method.
"""

import numpy as np
import networkx
import cantera as ct

from . import soln2cti
from .sampling import sample, sample_metrics, calculate_error, SamplingInputs
from .reduce_model import trim, ReducedModel
from .dijkstra import ss_dijkstra_path_length_modified


def create_drgep_matrix(state, solution):
    """Creates DRGEP graph adjacency matrix

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
        base_rates = (
            net_stoich[:, valid_reactions] *
            solution.net_rates_of_progress[valid_reactions]
            )
        
        denominator_dest = np.sum(np.maximum(0.0, -base_rates), axis=1)
        denominator_prod = np.sum(np.maximum(0.0, base_rates), axis=1)
        denominator = np.maximum(denominator_prod, denominator_dest)[:, np.newaxis]

        numerator = np.zeros((solution.n_species, solution.n_species))
        for sp_b in range(solution.n_species):
            numerator[:, sp_b] += np.sum(
                base_rates[:, np.where(flags[sp_b, valid_reactions])[0]], axis=1
                )
        numerator = np.abs(numerator)

        # May get divide by zero if an inert species is present, and denominator
        # entry is zero.
        with np.errstate(divide='ignore', invalid='ignore'):
            adjacency_matrix = np.where(denominator != 0, numerator / denominator, 0)

    else:
        adjacency_matrix = np.zeros((solution.n_species, solution.n_species))

    # set diagonals to zero, to avoid self-directing graph edges
    np.fill_diagonal(adjacency_matrix, 0.0)

    return adjacency_matrix


def graph_search_drgep(graph, target_species):
    """Searches graph to generate a dictionary of the greatest paths to all species from one of the targets.

    Parameters
    ----------
    graph : networkx.DiGraph
        Graph representing model
    target_species : list of str
        List of target species to search from

    Returns
    -------
    overall_coefficients : dict
        Overall interaction coefficients; maximum over all paths from all targets to each species

    """
    overall_coefficients = {}
    for target in target_species:
        coefficients = ss_dijkstra_path_length_modified(graph, target)        
        overall_coefficients = {
            sp:max(overall_coefficients.get(sp, 0.0), coefficients[sp]) for sp in coefficients
            }
    
    return overall_coefficients


def get_importance_coeffs(species_names, target_species, matrices):
    """Calculate importance coefficients for all species

    Parameters
    ----------
    species_names : list of str
        Species names
    target_species : list of str
        List of target species
    matrices : list of numpy.ndarray
        List of adjacency matrices

    Returns
    -------
    importance_coefficients : dict
        Maximum coefficients over all sampled states

    """
    importance_coefficients = {}
    name_mapping = {i: sp for i, sp in enumerate(species_names)}
    for matrix in matrices:
        graph = networkx.DiGraph(matrix)
        networkx.relabel_nodes(graph, name_mapping, copy=False)
        coefficients = graph_search_drgep(graph, target_species)

        importance_coefficients = {
            sp:max(importance_coefficients.get(sp, 0.0), coefficients[sp]) for sp in coefficients
            }
    
    return importance_coefficients


def reduce_drgep(solution, model_file, species_safe, threshold,
                 importance_coeffs, sample_inputs, sampled_metrics
               ):
    """Given a threshold and DRGEP coefficients, reduce the model and determine the error.

    Parameters
    ----------
    solution : cantera.Solution
        Model being reduced
    model_file : str
        Filename for model being reduced
    species_safe : list of str
        List of species to always be retained
    threshold : float
        DRG threshold for trimming graph
    importance_coeffs : dict
        Dictionary with species and their overall interaction coefficients.
    sample_inputs : SamplingInputs
        Filename information for sampling (e.g., autoignition inputs/outputs)
    sampled_metrics: numpy.ndarray
        Global metrics from original model used to evaluate error

    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

    """
    species_removed = [sp for sp in solution.species_names
                       if importance_coeffs[sp] < threshold 
                       and sp not in species_safe
                       ]

    # Cut the exclusion list from the model.
    new_solution = trim(solution, species_removed, model_file)
    new_model_file = soln2cti.write(new_solution)

    reduced_model_metrics = sample_metrics(sample_inputs, new_model_file)
    error = calculate_error(sampled_metrics, reduced_model_metrics)

    return ReducedModel(
        model=new_solution, error=error, filename=new_model_file
        )


def run_drgep(model_file, sample_inputs, error_limit, species_targets,
              species_safe, threshold_upper=None
              ):
    """Main function for running DRGEP reduction.
    
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
    threshold_upper : float, optional
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
        matrices.append(create_drgep_matrix((state[0], state[1], state[2:]), solution))

    # For DRGEP, find the overall interaction coefficients for all species 
    # using the maximum over all the sampled states
    importance_coeffs = get_importance_coeffs(
        solution.species_names, target_species, matrices
        )

    # begin reduction iterations
    logging.info('Beginning DRGEP reduction loop')
    logging.info(45 * '-')
    logging.info('Threshold | Number of species | Max error (%)')

    iterations = 0
    error_current = 0.0
    threshold = 0.01
    threshold_increment = 0.01
    while error_current <= error_limit:
        reduced_model = reduce_drgep(
            solution, model_file, species_safe, threshold, importance_coeffs, 
            sample_inputs, sampled_metrics
            )
        error_current = reduced_model.error
        num_species = reduced_model.model.n_species

        # reduce threshold if past error limit on first iteration
        if not iterations and error_current > error_limit:
            error_current = 0.0
            threshold /= 10
            threshold_increment /= 10
            if threshold <= 1e-6:
                raise SystemExit(
                    'Threshold value dropped below 1e-6 without producing viable reduced model'
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

    if threshold_upper:
        for sp in reduced_model.species_names:
            if importance_coeffs[sp] < threshold_upper and (sp not in species_safe):
                reduced_model.limbo_species.append(sp)
    
    logging.info(45 * '-')
    logging.info('DRGEP reduction complete.')
    logging.info(f'Skeletal model: {num_species} species and '
                 f'{reduced_model.model.n_reactions} reactions.'
                 )
    logging.info(f'Maximum error: {error_current:.2f}%')
    return reduced_model
