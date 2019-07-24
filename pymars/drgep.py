"""Module containing Directed Relation Graph with Error Propagation (DRGEP) reduction method.
"""
import os
import logging
from collections import deque
from heapq import heappush, heappop
from itertools import count

import numpy as np
import networkx
import cantera as ct

from . import soln2cti
from .sampling import sample, sample_metrics, calculate_error
from .reduce_model import trim, ReducedModel


def mod_dijkstra(G, source, get_weight, pred=None, paths=None, 
                 cutoff=None, target=None
                 ):
    """Modified implementation of Dijkstra's algorithm for DRGEP method.
    
    Multiples values along graph pathways instead of adding and returns a dictionary 
    with nodes as keys and values containing the greatest path to that node. 
    Each edge weight must be <= 1 so that the further they are from the source, 
    the less important they are.

    Parameters
    ----------
    G : networkx.Graph
        Graph to be considered
    source : str or int
       Starting node for path
    get_weight: function
        Function for getting edge weight
    pred: list, optional
        List of predecessors of a node
    paths: dict, optional
        Path from the source to a target node.
    target : str or int, optional
       Ending node for path
    cutoff : int or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    distance, path : dict
       Returns a tuple of two dictionaries keyed by node.
       The first dictionary stores distance from the source.
       The second stores the path from the source to that node.
    pred, distance : dict
       Returns two dictionaries representing a list of predecessors
       of a node and the distance to each node.
    distance : dict
       Dictionary of greatest lengths keyed by target.

    """
    G_succ = G.succ if G.is_directed() else G.adj #Adjaceny list
    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {source: 0}
    c = count()
    fringe = []  # use heapq with (distance,label) tuples
    push(fringe, (0, next(c), source))
    while fringe:
        cont = 0
        (d, _, v) = pop(fringe)
        if v in dist and d < dist[v]:
            continue  # already searched this node.
        if v == source:
            d = 1
        dist[v] = d 
        if v == target:
            break
        for u, e in G_succ[v].items(): #For all adjancent edges, get weights and multiply them by current path taken. 
            cost = get_weight(v, u, e)
            if cost is None:
                continue
            vu_dist = dist[v] * get_weight(v, u, e)
            if cutoff is not None:
                if vu_dist < cutoff:
                    continue
            #if v in dist:
                #if vu_dist > dist[v]:
                    #raise ValueError('Contradictory paths found:',
                                    #'weights greater than one?')
            elif u not in seen or vu_dist > seen[u]: #If this path to u is greater than any other path we've seen, push it to the heap to be added to dist.
                seen[u] = vu_dist
                push(fringe, (vu_dist, next(c), u))
                if paths is not None:
                    paths[u] = paths[v] + [u]
                if pred is not None:
                    pred[u] = [v]
            elif vu_dist == seen[u]:
                if pred is not None:
                    pred[u].append(v)

    if paths is not None:
        return (dist, paths)
    if pred is not None:
        return (pred, dist)
    return dist


def ss_dijkstra_path_length_modified(G, source, cutoff=None, weight='weight'):
    """Compute the greatest path length via multiplication between source and all other
    reachable nodes for a weighted graph with all weights <= 1.

    Parameters
    ----------
    G : networkx.Graph
        Graph to be considered
    source : node label
       Starting node for path
    weight: str, optional (default='weight')
       Edge data key corresponding to the edge weight.
    cutoff : int or float, optional
       Depth to stop the search. Only paths of length <= cutoff are returned.

    Returns
    -------
    length : dictionary
       Dictionary of shortest lengths keyed by target.

    Examples
    --------
    >>> G=networkx.path_graph(5)
    >>> length=networkx.ss_dijkstra_path_length_modified(G,0)
    >>> length[4]
    1
    >>> print(length)
    {0: 0, 1: 1, 2: 1, 3: 1, 4: 1}

    Notes
    -----
    Edge weight attributes must be numerical and <= 1.
    Distances are calculated as products of weighted edges traversed.
    Don't use a cutoff.

    See Also
    --------
    single_source_dijkstra()

    """
    if G.is_multigraph():
        get_weight = lambda u, v, data: min(
            eattr.get(weight, 1) for eattr in data.values())
    else:
        get_weight = lambda u, v, data: data.get(weight, 1)

    return mod_dijkstra(G, source, get_weight, cutoff=cutoff)


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
        # ensure target has importance of 1.0
        coefficients[target] = 1.0

        for sp in coefficients:
            overall_coefficients[sp] = max(overall_coefficients.get(sp, 0.0), coefficients[sp])
    
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
    importance_coefficients = {sp:0.0 for sp in species_names}
    name_mapping = {i: sp for i, sp in enumerate(species_names)}
    for matrix in matrices:
        graph = networkx.DiGraph(matrix)
        networkx.relabel_nodes(graph, name_mapping, copy=False)
        coefficients = graph_search_drgep(graph, target_species)
        
        importance_coefficients = {
            sp:max(coefficients.get(sp, 0.0), importance_coefficients[sp]) 
            for sp in importance_coefficients
        }
    
    return importance_coefficients


def reduce_drgep(model_file, species_safe, threshold, importance_coeffs, ignition_conditions, 
                 sampled_metrics, phase_name='', previous_model=None, num_threads=1, path=''
                 ):
    """Given a threshold and DRGEP coefficients, reduce the model and determine the error.

    Parameters
    ----------
    model_file : str
        Filename for model being reduced
    species_safe : list of str
        List of species to always be retained
    threshold : float
        DRG threshold for trimming graph
    importance_coeffs : dict
        Dictionary with species and their overall interaction coefficients.
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    sampled_metrics: numpy.ndarray
        Global metrics from original model used to evaluate error
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    previous_model : ReducedModel, optional
        Model produced at previous threshold level; used to avoid repeated work.
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
    species_removed = [sp for sp in solution.species_names
                       if importance_coeffs[sp] < threshold 
                       and sp not in species_safe
                       ]
    
    if (previous_model and 
        len(species_removed) == solution.n_species - previous_model.model.n_species
        ):
        return previous_model

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

    return ReducedModel(
        model=reduced_model, filename=reduced_model_filename, error=error
        )


def run_drgep(model_file, ignition_conditions, psr_conditions, flame_conditions, 
              error_limit, species_targets, species_safe, phase_name='',
              threshold_upper=None, num_threads=1, path=''
              ):
    """Main function for running DRGEP reduction.
    
    Parameters
    ----------
    model_file : str
        Original model file
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    psr_conditions : list of InputPSR
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame
        List of laminar flame simulation conditions.
    error_limit : float
        Maximum allowable error level for reduced model
    species_targets : list of str
        List of target species names
    species_safe : list of str
        List of species names to always be retained
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    threshold_upper : float, optional
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
        matrices.append(create_drgep_matrix((state[0], state[1], state[2:]), solution))

    # For DRGEP, find the overall interaction coefficients for all species 
    # using the maximum over all the sampled states
    importance_coeffs = get_importance_coeffs(
        solution.species_names, species_targets, matrices
        )

    # begin reduction iterations
    logging.info('Beginning DRGEP reduction loop')
    logging.info(45 * '-')
    logging.info('Threshold | Number of species | Max error (%)')

    # start with detailed (starting) model
    previous_model = ReducedModel(model=solution, filename=model_file, error=0.0)

    first = True
    error_current = 0.0
    threshold = 0.01
    threshold_increment = 0.01
    while error_current <= error_limit:
        reduced_model = reduce_drgep(
            model_file, species_safe, threshold, importance_coeffs, ignition_conditions, 
            sampled_metrics, phase_name=phase_name, previous_model=previous_model, 
            num_threads=num_threads, path=path
            )
        error_current = reduced_model.error
        num_species = reduced_model.model.n_species

        # reduce threshold if past error limit on first iteration
        if first and error_current > error_limit:
            error_current = 0.0
            threshold /= 10
            threshold_increment /= 10
            if threshold <= 1e-6:
                raise SystemExit(
                    'Threshold value dropped below 1e-6 without producing viable reduced model'
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
    
    if reduced_model.error > error_limit:
        threshold -= (2 * threshold_increment)
        reduced_model = reduce_drgep(
            model_file, species_safe, threshold, importance_coeffs, ignition_conditions, 
            sampled_metrics, phase_name=phase_name, num_threads=num_threads, path=path
            )
    else:
        soln2cti.write(reduced_model, f'reduced_{reduced_model.model.n_species}.cti', path=path)

    if threshold_upper:
        for sp in reduced_model.model.species_names:
            if importance_coeffs[sp] < threshold_upper and (sp not in species_safe):
                reduced_model.limbo_species.append(sp)
    
    logging.info(45 * '-')
    logging.info('DRGEP reduction complete.')
    logging.info(f'Skeletal model: {reduced_model.model.n_species} species and '
                 f'{reduced_model.model.n_reactions} reactions.'
                 )
    logging.info(f'Maximum error: {reduced_model.error:.2f}%')
    logging.info('Final reduced model saved as ' + reduced_model.filename)
    return reduced_model
