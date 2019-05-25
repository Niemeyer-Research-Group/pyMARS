"""Module containing sensitivity analysis reduction stage. """
import logging

import numpy as np
import cantera as ct

from . import soln2cti
from .sampling import sample_metrics, calculate_error, read_metrics, SamplingInputs
from .reduce_model import trim, ReducedModel


def evaluate_species_errors(starting_model, sample_inputs, metrics, species_limbo):
    """Calculate error induced by removal of each limbo species

    Parameters
    ----------
    starting_model : ReducedModel
        Container with model and file information
    sample_inputs : SamplingInputs
        Contains filenames for sampling
    metrics : numpy.ndarray
        Calculated metrics for starting model, used for evaluating error
    species_limbo:
        List of species to consider removal
    
    Returns
    -------
    species_errors : numpy.ndarray
        Maximum errors induced by removal of each limbo species

    """
    species_errors = np.zeros(len(species_limbo))
    for idx, species in enumerate(species_limbo):
        test_model = trim(starting_model.model, [species], starting_model.filename)
        test_model_file = soln2cti.write(test_model)
        reduced_model_metrics = sample_metrics(sample_inputs, test_model_file)
        species_errors[idx] = calculate_error(metrics, reduced_model_metrics)
    
    return species_errors


def run_sa(model_file, starting_error, sample_inputs, error_limit, 
           species_safe, algorithm_type='initial', species_limbo=[]):
    """Runs a sensitivity analysis to remove species on a given model.
    
    Parameters
    ----------
    model_file : str
        Model being analyzed
    starting_error : float
        Error percentage between the reduced and original models
    sample_inputs : SamplingInputs
        Contains filenames for sampling
    error_limit : float
        Maximum allowable error level for reduced model
    species_safe : list of str
        List of species names to always be retained
    algorithm_type : {'initial', 'greedy'}
        Type of sensitivity analysis: initial (order based on initial error), or 
        greedy (all species error re-evaluated after each removal)
    species_limbo:
        List of species to consider; if empty, consider all not in ``species_safe``
    
    Returns
    -------
    ReducedModel
        Return reduced model and associated metadata

    """
    current_model = ReducedModel(
        model=ct.Solution(model_file), error=starting_error, filename=model_file
        )

    # Read in already calculated metrics from the starting model
    initial_metrics = read_metrics(sample_inputs)

    if not species_limbo:
        species_limbo = [sp for sp in current_model.model.species_names if sp not in species_safe]

    logging.info(f'Beginning sensitivity analysis stage, using {algorithm_type} approach.')
    logging.info(53 * '-')
    logging.info('Number of species |  Species removed  | Max error (%)')

    # Need to first evaluate all induced errors of species; for the ``initial`` method,
    # this will be the only evaluation.
    species_errors = evaluate_species_errors(
        current_model, sample_inputs, initial_metrics, species_limbo
        )

    while species_limbo:
        # use difference between error and current error to find species to remove
        idx = np.argmin(np.abs(species_errors - current_model.error))
        species_errors = np.delete(species_errors, idx)
        species_remove = species_limbo.pop(idx)

        test_model = trim(current_model.model, [species_remove], current_model.filename)
        test_model_file = soln2cti.write(test_model)

        reduced_model_metrics = sample_metrics(sample_inputs, test_model_file)
        error = calculate_error(initial_metrics, reduced_model_metrics)

        logging.info(f'{test_model.n_species:^17} | {species_remove:^17} | {error:^.2f}')

        # Ensure new error isn't too high
        if error > error_limit:
            break
        else:
            current_model.model = test_model
            current_model.error = error
            current_model.filename = test_model_file

        # If using the greedy algorithm, now need to reevaluate all species errors
        if algorithm_type == 'greedy':
            species_errors = evaluate_species_errors(
                current_model.model, sample_inputs, initial_metrics, species_limbo
                )
            if min(species_errors) > error_limit:
                break
    
    # Final model; may need to rewrite
    reduced_model = current_model
    reduced_model.filename = soln2cti.write(reduced_model.model)

    logging.info(53 * '-')
    logging.info('Sensitivity analysis stage complete.')
    logging.info(f'Skeletal model: {reduced_model.model.n_species} species and '
                 f'{reduced_model.model.n_reactions} reactions.'
                 )
    logging.info(f'Maximum error: {reduced_model.error:.2f}%')
    return reduced_model
