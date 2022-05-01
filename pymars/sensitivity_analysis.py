"""Module containing sensitivity analysis reduction stage. """
import logging

import numpy as np
import cantera as ct

from . import soln2cti
from .sampling import sample_metrics, calculate_error, read_metrics
from .reduce_model import trim, ReducedModel

# Taken from http://stackoverflow.com/a/22726782/1569494
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from contextlib import contextmanager
    import shutil
    import tempfile
    import errno

    @contextmanager
    def TemporaryDirectory():
        name = tempfile.mkdtemp()
        try:
            yield name
        finally:
            try:
                shutil.rmtree(name)
            except OSError as e:
                # Reraise unless ENOENT: No such file or directory
                # (ok if directory has already been deleted)
                if e.errno != errno.ENOENT:
                    raise


def evaluate_species_errors(starting_model, metrics, species_limbo,ignition_conditions=[], flame_conditions=[], psr_conditions=[],
                            phase_name='', num_threads=1
                            ):
    """Calculate error induced by removal of each limbo species

    Parameters
    ----------
    starting_model : ReducedModel
        Container with model and file information
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    metrics : numpy.ndarray
        Calculated metrics for starting model, used for evaluating error
    species_limbo : list of str
        List of species to consider removal
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    num_threads : int, optional
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.
    
    Returns
    -------
    species_errors : numpy.ndarray
        Maximum errors induced by removal of each limbo species

    """
    species_errors = np.zeros(len(species_limbo))
    with TemporaryDirectory() as temp_dir:
        for idx, species in enumerate(species_limbo):

            test_model = trim(
                starting_model.filename, [species], f'reduced_model_{species}.cti', 
                phase_name=phase_name
                )
            test_model_file = soln2cti.write(
                test_model, f'reduced_model_{species}.cti', path=temp_dir
                )
            reduced_model_metrics = sample_metrics(
                test_model_file, ignition_conditions, flame_conditions, phase_name=phase_name, 
                num_threads=num_threads
                )
            species_errors[idx] = calculate_error(metrics, reduced_model_metrics)
    
    return species_errors


def run_sa(model_file, starting_error, error_limit, species_safe, ignition_conditions=[],  flame_conditions=[], psr_conditions=[],
           phase_name='', algorithm_type='greedy', species_limbo=[],
           num_threads=1, path=''
           ):
    """Runs a sensitivity analysis to remove species on a given model.
    
    Parameters
    ----------
    model_file : str
        Model being analyzed
    starting_error : float
        Error percentage between the reduced and original models
    ignition_conditions : list of InputIgnition
        List of autoignition initial conditions.
    psr_conditions : list of InputPSR, optional
        List of PSR simulation conditions.
    flame_conditions : list of InputLaminarFlame, optional
        List of laminar flame simulation conditions.
    error_limit : float
        Maximum allowable error level for reduced model
    species_safe : list of str
        List of species names to always be retained
    phase_name : str, optional
        Optional name for phase to load from CTI file (e.g., 'gas'). 
    algorithm_type : {'initial', 'greedy'}
        Type of sensitivity analysis: initial (order based on initial error), or 
        greedy (all species error re-evaluated after each removal)
    species_limbo : list of str, optional
        List of species to consider; if empty, consider all not in ``species_safe``
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
    current_model = ReducedModel(
        model=ct.Solution(model_file, phase_name), error=starting_error, filename=model_file
        )
    
    logging.info(f'Beginning sensitivity analysis stage, using {algorithm_type} approach.')

    # The metrics for the starting model need to be determined or read
    initial_metrics = sample_metrics(
        model_file, ignition_conditions, flame_conditions, reuse_saved=True, phase_name=phase_name,
        num_threads=num_threads, path=path
        )

    if not species_limbo:
        species_limbo = [
            sp for sp in current_model.model.species_names if sp not in species_safe
            ]

    logging.info(53 * '-')
    logging.info('Number of species |  Species removed  | Max error (%)')

    # Need to first evaluate all induced errors of species; for the ``initial`` method,
    # this will be the only evaluation.
    species_errors = evaluate_species_errors(
        current_model,  initial_metrics, species_limbo, ignition_conditions, flame_conditions,
        phase_name=phase_name, num_threads=num_threads
        )

    # Use a temporary directory to avoid cluttering the working directory with
    # all the temporary model files
    with TemporaryDirectory() as temp_dir:
        while species_limbo:
            # use difference between error and current error to find species to remove
            idx = np.argmin(np.abs(species_errors - current_model.error))
            species_errors = np.delete(species_errors, idx)
            species_remove = species_limbo.pop(idx)

            test_model = trim(
                current_model.filename, [species_remove], f'reduced_model_{species_remove}.cti', 
                phase_name=phase_name
                )
            test_model_file = soln2cti.write(
                test_model, output_filename=f'reduced_model_{species_remove}.cti', path=temp_dir
                )

            reduced_model_metrics = sample_metrics(
                test_model_file, ignition_conditions, flame_conditions, phase_name=phase_name, 
                num_threads=num_threads, path=path
                )
            error = calculate_error(initial_metrics, reduced_model_metrics)

            logging.info(f'{test_model.n_species:^17} | {species_remove:^17} | {error:^.2f}')

            # Ensure new error isn't too high
            if error > error_limit:
                break
            else:
                current_model = ReducedModel(model=test_model, filename=test_model_file, error=error)

            # If using the greedy algorithm, now need to reevaluate all species errors
            if algorithm_type == 'greedy':
                species_errors = evaluate_species_errors(
                    current_model, initial_metrics, species_limbo, ignition_conditions, flame_conditions,
                    phase_name=phase_name, num_threads=num_threads
                    )
                if min(species_errors) > error_limit:
                    break
    
    # Final model; may need to rewrite
    reduced_model = ReducedModel(
        model=current_model.model, filename=f'reduced_{current_model.model.n_species}.cti', 
        error=current_model.error
        )
    soln2cti.write(reduced_model.model, reduced_model.filename, path=path)

    logging.info(53 * '-')
    logging.info('Sensitivity analysis stage complete.')
    logging.info(f'Skeletal model: {reduced_model.model.n_species} species and '
                 f'{reduced_model.model.n_reactions} reactions.'
                 )
    logging.info(f'Maximum error: {reduced_model.error:.2f}%')
    return reduced_model
