"""Contains main driver function for pyMARS program."""
import logging
import os

# local imports
from .sampling import SamplingInputs, sample_metrics, check_inputs
from . import soln2cti
from .drgep import run_drgep
from .drg import run_drg
from .pfa import run_pfa
from .sensitivity_analysis import run_sa
from .tools import convert


def pymars(model_file, conditions, error_limit, method=None, 
           target_species=[], safe_species=[], 
           run_sensitivity_analysis=False, upper_threshold=None, sensitivity_type='greedy',
           path='', num_threads=1
           ):
    """Driver function for reducing a chemical kinetic model.

    Parameters
    ----------
    model_file : str
        Cantera-format model to be reduced (e.g., 'mech.cti').
    conditions : str
        File with list of autoignition initial conditions.
    error_limit : float
        Maximum error % for the reduced model.
    method : {'DRG', 'DRGEP', 'PFA'}, optional
        Skeletal reduction method to use.
    target_species: list of str, optional
        List of target species for graph-based reduction.
    safe_species : list of str, optional
        List of non-target species to always retain.
    run_sensitivity_analysis : bool, optional
        Flag to run sensitivity analysis after completing another method.
    upper_threshold : float, optional
        Upper threshold (epsilon^*) used to determine species for sensitivity analysis 
        in combination with DRG or DRGEP method
    sensitivity_type : {'initial', 'greedy'}, optional
        Type of sensitivity analysis
    path : str
        Path to directory for writing files
    num_threads : int, optional
        Number of CPU threads to use for performing simulations in parallel.
        Optional; default = 1, in which the multiprocessing module is not used.
        If 0, then use the available number of cores minus one. Otherwise,
        use the specified number of threads.

    Examples
    --------
    >>> pymars('gri30.cti', conditions_file, 10.0, 'DRGEP', ['CH4', 'O2'], safe_species=['N2'])

    """

    # TODO: allow specification of other sampling filenames
    sampling_inputs = SamplingInputs(input_ignition=conditions)

    # check validity of input file
    assert check_inputs(sampling_inputs)

    if method in ['DRG', 'DRGEP', 'PFA']:
        assert target_species, 'Need to specify at least one target species for graph-based reduction methods'
    
    if not method and not run_sensitivity_analysis:
        raise ValueError(
            'Either a graph-based method or sensitivity analysis (or both) nmust be specified.'
            )
    
    if method == 'DRG':
        reduced_model = run_drg(
            model_file, sampling_inputs, error_limit, target_species, safe_species, 
            threshold_upper=upper_threshold, num_threads=num_threads, path=path
            )    
    elif method == 'DRGEP':
        reduced_model = run_drgep(
            model_file, sampling_inputs, error_limit, target_species, safe_species, 
            threshold_upper=upper_threshold, num_threads=num_threads, path=path
            )
    elif method == 'PFA':
        reduced_model = run_pfa(
            model_file, sampling_inputs, error_limit, target_species, safe_species,
            num_threads=num_threads, path=path
            )
    
    error = 0.0
    limbo_species = []
    if method in ['DRG', 'DRGEP', 'PFA']:
        model_file = reduced_model.filename
        error = reduced_model.error
        limbo_species = reduced_model.limbo_species

    if run_sensitivity_analysis:
        if not sensitivity_type:
            sensitivity_type = 'greedy'

        reduced_model = run_sa(
            model_file, error, sampling_inputs, error_limit, target_species + safe_species, 
            algorithm_type=sensitivity_type, species_limbo=limbo_species, 
            num_threads=num_threads, path=path
            )
   
    return soln2cti.write(reduced_model.model, 'reduced_model.cti', path=path)
