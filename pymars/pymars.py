"""Contains main driver function for pyMARS program."""
import logging
import os

# local imports
from .sampling import SamplingInputs, sample_metrics
from . import soln2cti
from .drgep import run_drgep
from .drg import run_drg
from .pfa import run_pfa
from .sensitivity_analysis import run_sa
from .convert_chemkin_file import convert


def pymars(model_file, conditions, error_limit, method, 
           target_species=[], safe_species=[], 
           run_sensitivity_analysis=False, upper_threshold=None,
           path=''
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
    method : {'DRG', 'DRGEP', 'PFA'}
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
    path : str
        Path to directory for writing files

    Examples
    --------
    >>> pymars('gri30.cti', conditions_file, 10.0, 'DRGEP', ['CH4', 'O2'], retained_species=['N2'])

    """

    # TODO: allow specification of other sampling filenames
    sampling_inputs = SamplingInputs(input_ignition=conditions)

    if method in ['DRG', 'DRGEP', 'PFA']:
        assert target_species, 'Need to specify target species for graph-based reduction methods'

    if run_sensitivity_analysis:
        assert (upper_threshold, 
            'Need to specify upper threshold (epsilon^*) for sensitivity analysis'
            )
    
    if method == 'DRG':
        reduced_model = run_drg(
            model_file, sampling_inputs, error_limit, target_species, safe_species, upper_threshold
            )
    elif method == 'PFA':
        reduced_model = run_pfa(
            model_file, sampling_inputs, error_limit, target_species, safe_species
            )
    elif method == 'DRGEP':
        reduced_model = run_drgep(
            model_file, sampling_inputs, error_limit, target_species, safe_species, upper_threshold
            )
    
    error = 0.0
    limbo_species = []
    if method in ['DRG', 'DRGEP', 'PFA']:
        model_file = reduced_model.filename
        error = reduced_model.error
        limbo_species = reduced_model.limbo_species
    else:
        # The metrics for the starting model need to be determined
        sample_metrics(sampling_inputs, model_file, save_output=True)

    if run_sensitivity_analysis:
        reduced_model = run_sa(
            model_file, error, sampling_inputs, error_limit, 
            target_species + safe_species, limbo_species
            )
   
    return soln2cti.write(reduced_model.model)
