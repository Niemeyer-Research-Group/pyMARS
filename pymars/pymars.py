"""Contains main driver function for pyMARS program."""
import logging

# local imports
from .sampling import SamplingInputs
from . import soln2cti
from .drgep import run_drgep
from .drg import run_drg
from .pfa import run_pfa
from .sensitivity_analysis import run_sa
from .convert_chemkin_file import convert

def pymars(model_file, conditions, error_limit, method, 
           target_species, safe_species=[], 
           run_sensitivity_analysis=False, upper_threshold=None
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
    target_species: list of str
        List of target species for reduction.
    safe_species : list of str, optional
        List of non-target species to always retain.
    run_sensitivity_analysis : bool, optional
        Flag to run sensitivity analysis after completing another method.
    upper_threshold : float, optional
        Upper threshold (epsilon^*) used to determine species for sensitivity analysis

    Examples
    --------
    >>> pymars('gri30.cti', conditions_file, 10.0, 'DRGEP', ['CH4', 'O2'], retained_species=['N2'])

    """

    # TODO: allow specification of other sampling filenames
    sampling_inputs = SamplingInputs(input_ignition=conditions)

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
    
    if method in ['DRG', 'DRGEP', 'PFA']:
        model_file = reduced_model.filename

    if run_sensitivity_analysis and reduced_model.limbo_species:
        reduced_model = run_sa(
            model_file, reduced_model.error, conditions, target_species, 
            safe_species, error_limit, reduced_model.limbo_species
            )
   
    return soln2cti.write(reduced_model.model)
