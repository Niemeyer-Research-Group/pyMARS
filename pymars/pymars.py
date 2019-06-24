"""Contains main driver function for pyMARS program."""
import os
import sys
import logging

from argparse import ArgumentParser


# local imports
from .sampling import SamplingInputs, sample_metrics, check_inputs
from . import soln2cti
from .drgep import run_drgep
from .drg import run_drg
from .pfa import run_pfa
from .sensitivity_analysis import run_sa
from .tools import convert


def main(model_file, conditions, error_limit, method=None, 
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
   
    return reduced_model


def pymars(argv):
    """
    """
    parser = ArgumentParser(description='pyMARS: Reduce chemical kinetic models.')

    parser.add_argument('-m', '--model',
                        help='input model filename (e.g., "mech.cti").',
                        type=str,
                        )
    parser.add_argument('-e', '--error',
                        help='Maximum error percentage for the reduced model.',
                        type=float,
                        )
    parser.add_argument('--method',
                        help='skeletal reduction method to use.',
                        type=str,
                        choices=['DRG', 'DRGEP', 'PFA'],
                        default=None
                        )
    parser.add_argument('--targets',
                        help='List of target species (e.g., "CH4 O2").',
                        type=str,
                        nargs='+',
                        default=[]
                        )
    parser.add_argument('--conditions',
                        help='File with list of autoignition initial conditions.',
                        type=str,
                        default='ignition_input.yaml'
                        )
    parser.add_argument('--retained_species',
                        help='List of non-target species to always retain (e.g., "N2 Ar")',
                        type=str,
                        nargs='*',
                        default=[]
                        )
    parser.add_argument('--sensitivity_analysis',
                        help='Run sensitivity analysis after completing another method.',
                        action='store_true',
                        default=False,
                        )
    parser.add_argument('--sensitivity_type',
                        help='Sensitivity analysis algorithm to be used',
                        choices=['initial', 'greedy'],
                        type=str,
                        default='greedy'
                        )
    parser.add_argument('--upper_threshold',
                        help='Upper threshold value used to determine species for sensitivity analysis.',
                        type=float,
                        default=1.0
                        )
    parser.add_argument('--path',
                        help='Path to directory for writing files.',
                        type=str,
                        default=''
                        )
    parser.add_argument('--num_threads',
                        help=(
                            'Number of CPU cores to use for running simulations in parallel. '
                            'If no number, then use available number of cores minus 1.'
                            ),
                        nargs='?',
                        const=0,
                        default=1,
                        type=int
                        )

    # Specifying conversion requires its own set of options
    parser.add_argument('--convert',
                        help='Convert files between Cantera and Chemkin formats (.cti <=> .inp)',
                        action='store_true',
                        default=False,
                        )
    parser.add_argument('--thermo',
                        help='thermodynamic data filename (only necessary for Chemkin files).',
                        type=str,
                        default=None
                        )
    parser.add_argument('--transport',
                        help='transport data filename (only necessary for Chemkin files).',
                        type=str,
                        default=None
                        )

    parser.add_argument('-V', '--version',
                        action='store_true',
                        help='Show the version of pyMARS and quit')
    
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if args.version:
        from ._version import __version__
        print('pyMARS {version} from {path} ()'.format(
            version=__version__,
            path=os.path.abspath(os.path.dirname(__file__))))
        sys.exit(0)
    
    if not args.model:
        parser.error('A model file is a required input using -m or --model')

    if args.convert:
        # Convert model and exit
        files = convert(args.model, args.thermo, args.transport, args.path)
        if isinstance(files, list):
            logging.info('Converted files: ' + ' '.join(files))
        else:
            logging.info('Converted file: ' + files)
    else:
        if not args.error:
            parser.error('An error limit is required using -e or --error.')

        if not args.method and not args.sensitivity_analysis:
            parser.error('Either --method or --sensitivity_analysis (or both) must be given.')

        if args.method and not args.targets:
            parser.error(
                'At least one target species must be specified for the graph-based '
                'reduction methods.'
                )

        # Check for Chemkin format and convert if needed
        if os.path.splitext(args.model)[1] != '.cti':
            logging.info('Chemkin file detected; converting before reduction.')
            args.model = convert(args.model, args.thermo, args.transport, args.path)

        main(
            args.model, args.conditions, args.error, 
            method=args.method, target_species=args.targets,
            safe_species=args.retained_species, run_sensitivity_analysis=args.sensitivity_analysis, 
            upper_threshold=args.upper_threshold, sensitivity_type=args.sensitivity_type, 
            path=args.path, num_threads=args.num_threads
            )
