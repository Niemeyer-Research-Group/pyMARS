"""Main module"""
from argparse import ArgumentParser
from os.path import splitext
from warnings import warn
import logging

from .pymars import pymars
from .tools import convert

parser = ArgumentParser(description='pyMARS: Reduce chemical kinetic models.')

parser.add_argument('-m', '--model',
                    help='input model filename (e.g., "mech.cti").',
                    type=str,
                    required=True
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
                    )
parser.add_argument('--upper_threshold',
                    help='Upper threshold value used to determine species for sensitivity analysis.',
                    type=float
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

args = parser.parse_args()

if args.convert and not args.error:
    parser.error('An error limit is required for performing model reduction.')

if not args.method and not args.sensitivity_analysis:
    parser.error('Either --method or --sensitivity_analysis (or both) must be given.')

if args.method and not args.targets:
    parser.error('At least one target species must be specified for the graph-based reduction methods.')

if args.convert:
    # Convert model and exit
    files = convert(args.model, args.thermo, args.transport, args.path)
    if isinstance(files, list):
        logging.info('Converted files: ' + ' '.join(files))
    else:
        logging.info('Converted file: ' + files)
else:
    # Check for Chemkin format and convert if needed
    if splitext(args.model)[1] != '.cti':
        logging.info('Chemkin file detected; converting before reduction.')
        args.model = convert(args.model, args.thermo, args.transport, args.path)

    pymars(
        args.model, args.conditions, args.error, 
        method=args.method, target_species=args.targets,
        safe_species=args.retained_species, run_sensitivity_analysis=args.sensitivity_analysis, 
        upper_threshold=args.upper_threshold, sensitivity_type=args.sensitivity_type, 
        path=args.path, num_threads=args.num_threads
        )
