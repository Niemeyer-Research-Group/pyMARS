"""Main module"""
from argparse import ArgumentParser
from os.path import splitext
from warnings import warn

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
                    required=True
                    )
parser.add_argument('--method',
                    help='skeletal reduction method to use.',
                    type=str,
                    choices=['DRG', 'DRGEP', 'PFA'],
                    required=True,
                    )
parser.add_argument('--targets',
                    help='List of target species (e.g., "CH4 O2").',
                    type=str,
                    nargs='+',
                    required=True,
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
parser.add_argument('--run_sa',
                    help='Run sensitivity analysis after completing another method.',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('--epsilon_star',
                    help='Epsilon^* value used to determine species for sensitivity analysis.',
                    type=float
                    )
parser.add_argument('--path',
                    help='Path to directory for writing files.',
                    type=str,
                    default=''
                    )
parser.add_argument('--num_threads',
                    help='Number of CPU cores to use for running simulations in parallel.',
                    default=None,
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
                    default=''
                    )
parser.add_argument('--transport',
                    help='transport data filename (only necessary for Chemkin files).',
                    type=str,
                    default=''
                    )

args = parser.parse_args()

if args.run_sa and args.epsilon_star is None:
    parser.error('--run_sa requires --epsilon_star.')

if args.convert:
    # Convert model and exit
    convert(args.model, args.thermo, args.transport, args.path)
else:
    # Check for Chemkin format and convert if needed
    if splitext(args.model)[1] != '.cti':
        warn('Chemkin file detected; converting before reduction.')
        args.model = convert(args.model, args.thermo, args.transport, args.path)

    pymars(args.model, args.conditions, args.error, args.method, args.targets,
           args.retained_species, args.run_sa, args.epsilon_star, args.path,
           args.num_threads
           )
