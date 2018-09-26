"""Main module"""
from argparse import ArgumentParser
from os.path import splitext
from warnings import warn

from pymars import pymars
from convert_chemkin_file import convert

parser = ArgumentParser(description='pyMARS: Reduce chemical kinetic models.')

parser.add_argument('-m', '--model',
                    help='input model filename (e.g., mech.cti).',
                    type=str,
                    required=True
                    )
parser.add_argument('--conditions',
                    help='File with list of autoignition initial conditions.',
                    type=str,
                    required=True
                    )
parser.add_argument('-e', '--error',
                    help='Maximum error % for the reduced model.',
                    type=float,
                    )
parser.add_argument('--method',
                    help='skeletal reduction method to use.',
                    type=str,
                    choices=['DRG', 'DRGEP', 'PFA']
                    )
parser.add_argument('--targets',
                    help="List of target species (e.g., 'CH4 O2').",
                    type=str,
                    nargs='+',
                    required=True,
                    )
parser.add_argument('--retained_species',
                    help="List of non-target species to always retain (e.g., 'N2 Ar')",
                    type=str,
                    nargs='*',
                    default=[]
                    )
parser.add_argument('--run_sa',
                    help='run sensitivity analysis after completing another method.',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('--epsilon_star',
                    help='Epsilon^* value used to determine species for sensitivity analysis.',
                    type=float
                    )

# Specifying conversion requires its own set of options
parser.add_argument('--convert',
                    help='Convert files between Cantera and Chemkin formats (.cti <=> .inp)',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('--thermo',
                    help='thermodynamic data filename (only necessary for CHEMKIN files).',
                    type=str
                    )
parser.add_argument('--transport',
                    help='transport data filename (only necessary for CHEMKIN files).',
                    type=str
                    )

args = parser.parse_args()

if args.run_sa and args.epsilon_star is None:
    parser.error('--run_sa requires --epsilon_star.')

if args.convert:
    # Convert model and exit
    convert(args.model, args.thermo, args.transport)
else:
    # Check for Chemkin format and convert if needed
    if splitext(args.model)[1] in ['.inp', '.dat', '.txt']:
        warn('Chemkin file detected; converting before reduction.')
        args.model = convert(args.model, args.thermo_file, args.transport_file)

    pymars(args.model, args.conditions, args.error, args.method, args.targets,
           args.retained_species, args.run_sa, args.epsilon_star
           )
