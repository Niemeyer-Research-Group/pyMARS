from argparse import ArgumentParser

from .pymars import pymars

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
parser.add_argument('--species',
                    help="comma-separated list of species to always retain (e.g., 'CH4,O2')",
                    type=str
                    )
parser.add_argument('--thermo',
                    help='thermodynamic data filename (only necessary for CHEMKIN files).',
                    type=str
                    )
parser.add_argument('--transport',
                    help='transport data filename (only necessary for CHEMKIN files).',
                    type=str
                    )
parser.add_argument('--convert',
                    help='Convert files between Cantera and Chemkin formats (.cti <=> .inp)',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('--run_drg',
                    help='flag to use DRG reduction method.',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('--run_drgep',
                    help='flag to use DRGEP reduction method.',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('--run_pfa',
                    help='flag to use PFA reduction method.',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('--run_sa',
                    help='run sensitivity analysis after completing another method.',
                    action='store_true',
                    default=False,
                    )
parser.add_argument('-e', '--error',
                    help='Maximum error % for the reduced model.',
                    type=float,
                    )
parser.add_argument('--target',
                    help="Comma-separated list of target species (e.g., 'CH4,O2').",
                    type=str,
                    required=True,
                    )
parser.add_argument('--ep_star',
                    help='Epsilon^* value used in sensitivity analysis.',
                    type=float
                    )

args = parser.parse_args()

pymars(args)
