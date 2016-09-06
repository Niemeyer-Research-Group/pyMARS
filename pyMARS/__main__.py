#!/usr/bin/env python
from pyMARS import readin
import argparse



"""-------------------------------------------------------------------------
Get details from command line
-------------------------------------------------------------------------"""
def main(args=None):
    """
    Arguments
        --file: Input mechanism file (ex. --file=gri30.cti)
        --species: Species to eliminate (ex. --species='H, OH')
        --thermo: Thermo data file if Chemkin format (ex. --thermo= thermo.dat)
        --transport: Transport data file if Chemkin format
        --plot
        --points
        --writecsv
        --writehdf5
        --run_drg

    """
    #gets arguments from terminal
    parser=argparse.ArgumentParser(description='pyMARS main: \
                converts and trims mechanism files \n')
    parser.add_argument('--file', \
                        help='input mechanism file name', \
                        type=str)

    parser.add_argument('--species',  \
                        help='comma separated list input of species to exclude',\
                        type=str)

    parser.add_argument('--thermo', \
                        help='thermodynamic data file', \
                        type=str)

    parser.add_argument('--transport', \
                        help='transport data file', \
                        type=str)

    parser.add_argument('--plot', \
                        help='plots ignition curve', \
                        action="store_true")

    parser.add_argument('--points', \
                        help='print sim sample points', \
                        action="store_true")

    parser.add_argument('--writecsv', \
                        help='write species data to csv', \
                        action="store_true")

    parser.add_argument('--writehdf5', \
                        help='write species data to hdf5', \
                        action="store_true")
    parser.add_argument('--run_drg', \
                        help='run Direct Relation Graph method to reduce', \
                        action="store_true")
    parser.add_argument('--iterate', \
                        help='Iterate DRG up to acceptable error limit', \
                        action="store_true")


    args=parser.parse_args()
    if args.file is not None:
        readin(args)
    if args.file is None:
        #readin(args)
        print """
        Python Model Automated Reduction Software (pyMARS)

        Arguments
            --file: Input mechanism file (ex. --file=gri30.cti)
            --species: Species to eliminate (ex. --species='H, OH')
            --thermo: Thermo data file if Chemkin format(ex.--thermo= thermo.dat)
            --transport: Transport data file if Chemkin format
            --plot
            --points
            --writecsv
            --writehdf5
            --run_drg
            --iterate: Insead, iterate DRG up to acceptable error limit

        """
