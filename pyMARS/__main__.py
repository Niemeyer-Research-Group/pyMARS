#!/usr/bin/env python
from pyMARS import readin
import argparse


def main(args=None):
    """
    Parameters
    ----------
        file : str
            Input mechanism file (ex. --file=gri30.cti)
        conditions_file : str
            File of initial conditions for autoignition
        species : str
            Species to eliminate (ex. --species='H, OH')
        thermo : str
            Thermo data file if Chemkin format (ex. --thermo= thermo.dat)
        transport : str
            Transport data file if Chemkin format
        conditions :
            Text file of initial conditions for autoignition
        thresholds :
            csv file containing threshold values to test (usr prompted otherwise)
	error:
	    A float representing the maximum ammount of error allowed to be introcued to the model (0-100).
        run_drg :
            Run Direct Relation Graphing model reduction based on
            a given threshold value
        run_drgep :
            Run Direct Relation Graphing with Error Propogation model reduction based on
            a given allowed error.

        plot :
            Plot a temperature profile of autoignition
        points :
            Return sampling and ignition points
        writecsv :
            Write autoignition data to a csv file
        writehdf5 :
            Write autoignition to a hdf5 file

    """
    #gets arguments from terminal
    parser=argparse.ArgumentParser(description='pyMARS main: \
                converts and trims mechanism files \n')
    parser.add_argument('--file', \
                        help='input mechanism file name', \
                        type=str)
    parser.add_argument('--conditions', \
                        help='initial conditions file name', \
                        type=str)

    parser.add_argument('--species',  \
                        help='comma separated list input of species to never exclude',\
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
    parser.add_argument('--thresholds', \
                        help='csv file containing threshold values to test (usr prompted otherwise)', \
                        type=str)
    parser.add_argument('--convert', \
                        help='Only convet selected file from .cti <====> .inp', \
                        action="store_true")
    parser.add_argument('--run_drgep', \
			help='run Direct Relation Graph with Error Propigation method to reduce', \
			action="store_true")
    parser.add_argument('--error', \
                        help='Maximum allowed error indtroducted by the reduced simulation', \
                        type=float)
    parser.add_argument('--target', \
			help='target species for the model reduction.  Comma seperated list.', \
			type=str)

    args=parser.parse_args()
    if args.file is not None:
        readin(args)
    if args.file is None:
        #readin(args)
        print """
        Python Model Automated Reduction Software (pyMARS)

        Arguments
            --file:
                Input mechanism file (ex. --file=gri30.cti)
            --run_drg:
                Run Direct Relation Graphing model reduction based on
                a given threshold value
            --conditions:
                Text file of initial conditions for autoignition
            --thresholds:
                csv file containing threshold values to test (usr prompted otherwise)
            --convert:
                Convert .cti <==> .inp
            --thermo:
                Thermo data file if Chemkin format (ex. --thermo= thermo.dat)
            --transport:
                Transport data file if Chemkin format
            --species:
                Specific species to not eliminate (ex. --species='H, OH')
            --plot:
                Plot a temperature profile of autoignition
            --points:
                Return sampling and ignition points
            --writecsv:
                Write autoignition data to a csv file
            --writehdf5:
                Write autoignition to a hdf5 file
            --run_drgep :
                Run Direct Relation Graphing with Error Propogation model reduction based on
                a given allowed error.
	    --error:
	        A float representing the maximum ammount of error allowed to be introcued to the model (0-100).
	    --target:
		Comma seperated list of species to use as targets for model reduction.


        """
main()
