#!/usr/bin/env python
from home import readin
import argparse


def m(args=None):
	'''
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

	'''
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
    
	parser.add_argument('--run_drg', \
		help='run Direct Relation Graph method to reduce', \
		action="store_true")
    
	parser.add_argument('--convert', \
		help='Only convet selected file from .cti <====> .inp', \
		action="store_true")
    
	parser.add_argument('--run_drgep', \
		help='run Direct Relation Graph with Error Propigation method to reduce', \
		action="store_true")
	
	parser.add_argument('--run_sa', \
		help='run a sensativity analysis after completing another method', \
		action="store_true")

	parser.add_argument('--error', \
                       help='Maximum allowed error indtroducted by the reduced simulation', \
                       type=float)

	parser.add_argument('--target', \
		help='target species for the model reduction.  Comma seperated list', \
		type=str)
	
	parser.add_argument('--ep_star', \
                       help='Ep star value used in sensativity analysis.', \
                       type=float)

	args=parser.parse_args()
	if args.file is not None:
		readin(args)
	if args.file is None:
		#readin(args)
		string = '''
        Python Model Automated Reduction Software (pyMARS)

        Arguments
            --file:
                Input mechanism file (ex. --file=gri30.cti)
            --run_drg:
                Run Direct Relation Graphing model reduction based on
                a given threshold value
            --conditions:
                Text file of initial conditions for autoignition
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
            --run_drgep :
                Run Direct Relation Graphing with Error Propogation model reduction based on
                a given allowed error.
	    --error:
	        A float representing the maximum ammount of error allowed to be introcued to the model (0-100).
	    --target:
		Comma seperated list of species to use as targets for model reduction.
	    --run_sa:
		Run a senativity analysis after completing a method.
	    --ep_star:
		The ep_star value for completing a sensativity analysis.  


        '''
		print(string)
