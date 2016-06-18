#trims reaction mechanism files

import cantera as ct
import os
import sys
import argparse

def readin(args):
    """Function to import data file and identify format.

    Parameters
    ----------
    data_file:
        Local Chemkin or Cantera data file containing mechanism information
    exclusion_list:
        list of species to trim
    Returns
    -------
        Converted mechanism file
        Trimmed Solution Object
        Trimmed Mechanism file
    """
    data_file=args.file

    global solution_objects

    #import working functions
    from create_trimmed_model import create_trimmed_model
    from convert_chemkin_file import convert
    from write_to_cti import write
    from autoignition_module import run_sim
    ext= os.path.splitext(data_file)[1]

    #if no species are to be trimmed
    if args.species == None:
        exclusion_list=[]
    else:
        exclusion_list=[str(item) for item in args.species.split(',')]
    print(exclusion_list)

    thermo_file = args.thermo
    transport_file=args.transport

    if ext == ".cti" or ext == ".xml":
        print("This is an Cantera xml or cti file")
        #trims file
        solution_objects=create_trimmed_model(data_file, exclusion_list)
        write(data_file, solution_objects)
        run_sim(solution_objects, args)

    elif ext == ".inp" or ext == ".dat" or ext == ".txt":
        print("This is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(data_file, thermo_file, transport_file)

        #trims newly converted file
        solution_objects=create_trimmed_model(converted_file_name, \
                                    exclusion_list)
        write(data_file, solution_objects)
        run_sim(solution_objects, args)


    else:
        print("File type not supported")




"""-------------------------------------------------------------------------
Get details from command line
-------------------------------------------------------------------------"""


#gets arguments from terminal
parser=argparse.ArgumentParser(description='pyMARS main: \
            converts and trims mechanism files \n')
parser.add_argument('--file', help='input mechanism file name', type=str)
parser.add_argument('--species',  help='comma separated list input of species\
                                to exclude', type=str)
parser.add_argument('--thermo', help='thermodynamic data file', type=str)
parser.add_argument('--transport', help='transport data file', type=str)
parser.add_argument('--plot', help='plots ignition curve', action="store_true")
parser.add_argument('--points', help='print sim sample points', action="store_true")
parser.add_argument('--writecsv', help='write species data to csv', action="store_true")
parser.add_argument('--writehdf5', help='write species data to hdf5', action="store_true")
args=parser.parse_args()
data_file=args.file


if __name__ == '__main__':
    readin(args)
