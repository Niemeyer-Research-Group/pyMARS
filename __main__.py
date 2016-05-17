#trims reaction mechanism files

import cantera as ct
import os
import sys
import argparse

def readin(data_file, exclusion_list, thermo_file, transport_file):
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
    print(exclusion_list)
    global solution_objects
    #import working functions
    from create_trimmed_model import create_trimmed_model
    from convert_chemkin_file import convert
    from write_to_cti import write

    if data_file.endswith(".xml") or data_file.endswith(".cti"):
        print("This is an Cantera xml or cti file")
        #trims file
        solution_objects=create_trimmed_model(data_file, exclusion_list)

    if data_file.endswith(".inp") or data_file.endswith('.dat') \
                or data_file.endswith('.txt'):
        print("This is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(data_file, thermo_file, transport_file)

        #trims newly converted file
        solution_objects=create_trimmed_model(converted_file_name, \
                                    exclusion_list)
    else:
        print("File type not supported")

    write(data_file, solution_objects)

    return (solution_objects)



#gets arguments from terminal
parser=argparse.ArgumentParser(description='pyMARS main: \
            converts and trims mechanism files \n')
parser.add_argument('--file', help='input mechanism file name', type=str)
parser.add_argument('--species',  help='comma separated list input of species\
                                to exclude', type=str)
parser.add_argument('--thermo', help='thermodynamic data file', type=str)
parser.add_argument('--transport', help='transport data file', type=str)
args=parser.parse_args()

data_file=args.file

#if no species are to be trimmed
if args.species == None:
    exclusion_list=[]
else:
    exclusion_list=[str(item) for item in args.species.split(',')]
thermo_file=args.thermo
transport_file=args.transport

if __name__ == '__main__':
    readin(data_file, exclusion_list, thermo_file, transport_file)
