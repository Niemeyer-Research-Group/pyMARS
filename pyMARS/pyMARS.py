
import os, sys, argparse


from create_trimmed_model import create_trimmed_model
from convert_chemkin_file import convert
from write_to_cti import write
from autoignition_module import run_sim


def readin(args):
    """Function to import data file and identify format.

    Arguments
        --file: Input mechanism file (ex. --file=gri30.cti)
        --species: Species to eliminate (ex. --species='H, OH')
        --thermo: Thermo data file if Chemkin format (ex. --thermo= thermo.dat)
        --transport: Transport data file if Chemkin format
        --plot
        --points
        --writecsv
        --writehdf5

    ----------
    Returns
        Converted mechanism file
        Trimmed Solution Object
        Trimmed Mechanism file
    """

    data_file= args.file
    ext= os.path.splitext(data_file)[1]
    thermo_file = args.thermo
    transport_file=args.transport


    #if no species are to be trimmed
    if args.species == None:
        exclusion_list=[]
    else:
        exclusion_list=[str(item) for item in args.species.split(',')]
        print(exclusion_list)


    if ext == ".cti" or ext == ".xml":

        print("\nThis is an Cantera xml or cti file\n")
        #trims file
        solution_objects=create_trimmed_model(data_file, exclusion_list)
        write(solution_objects[1])
        run_sim(data_file, args)

    elif ext == ".inp" or ext == ".dat" or ext == ".txt":

        print("This is a Chemkin file")
        #convert file to cti
        converted_file_name = convert(data_file, thermo_file, transport_file)

        #trims newly converted file
        solution_objects=create_trimmed_model(converted_file_name, \
                                    exclusion_list)
        print(converted_file_name)
        run_sim(converted_file_name, args)


    else:
        print("File type not supported")
