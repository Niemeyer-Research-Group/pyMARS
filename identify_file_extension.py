#trims reaction mechanism files

import cantera as ct
import os

def readin(data_file, exclusion_list):
    """Function to import data file and identify format.

    Parameters
    ----------
    data_file:
        Local Chemkin or Cantera data file containing mechanism information

    Returns
    -------
        Converted mechanism file
        Solution object index with species and reactions removed
    """
    
    global Result
    #import working functions
    from create_trimmed_model import create_trimmed_model
    from convert_chemkin_file import convert

    if data_file.endswith(".xml") or data_file.endswith(".cti"):
        print("This is an Cantera xml or cti file")

        #trims file
        Result=create_trimmed_model(data_file, exclusion_list)

    elif data_file.endswith(".inp"):
        print("This is a Chemkin inp file")

        #convert file to cti
        converted_file_name = convert(data_file)

        #trims newly converted file
        Result=create_trimmed_model(converted_file_name, exclusion_list)

    else:
        print("File type not supported")

    return (Result)
