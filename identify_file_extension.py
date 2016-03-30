

import cantera as ct
import os

def readin(data_file):
    """Function to import data file and identify format.

    Parameters
    ----------
    data_file:
        Local Chemkin or Cantera data file containing mechanism information

    Returns
    -------
        Calls function corresponding to data file
    """

    if data_file.endswith(".xml") or data_file.endswith(".cti"):
        print("This is an Cantera xml or cti file")
        from create_trimmed_model import create_trimmed_model
        create_trimmed_model(data_file, exclusion_list)
    elif data_file.endswith(".inp"):
        print("This is a Chemkin inp file")
        #convert file to cti
        converted_file_name = os.path.splitext(data_file)[0] + '_converted'
        input_line= "ck2cti --input=%s --output=%s" %(data_file, converted_file_name)
        os.system(input_line)
        #trims newly converted file
        from create_trimmed_model import create_trimmed_model
        create_trimmed_model(converted_file_name, exclusion_list)
    else:
        print("File type not supported")

readin('gri30.txt')
