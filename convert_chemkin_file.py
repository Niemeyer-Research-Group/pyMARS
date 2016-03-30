#converts Chemkin file to Cantera format


import os

def convert(data_file):
    converted_file_name=os.path.splitext(data_file)[0]+'_converted'
    input_line= "ck2cti --input=%s --output=%s" %(data_file, converted_file_name)
    os.system(input_line)
    return converted_file_name + '.cti'
