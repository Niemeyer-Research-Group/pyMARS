#converts Chemkin file to Cantera format


import os

def convert(data_file):
    output_dir= 'Output_Data_Files'
    current_dir=os.path.dirname(os.path.abspath("~"))
    save_path=os.path.join(current_dir, output_dir)

    converted_file_name=os.path.splitext(data_file)[0]+'_converted'
    converted_file_path=os.path.join(save_path, converted_file_name)
    input_line= "ck2cti --input=%s --output=%s" %(data_file, converted_file_path)
    os.system(input_line)
    return converted_file_name + '.cti'
