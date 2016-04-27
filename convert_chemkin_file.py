#converts Chemkin file to Cantera format

import os

def convert(mech_file, thermo_file, transport_file):
    input_dir='Input_Data_Files'
    output_dir= 'Output_Data_Files'

    #make file save path
    current_dir=os.path.dirname(os.path.abspath("~"))
    save_path=os.path.join(current_dir, output_dir)

    #make file save name
    converted_file_name=os.path.splitext(mech_file)[0]+'_converted'
    converted_file_path=os.path.join(save_path, converted_file_name)

    #convert and save file
    mech_file=os.path.join(input_dir, mech_file)
    thermo_file=os.path.join(input_dir, thermo_file)
    transport_file=os.path.join(input_dir, transport_file)
    input_line= "ck2cti -i=%s -t=%s -tr=%s --output=%s" %(mech_file, thermo_file, transport_file, converted_file_path)
    os.system(input_line)

    local_path=os.path.join(output_dir, converted_file_name)
    return local_path + '.cti'
