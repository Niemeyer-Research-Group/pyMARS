#converts Chemkin file to Cantera format

import os
from chemkin_user_prompt import ask

def convert(mech_file):
    input_dir='Data_Files'
    output_dir= 'Output_Data_Files'

    #make file save path
    current_dir=os.path.dirname(os.path.abspath("~"))
    save_path=os.path.join(current_dir, output_dir)

    #make file save name
    converted_file_name=os.path.splitext(mech_file)[0]+'_converted'
    converted_file_path=os.path.join(save_path, converted_file_name)

    #determine if thermo and transport files are Separate
    external_files=ask()
    thermo_file=external_files[0]
    transport_file=external_files[1]

    mech_file=os.path.join(input_dir, mech_file)
    thermo_file_path=os.path.join(input_dir, thermo_file)
    transport_file_path=os.path.join(input_dir, transport_file)

    if thermo_file == 'none':
        if transport_file == 'none':
            input_line= "ck2cti --input=%s --output=%s" \
                %(mech_file, converted_file_path)
            print('no thermo, no transport')
        else:
            input_line= "ck2cti --input=%s --transport=%s --output=%s" \
                %(mech_file, transport_file_path, converted_file_path)
            print('no thermo, yes transport')

    if thermo_file != 'none':
        if transport_file == 'none':
            input_line= "ck2cti --input=%s --thermo=%s  --output=%s" \
                    %(mech_file, thermo_file_path, converted_file_path)
            print('yes thermo, no transport')
        else:
            input_line= "ck2cti --input=%s --thermo=%s --transport=%s --output=%s" \
                    %(mech_file, thermo_file_path, \
                    transport_file_path, converted_file_path)
            print('yes thermo, yes transport')
    print(input_line)


    #convert and save file
    os.system(input_line)
    local_path=os.path.join(output_dir, converted_file_name)
    return local_path + '.cti'
