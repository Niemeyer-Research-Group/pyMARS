"""Converts Chemkin file to Cantera format."""
import os

def convert(mech_file, thermo_file=None, transport_file=None):
    """Function to convert chemkin mechanism files to CTI.

    Parameters
    ----------
    mech_file : str
        Chemkin mechanism file
    thermo_file : str
        Chemkin thermodynamic properties file
    transport_file : str
        Chemkin transport data file
    Returns
    -------
    str
        Converted Cantera mechanism file

    Example
    -------
    >>> convert('gri30.inp')

    """

    # make file save path
    current_dir = os.path.dirname(os.path.abspath("~"))
    save_path = current_dir

    # make file save name
    converted_file_name = os.path.splitext(
        mech_file)[0] + '_converted' + '.cti'
    converted_file_path = os.path.join(save_path, converted_file_name)

    mech_file = os.path.join(current_dir, mech_file)
    try:
        thermo_file_path = os.path.join(current_dir, thermo_file)
    except:
        pass
    try:
        transport_file_path = os.path.join(current_dir, transport_file)
    except:
        pass

    # calls ck2cti based on given files
    if thermo_file is None:
        if transport_file == None:
            input_line = "ck2cti --input=%s --output=%s" \
                % (mech_file, converted_file_path)
            print('no thermo, no transport')
        else:
            input_line = "ck2cti --input=%s --transport=%s --output=%s" \
                % (mech_file, transport_file_path, converted_file_path)
            print('no thermo, yes transport')

    if thermo_file is not None:
        if transport_file is None:
            input_line = "ck2cti --input=%s --thermo=%s  --output=%s" \
                % (mech_file, thermo_file_path, converted_file_path)
            print('yes thermo, no transport')
        else:
            input_line = "ck2cti --input=%s --thermo=%s --transport=%s --output=%s" \
                % (mech_file, thermo_file_path,
                   transport_file_path, converted_file_path)
            print('yes thermo, yes transport')
    print(input_line + '\n')

    # convert and save file
    os.system(input_line)
    local_path = os.path.join(converted_file_name)
    return local_path
