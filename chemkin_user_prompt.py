#prompts user for chemkin thermo and transport data files

def ask():
    thermo_name=str(raw_input('Enter thermo data file with extension: '))
    transport_name=str(raw_input('Enter transport data file with extension: '))
    return (thermo_name, transport_name)
