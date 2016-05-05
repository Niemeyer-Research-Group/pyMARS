#prompts user for chemkin thermo and transport data files

def ask():
    q1=str(raw_input('Separate thermo file? [y/n]: '))
    q2=str(raw_input('Separate transport file? [y/n]: '))
    if q1 == 'y':
        thermo_name=str(raw_input('Enter thermo data file with extension: '))
    else:
        thermo_name='none'
    if q2 == 'y':
        transport_name=str(raw_input('Enter transport data file with extension: '))
    else: transport_name='none'
    return (thermo_name, transport_name)
