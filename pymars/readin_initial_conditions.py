import string

class Condition:
    """Class for storing initial conditions for autoignition simulations.
    """
    def __init__(self, pressure, temperature, equi, moles, name, species,
                 fuel, oxid, fuel_d, oxid_d
                 ):
        self.name = name
        self.pressure = pressure
        self.temperature = temperature
        self.moles = moles
        self.species = species
        self.fuel = fuel
        self.oxid = oxid
        self.equi = equi
        self.fuel_d = fuel_d
        self.oxid_d = oxid_d


def readin_conditions(initial_conditions_text_file):
    """Get list of initial conditions from text file.

    Parameters
    ----------
    initial_conditions_text_file : str
        Text file containing initial conditions for autoignition simulations

    Returns
    -------
    instance : list of Condition
        List of objects containing autoignition simulation initial conditions.
    
    Notes
    -----
    Input files are currently in a plaintext format, with the expected keywords and
    values in the format:
    CONV
    PRES 1.0
    TEMP 800
    EQUI 1.0
    FUEL nc7h16 1.0
    OXID o2 1.0
    OXID n2 3.76
    END

    Multiple initial conditions can be entered with multiple blocks, separated by an 
    empty line.

    """
    name_index = 0
    instance =[]
    with open(initial_conditions_text_file, 'r') as condition_file:
        lines = condition_file.readlines()

    for line in lines:
        name_index += 1
        if "CONV" in line:
            species_list = {}
            fuel_list = {}
            fuel = ""
            oxid_list = {}
            oxid = ""
            start = True
            reactants =''
            condition_name = 'Condition #' + str(name_index)
        if start is True:
            if "PRES" in line:
                p = str(line).strip(string.ascii_letters)
            if 'TEMP' in line:
                t = str(line).strip(string.ascii_letters)
            if 'EQUI' in line:
                e = str(line).strip(string.ascii_letters)

            if 'FUEL' in line:
                reactants += line.replace('FUEL ', '').replace(' ', ':').rstrip() + ','
                fuel += line.replace('FUEL ', '').replace(' ', ':').rstrip() + ','
                species = line.replace('FUEL ', '').rsplit(' ', 1)[0]
                moles = line.replace('FUEL ', '').rsplit(' ', 1)[1].rstrip()
                if species not in species_list.keys():
                    species_list[species] = moles

                if species not in fuel_list.keys():
                    fuel_list[species] = moles

            if 'OXID' in line:
                reactants += line.replace('OXID ', '').replace(' ', ':').rstrip() + ','
                oxid += line.replace('OXID ', '').replace(' ', ':').rstrip() + ','
                species = line.replace('OXID ', '').rsplit(' ', 1)[0]
                moles = line.replace('OXID ', '').rsplit(' ', 1)[1].rstrip()
                if species not in species_list.keys():
                    species_list[species] = moles
                
                if species not in oxid_list.keys():
                    oxid_list[species] = moles
            
            if 'END' in line:
                reactants = reactants[:-1]
                fuel = fuel[:-1]
                oxid = oxid[:-1]
                #normalize both 
                total = 0

                for species in fuel_list.keys():
                    total = total + float(fuel_list[species])
                
                for species in fuel_list.keys():
                    old = str(fuel_list[species])
                    #fuel_list[species] = float(fuel_list[species]) / total
                    #fuel = fuel.replace(old,str(fuel_list[species]))
                    #fuel_list[species] = float(fuel_list[species])
                    #fuel_list[species] = round(fuel_list[species],3)
                    
                total = 0

                for species in oxid_list.keys():
                    total = total + float(oxid_list[species])
                
                for species in oxid_list.keys():
                    old = str(oxid_list[species])
                    #oxid_list[species] = float(oxid_list[species]) / total
                    #oxid = oxid.replace(old,str(oxid_list[species]))
                    #oxid_list[species] = float(oxid_list[species])
                    #oxid_list[species] = round(oxid_list[species],3)


                instance.append(Condition(p, t, e, reactants, condition_name, species_list, fuel, oxid, fuel_list, oxid_list))
                start = False
    return instance
