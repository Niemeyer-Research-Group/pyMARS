import string

def readin_conditions(initial_conditions_text_file):
    """
        Get initial conditions from text file
        :param initial_conditoins_text_file:
            Text file containing conditions in the format
            CONV
            PRES 1.0
            TEMP 600
            REAC nc7h16 1.0
            REAC o2 11.0
            REAC n2 41.36
            END
        :return instance:
            Array of objects representing initial state
        """
    with open(initial_conditions_text_file, 'r') as condition_file:

        class condition:
            def __init__(self, pressure, temperature, moles, name, species):
                self.name = name
                self.pressure = pressure
                self.temperature = temperature
                self.moles = moles
                self.species = species

        name_index = 0
        instance =[]

        for line in condition_file:
                    name_index += 1
                    if "CONV" in line:
                        species_list = {}
                        start = True
                        reactants =''
                        condition_name = 'Condition #' + str(name_index)
                    if start is True:
                        if "PRES" in line:
                            p = str(line).translate(None, string.letters)
                        if 'TEMP' in line:
                            t = str(line).translate(None, string.letters)
                        if 'REAC' in line:
                            reactants += line.replace('REAC ', '').replace(' ', ':').rstrip() + ','
                            species = line.replace('REAC ', '').rsplit(' ', 1)[0]
                            moles = line.replace('REAC ', '').rsplit(' ', 1)[1].rstrip()
                            if species not in species_list.keys():
                                species_list[species] = moles
                        if 'END' in line:
                            reactants = reactants[:-1]
                            instance.append(condition(p, t, reactants, condition_name, species_list))
                            start = False
    return instance

readin_conditions('input.txt')
