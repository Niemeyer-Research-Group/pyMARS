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
            def __init__(self, pressure, temperature, moles, name):
                self.name = name
                self.pressure = pressure
                self.temperature = temperature
                self.moles = moles

        name_index = 0
        instance =[]
        for line in condition_file:
                    name_index += 1
                    if "CONV" in line:
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

                        if 'END' in line:
                            reactants = reactants[:-1]
                            instance.append(condition(p, t, reactants, condition_name))
                            start = False

    return instance
