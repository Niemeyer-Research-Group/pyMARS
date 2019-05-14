"""writes a solution object to a chemkin inp file

currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required
"""

import os
from string import Template

import cantera as ct
from cantera import ck2cti

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0


def eliminate(input_string, char_to_replace, spaces='single'):
    """Eliminate characters from a string

    Parameters
    ----------
    input_string : str
        string to be modified
    char_to_replace : list of str
        array of character strings to be removed
    
    Returns
    -------
    input_string : str
        string with characters eliminated
    
    """
    for char in char_to_replace:
        input_string = input_string.replace(char, "")
    if spaces == 'double':
        input_string = input_string.replace(" ", "  ")
    return input_string

def section_break(file, section_title):
    """Insert break and new section title into file

    Parameters
    ----------
    file : io.TextIOWrapper
        Open file handle for writing
    section_title : str
        title string for next section_break

    """
    file.write('!'+ "-"*75 + '\n')
    file.write('!  ' + section_title +'\n')
    file.write('!'+ "-"*75 + '\n')

def replace_multiple(input_string, replace_list):
    """Replace multiple characters in a string

    Parameters
    ----------
    input_string : str
        string to be modified
    replace_list : dict
        dict containing items to be replaced (value replaces key)
    
    Returns
    -------
    input_string : str 
        string with characters replaced

    """
    for original_character, new_character in replace_list.items():
        input_string = input_string.replace(original_character,
                                            new_character)
    return input_string


def build_arrhenius(reaction):
    """Builds Arrhenius coefficient string

    Parameters
    ----------
    reaction : cantera.Reaction
        cantera reaction object
    
    Returns
    -------
    list of str
        list of strings with Arrhenius coefficients

    """
    coeff_sum = sum(reaction.reactants.values())
    pre_exponential_factor = reaction.rate.pre_exponential_factor
    temperature_exponent = '{:.3f}'.format(reaction.rate.temperature_exponent)
    activation_energy = '{:.2f}'.format(reaction.rate.activation_energy / CALORIES_CONSTANT)

    if type(reaction) == ct.ElementaryReaction:
        if coeff_sum == 1:
            pre_exponential_factor = str('{:.3E}'.format(pre_exponential_factor))
        elif coeff_sum == 2:
            pre_exponential_factor = str('{:.3E}'.format(pre_exponential_factor * 1e3))
        elif coeff_sum == 3:
            pre_exponential_factor = str('{:.3E}'.format(pre_exponential_factor * 1e6))
    elif type(reaction) == ct.ThreeBodyReaction:
        if coeff_sum == 1:
            pre_exponential_factor = str('{:.3E}'.format(pre_exponential_factor * 1e3))
        elif coeff_sum == 2:
            pre_exponential_factor = str('{:.3E}'.format(pre_exponential_factor * 1e6))
    elif type(reaction) not in [ct.ElementaryReaction, ct.ThreeBodyReaction]:
        pre_exponential_factor = str('{:.3E}'.format(pre_exponential_factor))
    
    return [pre_exponential_factor, temperature_exponent, activation_energy]


def build_modified_arrhenius(reaction, pres_range):
    """Builds Arrhenius coefficient strings for high and low pressure ranges

    Parameters
    ----------
    reaction : cantera.Reaction
        cantera reaction object
    pres_range : str
        simple string ('high' or 'low') to designate pressure range
    
    Returns
    -------
    str
        Arrhenius coefficient string
        
    """
    if pres_range == 'high':
        pre_exponential_factor = reaction.high_rate.pre_exponential_factor
        temperature_exponent = '{:.3f}'.format(reaction.high_rate.temperature_exponent)
        activation_energy = '{:.2f}'.format(reaction.high_rate.activation_energy/CALORIES_CONSTANT)
        if len(reaction.products) == 1:
            pre_exponential_factor = str('{:.5E}'.format(pre_exponential_factor * 1e3))
        else:
            pre_exponential_factor = str('{:.5E}'.format(pre_exponential_factor))

        return [pre_exponential_factor, temperature_exponent, activation_energy]
    elif pres_range == 'low':

        pre_exponential_factor = reaction.low_rate.pre_exponential_factor
        temperature_exponent = '{:.3f}'.format(reaction.low_rate.temperature_exponent)
        activation_energy = '{:.2f}'.format(reaction.low_rate.activation_energy/CALORIES_CONSTANT)
        if len(reaction.products) == 1:
            pre_exponential_factor = str('{:.5E}'.format(pre_exponential_factor * 1e6))
        else:
            pre_exponential_factor = str('{:.5E}'.format(pre_exponential_factor * 1e3))

        return [pre_exponential_factor, temperature_exponent, activation_energy]
    else:
        raise ValueError('Pressure range needs to be high or low')


def build_nasa(nasa_coeffs, row):
    """Creates string of NASA polynomial coefficients

    Parameters
    ----------
    nasa_coeffs : numpy.ndarray
        Array of thermodynamic coefficients
    row : int
        which row to write coefficients in

    Returns
    -------
    line_coeffs : str
        Line of thermodynamic coefficients

    """
    line_coeffs = ''
    lines = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14]]
    line_index = lines[row-2]
    for ix, c in enumerate(nasa_coeffs):
        if ix in line_index:
            if c >= 0:
                line_coeffs += ' '
            line_coeffs += str('{:.8e}'.format(c))
    return line_coeffs


def build_species_string(trimmed_solution):
    """Formats species list at top of mechanism file
    
    Parameters
    ----------
    trimmed_solution : cantera.Solution
        Solution object for reduced model

    Returns
    -------
    species_list_string : str
        String with species list

    """
    species_list_string = ''
    line = 1
    for sp_string in trimmed_solution.species_names:
        sp = ' '
        #get length of string next species is added
        length_new = len(sp_string)
        length_string = len(species_list_string)
        total = length_new + length_string + 3
        #if string will go over width, wrap to new line
        if total >= 70*line:
            species_list_string += '\n'
            line += 1
        species_list_string += sp_string + ((16-len(sp_string))*sp)
    return species_list_string


def write(solution):
    """Function to write cantera solution object to inp file.

    Parameters
    ----------
    solution : cantera.Solution
        Model to be printed

    Returns
    -------
    str
        Name of trimmed model file (.inp)

    Examples
    --------
    >>> soln2ck.write(gas)

    """
    trimmed_solution = solution
    input_file_name_stripped = trimmed_solution.name
    cwd = os.getcwd()
    output_file_name = os.path.join(cwd, 'pym_' + input_file_name_stripped + '.inp')
    
    with open(output_file_name, 'w+') as the_file:

        #Write title block to file
        section_break(the_file, 'Chemkin File converted from Solution Object by pyMARS')

        #Write phase definition to file
        element_names = eliminate(str(trimmed_solution.element_names),
                                  ['[', ']', '\'', ',']
                                  )
        element_string = Template('ELEMENTS\n' + '$element_names\n' + 'END\n')
        the_file.write(element_string.substitute(element_names=element_names))
        species_names = build_species_string(trimmed_solution)
        species_string = Template('SPECIES\n' + '$species_names\n' + 'END\n')
        the_file.write(species_string.substitute(species_names=species_names))

        #Write species to file
        section_break(the_file, 'Species data')
        the_file.write('THERMO ALL' + '\n' + '   300.000  1000.000  5000.000' +'\n')
        phase_unknown_list = []

        # write data for each species in the Solution object
        for species in trimmed_solution.species():
            t_low = '{0:.3f}'.format(species.thermo.min_temp)
            t_max = '{0:.3f}'.format(species.thermo.max_temp)
            t_mid = '{0:.3f}'.format(species.thermo.coeffs[0])
            temp_range = str(t_low) + '  ' + str(t_max) + '  ' + t_mid
            species_comp = ''
            for atom in species.composition:
                species_comp += '{:<4}'.format(atom)
                species_comp += str(int(species.composition[atom]))
            if type(species.transport).__name__ == 'GasTransportData':
                species_phase = 'G'
            else:
                phase_unknown_list.append(species.name)
                species_phase = 'G'
            line_1 = (
                    '{:<18}'.format(species.name) +
                    '{:<6}'.format('    ') +
                    '{:<20}'.format(species_comp) +
                    '{:<4}'.format(species_phase) +
                    '{:<31}'.format(temp_range) +
                    '{:<1}'.format('1') +
                    '\n')
            the_file.write(line_1)
            line_2_coeffs = build_nasa(species.thermo.coeffs, 2)
            line_2 = line_2_coeffs  + '    2\n'
            the_file.write(line_2)
            line_3_coeffs = build_nasa(species.thermo.coeffs, 3)
            line_3 = line_3_coeffs + '    3\n'
            the_file.write(line_3)
            line_4_coeffs = build_nasa(species.thermo.coeffs, 4)
            line_4 = line_4_coeffs + '                   4\n'
            the_file.write(line_4)

        the_file.write('END\n')

        #Write reactions to file
        section_break(the_file, 'Reaction Data')
        the_file.write('REACTIONS\n')

        #write data for each reaction in the Solution Object
        for reac in trimmed_solution.reactions():
            equation_string = eliminate(reac.equation, ' ', 'single')
            
            if type(reac) == ct.ThreeBodyReaction:
                arrhenius = build_arrhenius(reac)
                main_line = (
                            '{:<51}'.format(equation_string) +
                            '{:>9}'.format(arrhenius[0]) +
                            '{:>9}'.format(arrhenius[1]) +
                            '{:>11}'.format(arrhenius[2]) +
                            '\n'
                            )
                the_file.write(main_line)
                #trimms efficiencies list
                efficiencies = reac.efficiencies
                trimmed_efficiencies = reac.efficiencies
                for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]
                replace_list_2 = {'{':'',
                                  '}':'/',
                                  '\'':'',
                                  ':':'/',
                                  ',':'/'
                                  }
                efficiencies_string = replace_multiple(
                                            str(trimmed_efficiencies),
                                            replace_list_2
                                            )
                secondary_line = str(efficiencies_string) + '\n'
                if bool(efficiencies) is True:
                    the_file.write(secondary_line)
            elif type(reac) == ct.ElementaryReaction:
                arrhenius = build_arrhenius(reac)
                main_line = ('{:<51}'.format(equation_string) +
                             '{:>9}'.format(arrhenius[0]) +
                             '{:>9}'.format(arrhenius[1]) +
                             '{:>11}'.format(arrhenius[2]) +
                             '\n'
                             )
                the_file.write(main_line)
            elif type(reac) == ct.FalloffReaction:
                arr_high = build_modified_arrhenius(reac, 'high')
                main_line = ('{:<51}'.format(equation_string) +
                             '{:>9}'.format(arr_high[0]) +
                             '{:>9}'.format(arr_high[1]) +
                             '{:>11}'.format(arr_high[2]) +
                             '\n'
                             )
                the_file.write(main_line)
                arr_low = build_modified_arrhenius(reac, 'low')
                second_line = ('     LOW  /' +
                               '  ' + arr_low[0] +
                               '  ' + arr_low[1] +
                               '  ' + arr_low[2] + '/\n'
                               )
                the_file.write(second_line)

                #If optional Arrhenius data included:
                try:
                    third_line = ('     TROE/' +
                                  '   ' + str(reac.falloff.parameters[0]) +
                                  '  ' + str(reac.falloff.parameters[1]) +
                                  '  ' + str(reac.falloff.parameters[2]) +
                                  '  ' + str(reac.falloff.parameters[3]) +' /\n'
                                  )
                    the_file.write(third_line)
                except IndexError:
                    pass
                #trims efficiencies list
                efficiencies = reac.efficiencies
                trimmed_efficiencies = reac.efficiencies
                for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]
                replace_list_2 = {'{':'',
                                  '}':'/',
                                  '\'':'',
                                  ':':'/',
                                  ',':'/'
                                  }
                efficiencies_string = replace_multiple(str(trimmed_efficiencies), replace_list_2)

                fourth_line = str(efficiencies_string) + '\n'
                if bool(efficiencies) is True:
                    the_file.write(fourth_line)
            #duplicate option
            if reac.duplicate:
                duplicate_line = ' DUPLICATE' +'\n'
                the_file.write(duplicate_line)

        the_file.write('END')

    #Test mechanism file

    original_solution = solution
    #convert written chemkin file to cti, and get solution
    parser = ck2cti.Parser()
    outName = 'test_file.cti'
    parser.convertMech(output_file_name, outName=outName)
    new_solution = ct.Solution(outName)

    #TODO: test new solution vs original solution
    #test(original_solution, new_solution)
    os.remove(outName)
    return output_file_name


    """
    def build_falloff(j):

        Creates falloff reaction Troe parameter string

        param j:
            Cantera falloff parameters object

        falloff_str = str(',\n        falloff = Troe(' +
                        'A = ' + str(j[0]) +
                        ', T3 = ' + str(j[1]) +
                        ', T1 = ' + str(j[2]) +
                        ', T2 = ' + str(j[3]) +')       )\n\n')
        return falloff_str
    """
