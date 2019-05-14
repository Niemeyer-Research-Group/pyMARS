"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff, Plog and ThreeBody Reactions
Cantera development version 2.3.0a2 required
"""

import os
import textwrap
from string import Template

import cantera as ct

from .name_trim import name_trim

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

#1 debye = d coulomb-meters
DEBEYE_CONVERSION = 3.33564e-30

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


def wrap_nasa(input_string):
    """Wrap string to cti NASA format width

    Parameters
    ----------
    output_string : str
        string to be wrapped
    
    Returns
    -------
    str
        String wrapped to appropriate width

    """
    return textwrap.fill(input_string, width=50, subsequent_indent=16*' ')


def section_break(file, section_title):
    """Insert break and new section title into cti file

    Parameters
    ----------
    file : io.TextIOWrapper
        Open file handle for writing
    section_title : str
        title string for next section_break

    """
    file.write('#' + '-' * 75 + '\n')
    file.write('#  ' + section_title + '\n')
    file.write('#' + '-' * 75 + '\n\n')


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
    str
        String with Arrhenius coefficients

    """
    if type(reaction) == ct.PlogReaction:
        coeff_sum = reaction[0]
        equation_object = reaction[1][1]
        if coeff_sum == 1:
            pre_exponential_factor = equation_object.pre_exponential_factor
        if coeff_sum == 2:
            pre_exponential_factor = equation_object.pre_exponential_factor * 1e3
        if coeff_sum == 3:
            pre_exponential_factor = equation_object.pre_exponential_factor * 1e6
        temperature_exponent = equation_object.temperature_exponent
        activation_energy = equation_object.activation_energy / CALORIES_CONSTANT
    else:
        coeff_sum = sum(reaction.reactants.values())
        pre_exponential_factor = reaction.rate.pre_exponential_factor
        temperature_exponent = str(reaction.rate.temperature_exponent)
        activation_energy = str(reaction.rate.activation_energy / CALORIES_CONSTANT)

    if type(reaction) == ct.ElementaryReaction:
        if coeff_sum == 1:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor))
        if coeff_sum == 2:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor * 1e3))
        if coeff_sum == 3:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor * 1e6))
    elif type(reaction) == ct.ThreeBodyReaction:
        if coeff_sum == 1:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor * 1e3))
        if coeff_sum == 2:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor * 1e6))
    else:
        pre_exponential_factor = str('{:.5E}'.format(pre_exponential_factor))

    arrhenius = [pre_exponential_factor, temperature_exponent, activation_energy]
    return '[' + ', '.join(arrhenius) + ']'


def build_modified_arrhenius(equation_object, pres_range):
    """Builds Arrhenius coefficient strings for high and low pressure ranges

    Parameters
    ----------
    equation_object : cantera.Reaction
        cantera reaction object
    pres_range : str
        simple string ('high' or 'low') to designate pressure range
    
    Returns
    -------
    str
        Arrhenius coefficient string

    """

    if pres_range == 'high':
        pre_exponential_factor = equation_object.high_rate.pre_exponential_factor
        temperature_exponent = equation_object.high_rate.temperature_exponent
        activation_energy = (equation_object.high_rate.activation_energy /
                            CALORIES_CONSTANT)
        if len(equation_object.products) == 1:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor * 1e3))
        else:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor))
        arrhenius_high = [pre_exponential_factor,
                            temperature_exponent,
                            activation_energy
                            ]
        return str(arrhenius_high).replace("\'", "")
    elif pres_range == 'low':
        pre_exponential_factor = equation_object.low_rate.pre_exponential_factor
        temperature_exponent = equation_object.low_rate.temperature_exponent
        activation_energy = (equation_object.low_rate.activation_energy /
                            CALORIES_CONSTANT)

        if len(equation_object.products) == 1:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor * 1e6))
        else:
            pre_exponential_factor = str(
                            '{:.5E}'.format(pre_exponential_factor * 1e3))
        arrhenius_low = [pre_exponential_factor,
                        temperature_exponent,
                        activation_energy
                        ]
        return str(arrhenius_low).replace("\'", "")
    else:
        raise ValueError('Pressure range needs to be high or low')


def build_falloff(falloff_params):
    """Creates falloff reaction Troe parameter string

    Parameters
    ----------
    falloff_params : numpy.ndarray
        Array of Troe falloff parameters 

    Returns
    -------
    falloff_string : str
        String of Troe falloff parameters
    """
    falloff_string = str(
                ',\n        falloff = Troe(' +
                'A = ' + str(falloff_params[0]) +
                ', T3 = ' + str(falloff_params[1]) +
                ', T1 = ' + str(falloff_params[2]) +
                ', T2 = ' + str(falloff_params[3]) + ')       )\n\n'
                )
    return falloff_string


def build_species_string(species_names):
    """Formats species list at top of file
    
    Parameters
    ----------
    species_names : list of str
        list of species names

    Returns
    -------
    species_list_string : str
        String with formatted species list

    """
    species_list_string = ''
    line = 1
    for sp_str in species_names:
        #get length of string next species is added
        length_new = len(sp_str)
        length_string = len(species_list_string)
        total = length_new +length_string +3
        #if string will go over width, wrap to new line
        if line == 1:
            if total >= 55:
                species_list_string += '\n'
                species_list_string += ' ' * 17
                line += 1
        if line > 1:
            if total >= 70 * line:
                species_list_string += '\n'
                species_list_string += ' ' * 17
                line += 1
        species_list_string += sp_str + ' '
    return species_list_string


def write(solution):
    """Function to write cantera solution object to cti file.

    Parameters
    ----------
    solution : cantera.Solution
        Model to be written

    Returns
    -------
    output_file_name : str
        Name of output model file (.cti)

    Examples
    --------
    >>> soln2cti.write(gas)

    """
    # Remove extension from filename
    input_file_name = name_trim(solution.name)
    output_file_name = os.path.join(os.getcwd(), 'pym_' + input_file_name + '.cti')

    with open(output_file_name, 'w') as the_file: 

        # Write title block to file
        section_break(the_file, 'CTI File converted from solution object')
        unit_string = ('units(length = "cm", time = "s",' +
                       ' quantity = "mol", act_energy = "cal/mol")'
                       )
        the_file.write(unit_string + '\n\n')

        # Write Phase definition to file
        element_names = ' '.join(solution.element_names)        
        species_names = build_species_string(solution.species_names)
        the_file.write(
            f'ideal_gas(name = "{input_file_name}", \n' +
            f'     elements = "{element_names}", \n' +
            f'     species = """ {species_names} """, \n' +
            f'     reactions = "all", \n' +
            f'     initial_state = state(temperature = {solution.T}, ' +
            f'pressure = {solution.P})   )\n\n'
            )

        # Write species data to file
        section_break(the_file, 'Species data')
        
        for species in solution.species():
            # build strings with low- and high-temperature 7 NASA coefficients
            nasa_range_low = f'[{species.thermo.min_temp}, {species.thermo.coeffs[0]}]'
            nasa_coeffs_low = ["{:.9e}".format(c) for c in species.thermo.coeffs[8:15]]
            nasa_coeffs_low = wrap_nasa('[' + ', '.join(nasa_coeffs_low) + ']')
            
            nasa_range_high = f'[{species.thermo.coeffs[0]}, {species.thermo.max_temp}]'
            nasa_coeffs_high = ["{:.9e}".format(c) for c in species.thermo.coeffs[1:8]]
            nasa_coeffs_high = wrap_nasa('[' + ', '.join(nasa_coeffs_high) + ']')

            composition = ', '.join([f'{s}:{int(v)}' for s, v in species.composition.items()])

            # start writing composition and thermo data
            the_file.write(
                f'species(name = "{species.name}",\n' +
                f'    atoms = {composition}, \n' +
                f'    thermo = (\n' +
                f'       NASA(   {nasa_range_low}, {nasa_coeffs_low}  ),\n' +
                f'       NASA(   {nasa_range_high}, {nasa_coeffs_high}  )\n' +
                f'               ),\n'
                )

            #check if species has defined transport data, and write that if so
            if species.transport:
                indent = '                   '
                the_file.write(
                    f'    transport = gas_transport(\n' +
                    indent + f'geom = "{species.transport.geometry}",\n' +
                    indent + f'diam = {species.transport.diameter * 1e10}, \n' +
                    indent + f'well_depth = {species.transport.well_depth / ct.boltzmann}, \n' +
                    indent + f'polar = {species.transport.polarizability * 1e30}, \n' +
                    indent + f'rot_relax = {species.transport.rotational_relaxation}, \n'
                    )
                if species.transport.dipole != 0:
                    dipole = species.transport.dipole / DEBEYE_CONVERSION
                    the_file.write(indent + f'dipole= {dipole}) \n')
            
            the_file.write('        )\n\n')

        # Write reactions to file
        section_break(the_file, 'Reaction Data')

        # write data for each reaction
        for idx, reaction in enumerate(solution.reactions()):

            reaction_string = f'#  Reaction {idx + 1}\n'

            if type(reaction) == ct.ElementaryReaction:
                arrhenius = build_arrhenius(reaction)
                reaction_string += f'reaction( "{reaction.equation}",  {arrhenius}'

            elif type(reaction) == ct.ThreeBodyReaction:
                arrhenius = build_arrhenius(reaction)
                reaction_string += f'three_body_reaction( "{reaction.equation}",  {arrhenius}'

                # trims efficiencies list
                reduced_efficiencies = {s:reaction.efficiencies[s] 
                                        for s in reaction.efficiencies 
                                        if s in solution.species_names
                                        }
                efficiencies_str = '  '.join([f'{s}:{v}' for s, v in reduced_efficiencies.items()])
                if efficiencies_str:
                    reaction_string += f',\n         efficiencies = " {efficiencies_str} "'
                                    
            elif type(reaction) == ct.FalloffReaction:
                arrhenius_high = build_modified_arrhenius(reaction, 'high')
                arrhenius_low = build_modified_arrhenius(reaction, 'low')

                reaction_string += (f'falloff_reaction( "{reaction.equation}",\n' +
                                    f'         kf = {arrhenius_high},\n' +
                                    f'         kf0 = {arrhenius_low}'
                                    )
                
                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters)
                    reaction_string += (',\n' + 
                                        '         falloff = ' + falloff_str
                                        )
                
                # trims efficiencies list
                reduced_efficiencies = {s:reaction.efficiencies[s] 
                                        for s in reaction.efficiencies 
                                        if s in solution.species_names
                                        }
                efficiencies_str = '  '.join([f'{s}:{v}' for s, v in reduced_efficiencies.items()])
                if efficiencies_str:
                    reaction_string += f',\n         efficiencies = " {efficiencies_str} "'

            elif type(reaction) == ct.PlogReaction:
                reaction_string += f'pdep_arrhenius( "{reaction.equation}",\n'

                
                sum_coeffs = sum(equation_object.reactants.values())

                for rate in reaction.rates:
                    pressure = str('{:f}'.format(rate_line[0]/ct.one_atm))
                    arrhenius = build_arrhenius((sum_coeffs,rate_line))
                    arrhenius = arrhenius[1:-1]
                    reaction_string = Template(
                            reaction_string + 
                            "               [($pressure, 'atm'), $Arr],\n"
                            )
                    reaction_string = reaction_string.substitute(
                            pressure=pressure,
                            Arr = arrhenius
                            )
                if equation_object.duplicate is True:
                    reaction_string = (reaction_string +
                                    '               options=\'duplicate\')\n\n'
                                    )
                else:
                    reaction_string = reaction_string[:-1]
                    reaction_string = reaction_string + ')\n\n'
                the_file.write(reaction_string)
            else:
                raise NotImplementedError(f'Unsupported reaction type: {type(reaction)}')
            
            if reaction.duplicate:
                reaction_string += ',\n        options = "duplicate"'
                                
            reaction_string += ')\n\n'
            the_file.write(reaction_string)
    
    return output_file_name
