"""writes a solution object to a chemkin inp file

currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required
"""

import os

import cantera as ct

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0


def section_break(section_title):
    """Return string with break and new section title

    Parameters
    ----------
    section_title : str
        title string for next section_break
    
    Returns
    -------
    str
        String with section title and breaks

    """
    return('!'+ '-' * 75 + '\n' +
           f'!  {section_title}\n' +
           '!'+ '-' * 75 + '\n'
           )


def build_arrhenius(rate, reaction_order, reaction_type):
    """Builds Arrhenius coefficient string based on reaction type.

    Parameters
    ----------
    rate : cantera.Arrhenius
        Arrhenius-form reaction rate coefficient
    reaction_order : int or float
        Order of reaction (sum of reactant stoichiometric coefficients)
    reaction_type : {cantera.ElementaryReaction, cantera.ThreeBodyReaction, cantera.PlogReaction}
        Type of reaction

    Returns
    -------
    str
        String with Arrhenius coefficients

    """
    if reaction_type in [ct.ElementaryReaction, ct.PlogReaction]:
        pre_exponential_factor = rate.pre_exponential_factor * 1e3**(reaction_order - 1)

    elif reaction_type == ct.ThreeBodyReaction:
        pre_exponential_factor = rate.pre_exponential_factor * 1e3**reaction_order

    elif reaction_type in [ct.FalloffReaction, ct.ChemicallyActivatedReaction]:
        raise ValueError('Function does not support falloff or chemically activated reactions')
    else:
        raise NotImplementedError('Reaction type not supported: ', reaction_type)
    
    arrhenius = [f'{pre_exponential_factor:.4e}', 
                 f'{rate.temperature_exponent:.3f}', 
                 f'{(rate.activation_energy / CALORIES_CONSTANT):.2f}'
                 ]
    return '  '.join(arrhenius)


def build_falloff_arrhenius(rate, reaction_order, reaction_type, pressure_limit):
    """Builds Arrhenius coefficient strings for falloff and chemically-activated reactions.

    Parameters
    ----------
    rate : cantera.Arrhenius
        Arrhenius-form reaction rate coefficient
    reaction_order : int or float
        Order of reaction (sum of reactant stoichiometric coefficients)
    reaction_type : {ct.FalloffReaction, ct.ChemicallyActivatedReaction}
        Type of reaction
    pressure_limit : {'high', 'low'}
        string designating pressure limit
    
    Returns
    -------
    str
        Arrhenius coefficient string

    """
    assert pressure_limit in ['low', 'high'], 'Pressure range needs to be high or low'

    # Each needs more complicated handling due if high- or low-pressure limit
    if reaction_type == ct.FalloffReaction:
        if pressure_limit == 'low':
            pre_exponential_factor = rate.pre_exponential_factor * 1e3**(reaction_order)
        elif pressure_limit == 'high':
            pre_exponential_factor = rate.pre_exponential_factor * 1e3**(reaction_order - 1)

    elif reaction_type == ct.ChemicallyActivatedReaction:
        if pressure_limit == 'low':
            pre_exponential_factor = rate.pre_exponential_factor * 1e3**(reaction_order - 1)
        elif pressure_limit == 'high':
            pre_exponential_factor = rate.pre_exponential_factor * 1e3**(reaction_order - 2)
    else:
        raise ValueError('Reaction type not supported: ', reaction_type)

    arrhenius = [f'{pre_exponential_factor:.4e}', 
                 f'{rate.temperature_exponent:.3f}', 
                 f'{(rate.activation_energy / CALORIES_CONSTANT):.:.3e}'
                 ]
    return '  '.join(arrhenius)


def build_falloff(parameters, falloff_function):
    """Creates falloff reaction Troe parameter string

    Parameters
    ----------
    parameters : numpy.ndarray
        Array of falloff parameters; length varies based on ``falloff_function``
    falloff_function : {'Troe', 'SRI'}
        Type of falloff function

    Returns
    -------
    falloff_string : str
        String of falloff parameters

    """
    if falloff_function == ct.TroeFalloff:
        falloff_string = ('TROE / ' +
                          f'{parameters[0]}  {parameters[1]}  '
                          f'{parameters[2]}  {parameters[3]} /\n'
                          )
    elif falloff_function == ct.SriFalloff:
        falloff_string = ('SRI / ' + 
                          f'{parameters[0]}  {parameters[1]}  ' +
                          f'{parameters[2]}  {parameters[3]}  {parameters[4]} /\n'
                          )
    else:
        raise NotImplementedError(f'Falloff function not supported: {falloff_function}')

    return falloff_string


def build_thermo_data(nasa_coeffs, row):
    """Creates string of NASA polynomial coefficients for thermodynamic data.

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


def write(solution):
    """Function to write cantera solution object to inp file.

    Parameters
    ----------
    solution : cantera.Solution
        Model to be written

    Returns
    -------
    output_file_name : str
        Name of output model file (.inp)

    Examples
    --------
    >>> gas = cantera.Solution('gri30.cti')
    >>> soln2ck.write(gas)
    reduced_gri30.inp

    """
    # Remove extension from filename
    input_file_name = os.path.splitext(os.path.basename(solution.name))[0]
    output_file_name = f'reduced_{input_file_name}.inp'
    
    with open(output_file_name, 'w') as the_file:

        # Write title block to file
        the_file.write(section_break('Chemkin file converted from solution object'))

        # write species and element lists to file
        element_names = '  '.join(solution.element_names)
        the_file.write(
            'ELEMENTS\n' + 
            f'{element_names}\n' +
            'END\n\n'
            )
        species_names = fill('  '.join(solution.species_names), width=60)
        the_file.write(
            'SPECIES\n' + 
            f'{species_names}\n'
            'END\n\n'
            )

        # Write reactions to file
        the_file.write('REACTIONS\n')

        #write data for each reaction in the Solution Object
        for reaction in solution.reactions():

            reaction_string = f'{reaction.equation:<51}'

            # The Arrhenius parameters that follow the equation string on the main line 
            # depend on the type of reaction.
            if type(reaction) in [ct.ElementaryReaction, ct.ThreeBodyReaction]:
                arrhenius = build_arrhenius(
                    reaction.rate, 
                    sum(reaction.reactants.values()), 
                    type(reaction)
                    )

            elif type(reaction) == ct.FalloffReaction:
                # high-pressure limit is included on the main reaction line
                arrhenius = build_falloff_arrhenius(
                    reaction.high_rate, 
                    sum(reaction.reactants.values()), 
                    ct.FalloffReaction,
                    'high'
                    )

            elif type(reaction) == ct.ChemicallyActivatedReaction:
                # low-pressure limit is included on the main reaction line
                arrhenius = build_falloff_arrhenius(
                    reaction.low_rate, 
                    sum(reaction.reactants.values()), 
                    ct.ChemicallyActivatedReaction,
                    'low'
                    )

            elif type(reaction) == ct.ChebyshevReaction:
                arrhenius = '1.0e0  0.0  0.0'

            elif type(reaction) == ct.PlogReaction:
                arrhenius = build_arrhenius(
                    reaction.rates[0][1],
                    sum(reaction.reactants.values()), 
                    ct.PlogReaction
                    )
            reaction_string += arrhenius + '\n'
            
            # need to trim and print third-body efficiencies, if present
            if type(reaction) in [ct.ThreeBodyReaction, ct.FalloffReaction, 
                                  ct.ChemicallyActivatedReaction
                                  ]:
                # trims efficiencies list
                reduced_efficiencies = {s:reaction.efficiencies[s] 
                                        for s in reaction.efficiencies
                                        if s in solution.species_names
                                        }
                efficiencies_str = '  '.join([f'{s}/{v}/' for s, v in reduced_efficiencies.items()])
                if efficiencies_str:
                    reaction_string += efficiencies_str + '\n'
            
            # 

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

    # TODO
    # write thermo data
    
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

    return output_file_name
