"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff, Plog and ThreeBody Reactions
Cantera development version 2.3.0a2 required
"""

import os
import math
from textwrap import fill

import cantera as ct

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
DEBEYE_CONVERSION = 3.33564e-30

indent = ['',
          ' ',
          '  ',
          '   ',
          '    ',
          '     ',
          '      ',
          '       ',
          '        ',
          '          ',
          '           ',
          '            ',
          '             ',
          '              ',
          '               ',
          '                '
          ]


def section_break(section_title):
    """Return string with break and new section title

    Parameters
    ----------
    section_title : str
        title string for next section break
    
    Returns
    -------
    str
        String with section title and breaks

    """
    return('#' + '-' * 75 + '\n' + 
           f'#  {section_title}\n' +
           '#' + '-' * 75 + '\n\n'
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
    
    arrhenius = [f'{pre_exponential_factor:.6e}', 
                 str(rate.temperature_exponent), 
                 str(rate.activation_energy / CALORIES_CONSTANT)
                 ]
    return ', '.join(arrhenius)


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

    arrhenius = [f'{pre_exponential_factor:.6E}', 
                 str(rate.temperature_exponent), 
                 str(rate.activation_energy / CALORIES_CONSTANT)
                 ]
    return '[' + ', '.join(arrhenius) + ']'


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
    if falloff_function == 'Troe':
        falloff_string = ('Troe(' +
                          f'A = {parameters[0]}' +
                          f', T3 = {parameters[1]}' +
                          f', T1 = {parameters[2]}' +
                          f', T2 = {parameters[3]})'
                          )
    elif falloff_function == 'SRI':
        falloff_string = ('SRI(' + 
                          f'A = {parameters[0]}' +
                          f', B = {parameters[1]}' +
                          f', C = {parameters[2]}' +
                          f', D = {parameters[3]}' +
                          f', E = {parameters[4]})'
                          )
    else:
        raise NotImplementedError(f'Falloff function not supported: {falloff_function}')

    return falloff_string


def build_efficiencies(efficiencies, species_names, default_efficiency=1.0):
    """Creates line with list of third-body species efficiencies.

    Parameters
    ----------
    efficiencies : dict
        Dictionary of species efficiencies
    species_names : dict of str
        List of all species names
    default_efficiency : float, optional
        Default efficiency for all species; will be 0.0 for reactions with explicit third body

    Returns
    -------
    str
        Line with list of efficiencies

    """
    # Reactions with a default_efficiency of 0 and a single entry in the efficiencies dict
    # have an explicit third body specified.
    if len(efficiencies) == 1 and not default_efficiency:
        return ''

    reduced_efficiencies = {s:efficiencies[s] for s in efficiencies if s in species_names}
    return '  '.join([f'{s}:{v}' for s, v in reduced_efficiencies.items()])


def write(solution, output_filename='', path=''):
    """Function to write cantera solution object to cti file.

    Parameters
    ----------
    solution : cantera.Solution
        Model to be written
    output_filename : str, optional
        Name of file to be written; if not provided, use ``solution.name``
    path : str, optional
        Path for writing file.

    Returns
    -------
    output_filename : str
        Name of output model file (.cti)

    Examples
    --------
    >>> gas = cantera.Solution('gri30.cti')
    >>> soln2cti.write(gas, 'copy_gri30.cti')
    copy_gri30.cti

    """
    if output_filename:
        output_filename = os.path.join(path, output_filename)
    else:
        output_filename = os.path.join(path, f'{solution.name}.cti')
    
    if os.path.isfile(output_filename):
        os.remove(output_filename)

    with open(output_filename, 'w') as the_file: 

        # Write title block to file
        the_file.write(section_break('CTI file converted from solution object'))
        the_file.write('units(length = "cm", time = "s",' +
                       ' quantity = "mol", act_energy = "cal/mol")' +
                       '\n\n'
                       )

        # Write Phase definition to file
        element_names = ' '.join(solution.element_names)
        species_names = fill(
            ' '.join(solution.species_names), 
            width=55, 
            subsequent_indent=19*' ',
            break_long_words=False,
            break_on_hyphens=False
            )

        the_file.write(
            f'ideal_gas(name = "{os.path.splitext(os.path.basename(output_filename))[0]}", \n' +
            indent[5] + f'elements = "{element_names}", \n' +
            indent[5] + f'species = """ {species_names} """, \n' +
            indent[5] + f'reactions = "all", \n' +
            indent[5] + f'initial_state = state(temperature = {solution.T}, ' +
            f'pressure = {solution.P})   )\n\n'
            )

        # Write species data to file
        the_file.write(section_break('Species data'))
        
        for species in solution.species():
            # build strings with low- and high-temperature 7 NASA coefficients
            nasa_range_low = f'[{species.thermo.min_temp}, {species.thermo.coeffs[0]}]'
            nasa_coeffs_low = [f'{c:.9e}' for c in species.thermo.coeffs[8:15]]
            nasa_coeffs_low = fill(
                '[' + ', '.join(nasa_coeffs_low) + ']', 
                width=50, 
                subsequent_indent=16*' ',
                break_long_words=False,
                break_on_hyphens=False
                )
            
            nasa_range_high = f'[{species.thermo.coeffs[0]}, {species.thermo.max_temp}]'
            nasa_coeffs_high = [f'{c:.9e}' for c in species.thermo.coeffs[1:8]]
            nasa_coeffs_high = fill(
                '[' + ', '.join(nasa_coeffs_high) + ']', 
                width=50, 
                subsequent_indent=16*' ',
                break_long_words=False,
                break_on_hyphens=False
                )

            composition = ', '.join([f'{s}:{int(v)}' for s, v in species.composition.items()])

            # start writing composition and thermo data
            species_string = (
                f'species(name = "{species.name}",\n' +
                f'{indent[4]}atoms = "{composition}", \n' +
                f'{indent[4]}thermo = (\n' +
                f'{indent[7]}NASA(   {nasa_range_low}, {nasa_coeffs_low}  ),\n' +
                f'{indent[7]}NASA(   {nasa_range_high}, {nasa_coeffs_high}  )\n' +
                f'{indent[15]}),\n'
                )

            #check if species has defined transport data, and write that if so
            if species.transport:
                species_string += (
                    f'    transport = gas_transport(\n' +
                    indent[15] + f'geom = "{species.transport.geometry}",\n' +
                    indent[15] + f'diam = {species.transport.diameter * 1e10}, \n' +
                    indent[15] + f'well_depth = {species.transport.well_depth / ct.boltzmann}, \n' +
                    indent[15] + f'polar = {species.transport.polarizability * 1e30}, \n' +
                    indent[15] + f'rot_relax = {species.transport.rotational_relaxation}'
                    )
                if species.transport.dipole != 0:
                    dipole = species.transport.dipole / DEBEYE_CONVERSION
                    species_string += ', \n' + indent[15] + f'dipole= {dipole}'
                species_string += ')\n'
            
            species_string += '       )\n\n'
            the_file.write(species_string)

        # Write reactions to file
        the_file.write(section_break('Reactions'))

        # write data for each reaction
        for idx, reaction in enumerate(solution.reactions()):

            reaction_string = f'#  Reaction {idx + 1}\n'

            if type(reaction) == ct.ElementaryReaction:
                arrhenius = build_arrhenius(reaction.rate, 
                                            sum(reaction.reactants.values()), 
                                            ct.ElementaryReaction
                                            )
                reaction_string += f'reaction( "{reaction.equation}",  [{arrhenius}]'

            elif type(reaction) == ct.ThreeBodyReaction:
                arrhenius = build_arrhenius(reaction.rate, 
                                            sum(reaction.reactants.values()), 
                                            ct.ThreeBodyReaction
                                            )
                reaction_string += f'three_body_reaction( "{reaction.equation}",  [{arrhenius}]'

                # potentially trimmed efficiencies list
                efficiencies_str = build_efficiencies(reaction.efficiencies, solution.species_names)
                if efficiencies_str:
                    reaction_string += f',\n{indent[9]}efficiencies = "{efficiencies_str}"'
                                    
            elif type(reaction) == ct.FalloffReaction:
                arrhenius_high = build_falloff_arrhenius(
                    reaction.high_rate, sum(reaction.reactants.values()), 
                    ct.FalloffReaction, 'high'
                    )
                arrhenius_low = build_falloff_arrhenius(
                    reaction.low_rate, sum(reaction.reactants.values()), 
                    ct.FalloffReaction, 'low'
                    )

                reaction_string += (f'falloff_reaction( "{reaction.equation}",\n' +
                                    f'{indent[9]}kf = {arrhenius_high},\n' +
                                    f'{indent[9]}kf0 = {arrhenius_low}'
                                    )
                
                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    reaction_string += ',\n' + indent[9] + 'falloff = ' + falloff_str
                
                # potentially trimmed efficiencies list
                efficiencies_str = build_efficiencies(
                    reaction.efficiencies, solution.species_names, reaction.default_efficiency
                    )
                if efficiencies_str:
                    reaction_string += f',\n{indent[9]}efficiencies = "{efficiencies_str}"'
            
            elif type(reaction) == ct.ChemicallyActivatedReaction:
                arrhenius_high = build_falloff_arrhenius(
                    reaction.high_rate, sum(reaction.reactants.values()), 
                    ct.ChemicallyActivatedReaction, 'high'
                    )
                arrhenius_low = build_falloff_arrhenius(
                    reaction.low_rate, sum(reaction.reactants.values()), 
                    ct.ChemicallyActivatedReaction, 'low'
                    )

                reaction_string += (f'chemically_activated_reaction( "{reaction.equation}",\n' +
                                    f'{indent[15]} kLow = {arrhenius_low},\n' +
                                    f'{indent[15]} kHigh = {arrhenius_high}'
                                    )
                
                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    reaction_string += ',\n' + indent[9] + 'falloff = ' + falloff_str
                
                # potentially trimmed efficiencies list
                efficiencies_str = build_efficiencies(
                    reaction.efficiencies, solution.species_names, reaction.default_efficiency
                    )
                if efficiencies_str:
                    reaction_string += f',\n{indent[8]} efficiencies = "{efficiencies_str}"'

            elif type(reaction) == ct.PlogReaction:
                reaction_string += f'pdep_arrhenius( "{reaction.equation}",\n'

                rates = []
                for rate in reaction.rates:
                    pressure = f'{rate[0] / ct.one_atm}'
                    arrhenius = build_arrhenius(rate[1], 
                                                sum(reaction.reactants.values()), 
                                                ct.PlogReaction
                                                )
                    rates.append(f'{indent[15]}[({pressure}, "atm"), {arrhenius}]')
                # want to get the list of rates with a comma and newline between each entry, 
                # but not at the end.
                reaction_string += ',\n'.join(rates)
            
            elif type(reaction) == ct.ChebyshevReaction:
                reaction_string += f'chebyshev_reaction( "{reaction.equation}",\n'
                
                # need to modify first coefficient
                coeffs = reaction.coeffs[:]
                coeffs[0][0] -= math.log10(1e-3 ** (sum(reaction.reactants.values()) - 1))

                coeffs_strings = []
                for coeff_row in coeffs:
                    coeffs_strings.append('[' + ', '.join([f'{c:.6e}' for c in coeff_row]) + ']')
                coeffs_string = f',\n{indent[15] + indent[9]}'.join(coeffs_strings)
                
                reaction_string += (
                    f'{indent[15]} Tmin={reaction.Tmin}, Tmax={reaction.Tmax},\n' +
                    f'{indent[15]} Pmin=({reaction.Pmin / ct.one_atm}, "atm"), Pmax=({reaction.Pmax / ct.one_atm}, "atm"),\n' +
                    f'{indent[15]} coeffs=[{coeffs_string}]'
                    )

            else:
                raise NotImplementedError(f'Unsupported reaction type: {type(reaction)}')
            
            if reaction.duplicate:
                reaction_string += ',\n' + indent[8] + 'options = "duplicate"'
                                
            reaction_string += ')\n\n'
            the_file.write(reaction_string)
    
    return output_filename
