"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff, Plog and ThreeBody Reactions
Cantera development version 2.3.0a2 required
"""

import os
from textwrap import fill

import cantera as ct

from .name_trim import name_trim

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

#1 debye = d coulomb-meters
DEBEYE_CONVERSION = 3.33564e-30


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
    return '[' + ', '.join(arrhenius) + ']'


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
    if falloff_function == ct.TroeFalloff:
        falloff_string = ('Troe(' +
                          f'A = {parameters[0]}' +
                          f', T3 = {parameters[1]}' +
                          f', T1 = {parameters[2]}' +
                          f', T2 = {parameters[3]})'
                          )
    elif falloff_function == ct.SriFalloff:
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
    >>> gas = cantera.Solution('gri30.cti')
    >>> soln2cti.write(gas)
    reduced_gri30.cti

    """
    # Remove extension from filename
    input_file_name = os.path.splitext(os.path.basename(solution.name))[0]
    output_file_name = f'reduced_{input_file_name}.cti'

    with open(output_file_name, 'w') as the_file: 

        # Write title block to file
        the_file.write(section_break('CTI file converted from solution object'))
        the_file.write('units(length = "cm", time = "s",' +
                       ' quantity = "mol", act_energy = "cal/mol")' +
                       '\n\n'
                       )

        # Write Phase definition to file
        element_names = ' '.join(solution.element_names)
        species_names = fill(' '.join(solution.species_names), width=55, subsequent_indent=19*' ')

        the_file.write(
            f'ideal_gas(name = "reduced_{input_file_name}", \n' +
            f'     elements = "{element_names}", \n' +
            f'     species = """ {species_names} """, \n' +
            f'     reactions = "all", \n' +
            f'     initial_state = state(temperature = {solution.T}, ' +
            f'pressure = {solution.P})   )\n\n'
            )

        # Write species data to file
        the_file.write(section_break('Species data'))
        
        for species in solution.species():
            # build strings with low- and high-temperature 7 NASA coefficients
            nasa_range_low = f'[{species.thermo.min_temp}, {species.thermo.coeffs[0]}]'
            nasa_coeffs_low = ["{:.9e}".format(c) for c in species.thermo.coeffs[8:15]]
            nasa_coeffs_low = fill('[' + ', '.join(nasa_coeffs_low) + ']', width=50, subsequent_indent=16*' ')
            
            nasa_range_high = f'[{species.thermo.coeffs[0]}, {species.thermo.max_temp}]'
            nasa_coeffs_high = ["{:.9e}".format(c) for c in species.thermo.coeffs[1:8]]
            nasa_coeffs_high = fill('[' + ', '.join(nasa_coeffs_high) + ']', width=50, subsequent_indent=16*' ')

            composition = ', '.join([f'{s}:{int(v)}' for s, v in species.composition.items()])

            # start writing composition and thermo data
            species_string = (
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
                species_string += (
                    f'    transport = gas_transport(\n' +
                    indent + f'geom = "{species.transport.geometry}",\n' +
                    indent + f'diam = {species.transport.diameter * 1e10}, \n' +
                    indent + f'well_depth = {species.transport.well_depth / ct.boltzmann}, \n' +
                    indent + f'polar = {species.transport.polarizability * 1e30}, \n' +
                    indent + f'rot_relax = {species.transport.rotational_relaxation}'
                    )
                if species.transport.dipole != 0:
                    dipole = species.transport.dipole / DEBEYE_CONVERSION
                    species_string += ', \n' + indent + f'dipole= {dipole}'
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
                reaction_string += f'reaction( "{reaction.equation}",  {arrhenius}'

            elif type(reaction) == ct.ThreeBodyReaction:
                arrhenius = build_arrhenius(reaction.rate, 
                                            sum(reaction.reactants.values()), 
                                            ct.ThreeBodyReaction
                                            )
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
                arrhenius_high = build_falloff_arrhenius(
                    reaction.rate, sum(reaction.reactants.values()), 
                    ct.FalloffReaction, 'high'
                    )
                arrhenius_low = build_falloff_arrhenius(
                    reaction.rate, sum(reaction.reactants.values()), 
                    ct.FalloffReaction, 'low'
                    )

                reaction_string += (f'falloff_reaction( "{reaction.equation}",\n' +
                                    f'         kf = {arrhenius_high},\n' +
                                    f'         kf0 = {arrhenius_low}'
                                    )
                
                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    reaction_string += ',\n' + '         falloff = ' + falloff_str
                
                # trims efficiencies list
                reduced_efficiencies = {s:reaction.efficiencies[s] 
                                        for s in reaction.efficiencies 
                                        if s in solution.species_names
                                        }
                efficiencies_str = '  '.join([f'{s}:{v}' for s, v in reduced_efficiencies.items()])
                if efficiencies_str:
                    reaction_string += f',\n         efficiencies = " {efficiencies_str} "'
            
            elif type(reaction) == ct.ChemicallyActivatedReaction:
                arrhenius_high = build_falloff_arrhenius(
                    reaction.rate, sum(reaction.reactants.values()), 
                    ct.ChemicallyActivatedReaction, 'high'
                    )
                arrhenius_low = build_falloff_arrhenius(
                    reaction.rate, sum(reaction.reactants.values()), 
                    ct.ChemicallyActivatedReaction, 'low'
                    )

                reaction_string += (f'chemically_activated_reaction( "{reaction.equation}",\n' +
                                    f'                              kLow = {arrhenius_low},\n' +
                                    f'                              kHigh = {arrhenius_high}'
                                    )
                
                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    reaction_string += ',\n' + '         falloff = ' + falloff_str
                
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

                rates = []
                for rate in reaction.rates:
                    pressure = f'{rate[0] / ct.one_atm}'
                    arrhenius = build_arrhenius(rate[1], 
                                                sum(reaction.reactants.values()), 
                                                ct.PlogReaction
                                                )
                    rates.append(f'               [({pressure}, "atm"), {arrhenius}')
                # want to get the list of rates with a comma and newline between each entry, 
                # but not at the end.
                reaction_string += ',\n'.join(rates)
            
            elif type(reaction) == ct.ChebyshevReaction:
                reaction_string += f'chebyshev_reaction( "{reaction.equation}",\n'

                coeffs_strings = []
                for coeff_row in reaction.coeffs:
                    coeffs_strings.append('[' + ', '.join([f'{c:.6e}' for c in coeff_row]) + ']')
                coeffs_string = ',\n                           '.join(coeffs_strings)
                
                reaction_string += (
                    f'Tmin={reaction.Tmin}, Tmax={reaction.Tmax},\n' +
                    f'Pmin=({reaction.Pmin / ct.one_atm}, "atm"), Pmax=({reaction.Pmax / ct.one_atm}, "atm"),\n' +
                    f'coeffs=[{coeffs_string}]'
                    )

            else:
                raise NotImplementedError(f'Unsupported reaction type: {type(reaction)}')
            
            if reaction.duplicate:
                reaction_string += ',\n        options = "duplicate"'
                                
            reaction_string += ')\n\n'
            the_file.write(reaction_string)
    
    return output_file_name
