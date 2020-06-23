"""writes a solution object to a chemkin inp file

currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required
"""

import os
from textwrap import fill

import cantera as ct

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
DEBEYE_CONVERSION = 3.33564e-30


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
                 f'{(rate.activation_energy / CALORIES_CONSTANT):.3e}'
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
    if falloff_function == 'Troe':
        falloff_string = ('TROE / ' +
                          f'{parameters[0]}  {parameters[1]}  '
                          f'{parameters[2]}  {parameters[3]} /\n'
                          )
    elif falloff_function == 'SRI':
        falloff_string = ('SRI / ' + 
                          f'{parameters[0]}  {parameters[1]}  ' +
                          f'{parameters[2]}  {parameters[3]}  {parameters[4]} /\n'
                          )
    else:
        raise NotImplementedError(f'Falloff function not supported: {falloff_function}')

    return falloff_string


def write_thermo_data(species_list, filename='generated_thermo.dat'):
    """Writes thermodynamic data to Chemkin-format file.

    Parameters
    ----------
    species_list : list of cantera.Species
        List of species objects
    filename : str, optional
        Filename for new Chemkin thermodynamic database file

    """
    with open(filename, 'w') as the_file:

        the_file.write('THERMO\n' + 
                       '   300.000  1000.000  5000.000\n'
                       )

        # write data for each species in the Solution object
        for species in species_list:
            
            composition_string = ''.join([f'{s:2}{int(v):>3}' 
                                          for s, v in species.composition.items()
                                          ])

            # first line has species name, space for notes/date, elemental composition,
            # phase, thermodynamic range temperatures (low, high, middle), and a "1"
            # total length should be 80
            species_string = (
                f'{species.name:<18}' + 
                6*' ' + # date/note field
                f'{composition_string:<20}' +
                'G' + # only supports gas phase
                f'{species.thermo.min_temp:10.3f}' +
                f'{species.thermo.max_temp:10.3f}' +
                f'{species.thermo.coeffs[0]:8.2f}' +
                6*' ' + # unused atomic symbols/formula, and blank space
                '1\n'
                )
            
            # second line has first five coefficients of high-temperature range,
            # ending with a "2" in column 79
            species_string += (
                ''.join([f'{c:15.8e}' for c in species.thermo.coeffs[1:6]]) +
                '    ' +
                '2\n'
            )
            
            # third line has the last two coefficients of the high-temperature range,
            # first three coefficients of low-temperature range, and "3"
            species_string += (
                ''.join([f'{c:15.8e}' for c in species.thermo.coeffs[6:8]]) +
                ''.join([f'{c:15.8e}' for c in species.thermo.coeffs[8:11]]) +
                '    ' +
                '3\n'
            )

            # fourth and last line has the last four coefficients of the
            # low-temperature range, and "4"
            
            species_string += (
                ''.join([f'{c:15.8e}' for c in species.thermo.coeffs[11:15]]) +
                19*' ' +
                '4\n'
            )
            
            the_file.write(species_string)

        the_file.write('END\n')


def write_transport_data(species_list, filename='generated_transport.dat'):
    """Writes transport data to Chemkin-format file.

    Parameters
    ----------
    species_list : list of cantera.Species
        List of species objects
    filename : str, optional
        Filename for new Chemkin transport database file

    """
    geometry = {'atom': '0', 'linear': '1', 'nonlinear': '2'}

    with open(filename, 'w') as the_file:

        # write data for each species in the Solution object
        for species in species_list:
            
            # each line contains the species name, integer representing
            # geometry, Lennard-Jones potential well depth in K,
            # Lennard-Jones collision diameter in angstroms,
            # dipole moment in Debye,
            # polarizability in cubic angstroms, and
            # rotational relaxation collision number at 298 K.
            species_string = (
                f'{species.name:<16}' +
                f'{geometry[species.transport.geometry]:>4}' +
                f'{(species.transport.well_depth / ct.boltzmann):>10.3f}' + 
                f'{(species.transport.diameter * 1e10):>10.3f}' + 
                f'{(species.transport.dipole / DEBEYE_CONVERSION):>10.3f}' + 
                f'{(species.transport.polarizability * 1e30):>10.3f}' + 
                f'{species.transport.rotational_relaxation:>10.3f}' + 
                '\n'
            )
            
            the_file.write(species_string)


def write(solution, output_filename='', path='', 
          skip_thermo=False, skip_transport=False
          ):
    """Writes Cantera solution object to Chemkin-format file.

    Parameters
    ----------
    solution : cantera.Solution
        Model to be written
    output_filename : str, optional
        Name of file to be written; if not provided, use ``solution.name``
    path : str, optional
        Path for writing file.
    skip_thermo : bool, optional
        Flag to skip writing thermo data in separate file
    skip_transport : bool, optional
        Flag to skip writing transport data in separate file

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
    if output_filename:
        output_filename = os.path.join(path, output_filename)
    else:
        output_filename = os.path.join(path, f'{solution.name}.inp')
    
    if os.path.isfile(output_filename):
        os.remove(output_filename)
    
    with open(output_filename, 'w') as the_file:

        # Write title block to file
        the_file.write(f'!Chemkin file converted from solution object: {solution.name}\n\n')

        # write species and element lists to file
        element_names = '  '.join(solution.element_names)
        the_file.write(
            'ELEMENTS\n' + 
            f'{element_names}\n' +
            'END\n\n'
            )
        species_names = fill(
            '  '.join(solution.species_names), 
            width=60,
            break_long_words=False,
            break_on_hyphens=False
            )
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
            else:
                raise NotImplementedError(f'Unsupported reaction type: {type(reaction)}')

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
            
            # now write any auxiliary information for the reaction
            if type(reaction) == ct.FalloffReaction:
                # for falloff reaction, need to write low-pressure limit Arrhenius expression
                arrhenius = build_falloff_arrhenius(
                    reaction.low_rate, 
                    sum(reaction.reactants.values()), 
                    ct.FalloffReaction,
                    'low'
                    )
                reaction_string += f'LOW / {arrhenius} /\n'

                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    reaction_string += falloff_str

            elif type(reaction) == ct.ChemicallyActivatedReaction:
                # for chemically activated reaction, need to write high-pressure expression
                arrhenius = build_falloff_arrhenius(
                    reaction.low_rate, 
                    sum(reaction.reactants.values()), 
                    ct.ChemicallyActivatedReaction,
                    'high'
                    )
                reaction_string += f'HIGH / {arrhenius} /\n'

                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    reaction_string += falloff_str

            elif type(reaction) == ct.PlogReaction:
                # just need one rate per line
                for rate in reaction.rates:
                    pressure = f'{rate[0] / ct.one_atm}'
                    arrhenius = build_arrhenius(rate[1], 
                                                sum(reaction.reactants.values()), 
                                                ct.PlogReaction
                                                )
                    reaction_string += f'PLOG / {pressure} {arrhenius} /\n'

            elif type(reaction) == ct.ChebyshevReaction:
                reaction_string += (
                    f'TCHEB / {reaction.Tmin}  {reaction.Tmax} /\n' +
                    f'PCHEB / {reaction.Pmin / ct.one_atm}  {reaction.Pmax / ct.one_atm} /\n' +
                    f'CHEB / {reaction.nTemperature}  {reaction.nPressure} /\n'
                    )
                for coeffs in reaction.coeffs:
                    coeffs_row = ' '.join([f'{c:.6e}' for c in coeffs])
                    reaction_string += f'CHEB / {coeffs_row} /\n'

            if reaction.duplicate:
                reaction_string += 'DUPLICATE\n'
                                
            the_file.write(reaction_string)

        the_file.write('END')

    basename = os.path.splitext(output_filename)[0]
    outputs = [output_filename]

    # write thermo data
    if not skip_thermo:
        write_thermo_data(solution.species(), basename + '_thermo.dat')
        outputs.append(basename + '_thermo.dat')

    # TODO: more careful check for presence of transport data?
    if not skip_transport and all(sp.transport for sp in solution.species()):
        write_transport_data(solution.species(), basename + '_transport.dat')
        outputs.append(basename + '_transport.dat')

    return outputs
