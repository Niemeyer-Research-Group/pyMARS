# This file is part of Frhodo. Copyright © 2020, UChicago Argonne, LLC
# and licensed under BSD-3-Clause. See License.txt in the top-level 
# directory for license and copyright information.

''' 
Adapted from Kyle Niemeyer's pyMARS Jul 24, 2019

Writes a solution object to a chemkin inp file
currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required

KE Niemeyer, CJ Sung, and MP Raju. Skeletal mechanism generation for surrogate fuels using directed relation graph with error propagation and sensitivity analysis. Combust. Flame, 157(9):1760--1770, 2010. doi:10.1016/j.combustflflame.2009.12.022
KE Niemeyer and CJ Sung. On the importance of graph search algorithms for DRGEP-based mechanism reduction methods. Combust. Flame, 158(8):1439--1443, 2011. doi:10.1016/j.combustflflame.2010.12.010.
KE Niemeyer and CJ Sung. Mechanism reduction for multicomponent surrogates: A case study using toluene reference fuels. Combust. Flame, in press, 2014. doi:10.1016/j.combustflame.2014.05.001
TF Lu and CK Law. Combustion and Flame, 154:153--163, 2008. doi:10.1016/j.combustflame.2007.11.013

'''

import os, pathlib
from textwrap import fill
from collections import Counter

import cantera as ct

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

# number of calories in 1000 Joules
CALORIES_CONSTANT = 4184.0

# Conversion from 1 debye to coulomb-meters
DEBEYE_CONVERSION = 3.33564e-30

def reorder_reaction_equation(solution, reaction):
    # Split Reaction Equation
    rxn_eqn = reaction.equation
    for reaction_direction in [' <=> ', ' <= ', ' => ']:
        if reaction_direction in rxn_eqn:
            break
    for third_body in [' (+M)', ' + M', '']: # search eqn for third body
        if third_body in rxn_eqn:            # if reaches '', doesn't exist
            break

    # Sort and apply to reaction equation
    reaction_txt = []
    reaction_split = {'reactants': reaction.reactants, 
                      'products': reaction.products}
    for n, (reaction_side, species) in enumerate(reaction_split.items()):
        species_weights = []
        for key in species.keys():
            index = solution.species_index(key)
            species_weights.append(solution.molecular_weights[index])
        
        # Append coefficient to species
        species_list = []
        for species_text, coef in species.items():
            if coef == 1.0:
                species_list.append(species_text)
            else:
                species_list.append(f'{coef:.0f} {species_text}')
                
        species = species_list
        
        # Reorder species based on molecular weights
        species = [x for y, x in sorted(zip(species_weights, species))][::-1]
        reaction_txt.append(' + '.join(species) + third_body)
    
    reaction_txt = reaction_direction.join(reaction_txt)
    
    return reaction_txt


def match_reaction(solution, yaml_rxn):
    yaml_rxn = {'eqn': yaml_rxn}
    
    for reaction_direction in [' <=> ', ' <= ', ' => ']:
        if reaction_direction in yaml_rxn['eqn']:
            break
    for third_body in [' (+M)', ' + M', '']: # search eqn for third body
        if third_body in yaml_rxn['eqn']:    # if reaches '', doesn't exist
            break
    
    yaml_rxn_split = yaml_rxn['eqn'].split(reaction_direction)   
    for i, side in zip([0, 1], ['reac', 'prod']):
        yaml_rxn[side] = {}
        species = yaml_rxn_split[i].replace(third_body, '').split(' + ')
        yaml_rxn[side].update(Counter(species))
    
    for rxn in solution.reactions():
        if (rxn.reactants == yaml_rxn['reac'] and    
            rxn.products == yaml_rxn['prod'] and 
            third_body in str(rxn)):
                
            return str(rxn) # return rxn if match
    
    return yaml_rxn['eqn']  # returns yaml_str if no match


def get_notes(path=None, solution=None):
    """Get notes by parsing input mechanism in yaml format
    Parameters
    ----------
    path : path or str, optional
        Path of yaml file used as input in order to parse for notes
    solution : 
    """
    
    note = {'header': [], 'species_thermo': {}, 'species': {}, 'reaction': {}}
    
    if path is None: return note
    
    with open(path, 'r') as yaml_file:
        data = yaml.load(yaml_file, yaml.RoundTripLoader)
    
    # Header note
    if 'description' in data:
        note['header'] = data['description']
    else:
        note['header'] = ''
    
    # Species and thermo_species notes
    for species in data['species']:
        if 'note' in species:
            note['species'][species['name']] = species['note']
        else:
            note['species'][species['name']] = ''

        if 'note' in species['thermo']:
            note['species_thermo'][species['name']] = species['thermo']['note']
        else:
            note['species_thermo'][species['name']] = ''
    
    if 'reactions' in data:
        for rxn in data['reactions']:
            ct_rxn_eqn = match_reaction(solution, rxn['equation'])
            if 'note' in rxn:
                note['reaction'][ct_rxn_eqn] = '! ' + rxn['note'].replace('\n', '\n! ')
            else:
                note['reaction'][ct_rxn_eqn] = ''
    
    return note
    

def eformat(f, precision=7, exp_digits=3):
    s = f"{f: .{precision}e}"
    if s == ' inf' or s == '-inf':
        return s
    else:
        mantissa, exp = s.split('e')
        exp_digits += 1 # +1 due to sign
        return f"{mantissa}E{int(exp):+0{exp_digits}}" 
  
 
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
    
    activation_energy = rate.activation_energy / CALORIES_CONSTANT
    arrhenius = [f'{eformat(pre_exponential_factor)}',
                 f'{eformat(rate.temperature_exponent)}', 
                 f'{eformat(activation_energy)}']
    return '   '.join(arrhenius)


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

    activation_energy = rate.activation_energy / CALORIES_CONSTANT
    arrhenius = [f'{eformat(pre_exponential_factor)}', 
                 f'{eformat(rate.temperature_exponent)}', 
                 f'{eformat(activation_energy)}'
                 ]
    return '   '.join(arrhenius)


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
        falloff = [f'{eformat(f)}'for f in parameters]
        falloff_string = f"TROE / {'   '.join(falloff)} /\n"
    elif falloff_function == 'SRI':
        falloff = [f'{eformat(f)}'for f in parameters]
        falloff_string = f"SRI / {'   '.join(falloff)} /\n"
    else:
        raise NotImplementedError(f'Falloff function not supported: {falloff_function}')

    return falloff_string


def species_data_text(species_list, note):
    max_species_len = max([len(s) for s in species_list])
    if note:
        max_species_len = max([16, max_species_len])
        species_txt = []
        for species in species_list:
            text = f'{species:<{max_species_len}}   ! {note[species]}\n'
            species_txt.append(text)
        
        species_txt = ''.join(species_txt)
        
    else:
        species_names = [f"{s:<{max_species_len}}" for s in species_list]
        species_names = fill(
            '  '.join(species_names), 
            width=72,   # max length is 16, this gives 4 species per line
            break_long_words=False,
            break_on_hyphens=False
            )
        
        species_txt = f'{species_names}\n'
        
    text = ('SPECIES\n' + 
            species_txt + 
            'END\n\n\n')
    
    return text


def thermo_data_text(species_list, note, input_type='included'):
    """Returns thermodynamic data in Chemkin-format file.
    Parameters
    ----------
    species_list : list of cantera.Species
        List of species objects
    input_type : str, optional
        'included' if thermo will be printed in mech file, 'file' otherwise
    """
    
    if input_type == 'included':
        thermo_text = ['THERMO ALL\n' +  
                       '   300.000  1000.000  6000.000\n']
    else:
        thermo_text = ['THERMO\n' +  
                       '   300.000  1000.000  6000.000\n']

    # write data for each species in the Solution object
    for species in species_list:
        composition_string = ''.join([f'{s:2}{int(v):>3}' 
                                      for s, v in species.composition.items()
                                      ])

        # first line has species name, space for notes/date, elemental composition,
        # phase, thermodynamic range temperatures (low, high, middle), and a "1"
        # total length should be 80
        
        # attempt to split note and comment
        if len(note[species.name].split('\n', 1)) == 1:
            comment = ''
            comment_str = ''
            note_str = note[species.name]
        else:
            comment = '!\n'
            note_str, comment_str = note[species.name].split('\n', 1)
        
        if len(f'{species.name} {note_str}') > 24:
            comment_str += '\n' + note_str
            note_str = ''
            
        comment_str = comment_str.replace('\n', '\n! ')
        comment = f'{comment}! {comment_str}'
        
        name_and_note = f'{species.name} {note_str}'
        species_string = (comment + '\n' +
            f'{name_and_note:<24}' + # name and date/note field
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
        
        thermo_text.append(species_string)
    
    if input_type == 'included':
        thermo_text.append('END\n\n\n')
    else:
        thermo_text.append('END\n')
    
    return ''.join(thermo_text)
    

def write_transport_data(species_list, filename='generated_transport.dat'):
    """Writes transport data to Chemkin-format file.
    Parameters
    ----------
    species_list : list of cantera.Species
        List of species objects
    filename : path or str, optional
        Filename for new Chemkin transport database file
    """
    geometry = {'atom': '0', 'linear': '1', 'nonlinear': '2'}
    with open(filename, 'w') as trans_file:

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
            
            trans_file.write(species_string)


def write(solution, output_path='', input_yaml='',
          skip_thermo=False, same_file_thermo=True, 
          skip_transport=False):
    """Writes Cantera solution object to Chemkin-format file.
    Parameters
    ----------
    solution : cantera.Solution
        Model to be written
    output_path : path or str, optional
        Path of file to be written; if not provided, use cd / 'solution.name'
    input_yaml : path or str, optional
        Path of yaml file used as input in order to parse for notes
    skip_thermo : bool, optional
        Flag to skip writing thermo data
    same_file_thermo : bool, optional
        Flag to write thermo data in the mechanism file
    skip_transport : bool, optional
        Flag to skip writing transport data in separate file
    Returns
    -------
    output_file_name : str
        Name of output model file (.ck)
    Examples
    --------
    >>> gas = cantera.Solution('gri30.cti')
    >>> soln2ck.write(gas)
    reduced_gri30.ck
    """
    if output_path:
        if not isinstance(output_path, pathlib.PurePath):
            output_path = pathlib.Path(output_path)
    else:
        main_path = pathlib.Path.cwd()
        output_path = main_path / f'{solution.name}.ck'
    
    if output_path.is_file():
        output_path.unlink()
       
    main_path = output_path.parents[0]
    basename = output_path.stem
    output_files = [output_path]
    
    if input_yaml:
        if not isinstance(input_yaml, pathlib.PurePath):
            input_yaml = pathlib.Path(input_yaml)
            
        note = get_notes(input_yaml, solution)
    else:    
        note = get_notes()
    
    with open(output_path, 'w') as mech_file:
        # Write title block to file
        if note['header']:
            note["header"] = note['header'].replace('\n', '\n! ')
            mech_file.write(f'! {note["header"]}\n! \n')
        mech_file.write('! Chemkin file converted from Cantera solution object\n! \n\n')

        # write species and element lists to file
        element_names = '  '.join(solution.element_names)
        mech_file.write(
            'ELEMENTS\n' + 
            f'{element_names}\n' +
            'END\n\n\n'
            )
        
        mech_file.write(species_data_text(solution.species_names, note['species']))

        # Write thermo to file
        if not skip_thermo and same_file_thermo:
            mech_file.write(thermo_data_text(solution.species(), note['species_thermo'], 
                            input_type='included'))
            
        # Write reactions to file
        max_rxn_width = 3 + max([len(rxn.equation) for rxn in solution.reactions()] + [48])
        
        mech_file.write('REACTIONS  CAL/MOLE  MOLES\n')
        # Write data for each reaction in the Solution Object
        for n, reaction in enumerate(solution.reactions()):
            reaction_equation = str(reaction)
            
            reaction_string = ''
            if reaction_equation in note['reaction']:
                rxn_note = note['reaction'][reaction_equation]
                rxn_note = rxn_note.rsplit('\n! ', 1)
                if len(rxn_note) > 1:
                    reaction_string = f'{rxn_note[0]}\n'
                    after_eqn_text = rxn_note[-1].strip()
                    rxn_note[-1] = f'! {after_eqn_text}'
            else:
                rxn_note = ['']
            
            reaction_equation = reorder_reaction_equation(solution, reaction)
            reaction_string += f'{reaction_equation:<{max_rxn_width}}'

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

            reaction_string += f'{arrhenius}    {rxn_note[-1]}\n'
            
            # now write any auxiliary information for the reaction
            if type(reaction) == ct.FalloffReaction:
                # for falloff reaction, need to write low-pressure limit Arrhenius expression
                arrhenius = build_falloff_arrhenius(
                    reaction.low_rate, 
                    sum(reaction.reactants.values()), 
                    ct.FalloffReaction,
                    'low'
                    )
                reaction_string += f'{"LOW /   ".rjust(max_rxn_width)}{arrhenius} /\n'

                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    width = max_rxn_width - 10 - 15*(reaction.falloff.parameters.size - 3)
                    reaction_string += f'{"".ljust(width)}{falloff_str}'

            elif type(reaction) == ct.ChemicallyActivatedReaction:
                # for chemically activated reaction, need to write high-pressure expression
                arrhenius = build_falloff_arrhenius(
                    reaction.low_rate, 
                    sum(reaction.reactants.values()), 
                    ct.ChemicallyActivatedReaction,
                    'high'
                    )
                reaction_string += f'HIGH'
                reaction_string += f'{"HIGH /   ".rjust(max_rxn_width)}{arrhenius} /\n'

                # need to print additional falloff parameters if present
                if reaction.falloff.parameters.size > 0:
                    falloff_str = build_falloff(reaction.falloff.parameters, reaction.falloff.type)
                    width = max_rxn_width - 10 - 15*(reaction.falloff.parameters.size - 3)
                    reaction_string += f'{"".ljust(width)}{falloff_str}'

            elif type(reaction) == ct.PlogReaction:
                # just need one rate per line
                for rate in reaction.rates:
                    pressure = f'{eformat(rate[0] / ct.one_atm)}'
                    arrhenius = build_arrhenius(rate[1], 
                                                sum(reaction.reactants.values()), 
                                                ct.PlogReaction
                                                )
                    reaction_string += (f'{"PLOG / ".rjust(max_rxn_width-18)}'
                                        f'{pressure}   {arrhenius} /\n')

            elif type(reaction) == ct.ChebyshevReaction:
                reaction_string += (
                    f'TCHEB / {reaction.Tmin}  {reaction.Tmax} /\n' +
                    f'PCHEB / {reaction.Pmin / ct.one_atm}  {reaction.Pmax / ct.one_atm} /\n' +
                    f'CHEB / {reaction.nTemperature}  {reaction.nPressure} /\n'
                    )
                for coeffs in reaction.coeffs:
                    coeffs_row = ' '.join([f'{c:.6e}' for c in coeffs])
                    reaction_string += f'CHEB / {coeffs_row} /\n'
            
            # need to trim and print third-body efficiencies, if present
            if type(reaction) in [ct.ThreeBodyReaction, ct.FalloffReaction, 
                                  ct.ChemicallyActivatedReaction
                                  ]:
                # trims efficiencies list
                reduced_efficiencies = {s:reaction.efficiencies[s] 
                                        for s in reaction.efficiencies
                                        if s in solution.species_names
                                        }
                efficiencies_str = '  '.join([f'{s}/ {v:.3f}/' for s, v in reduced_efficiencies.items()])
                if efficiencies_str:
                    reaction_string += '   ' + efficiencies_str + '\n'
            
            if reaction.duplicate:
                reaction_string += '   DUPLICATE\n'
                                
            mech_file.write(reaction_string)

        mech_file.write('END')

    # write thermo data
    if not skip_thermo and not same_file_thermo:
        therm_path = main_path / f'{basename}.therm'
        with open(therm_path, 'w') as thermo_file:
            thermo_file.write(thermo_data_text(solution.species(), input_type='file'))
        output_files.append(therm_path)

    # TODO: more careful check for presence of transport data?
    if not skip_transport and all(sp.transport for sp in solution.species()):
        trans_path = main_path / f'{basename}_tranport.dat'
        write_transport_data(solution.species(), trans_path)
        output_files.append(trans_path)

    return output_files
