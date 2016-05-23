"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required"""

#from identify_file_extension import readin
import os
import textwrap
from string import Template
import cantera as ct


def write(data_file, solution_objects):
    """Function to write cantera solution object to cti file.

    Parameters
    ----------
    data_file:
        Original cantera mechanism file.
    solution_objects:
        Trimmed cantera solution object.

    Returns
    -------
        Trimmed Mechanism file
    """


    input_file_name_stripped=os.path.splitext(data_file)[0]
    output_file_name=os.path.abspath('Output_Data_Files/'+ 'trimmed_' + input_file_name_stripped + '.cti')
    try:
        os.remove(output_file_name)
    except OSError:
        pass

    f=open(output_file_name, 'w+')


    #Read In trimmed Solution to Cantera Object
    trimmed_solution=solution_objects[1]
    solution_T=trimmed_solution.T
    solution_P=trimmed_solution.P

    """-------------------------------------------------------------------------
    Work Functions
    -------------------------------------------------------------------------"""

    def eliminate(input_string, char_to_replace, spaces='single'):
        for char in char_to_replace:
                    input_string= input_string.replace(char, "")
        if spaces == 'double':
                    input_string=input_string.replace(" ", "  ")
        return input_string

    def wrap(input_string):
        output_string= textwrap.fill(input_string, width=60, \
                                subsequent_indent= '                        ')
        return output_string
    def wrap_nasa(input_string):
        output_string= textwrap.fill(input_string, width=50, \
                                    subsequent_indent= '                ')
        return output_string

    def section_break(title):
        f.write('#'+ "-"*75 + '\n')
        f.write('#  ' + title +'\n')
        f.write('#'+ "-"*75 + '\n\n')

    def replace_multiple(input_string, replace_list):
        for a, b in replace_list.items():
            input_string= input_string.replace(a, b)
        return input_string

    def build_Arr(equation_object):
            A=str("{:.5E}".format(equation_object.rate.pre_exponential_factor))
            b=equation_object.rate.temperature_exponent
            E=equation_object.rate.activation_energy
            Arr=[ A, b, E]
            return str(Arr).replace("\'", "")

    def build_mod_Arr(equation_object, t_range):
        if t_range =='high':
            A=str("{:.5E}".format(equation_object.high_rate.pre_exponential_factor))
            b=equation_object.high_rate.temperature_exponent
            E=equation_object.high_rate.activation_energy
            Arr_high=[ A, b, E]
            return str(Arr_high).replace("\'", "")
        if t_range == 'low':
            A=str("{:.5E}".format(equation_object.low_rate.pre_exponential_factor))
            b=equation_object.low_rate.temperature_exponent
            E=equation_object.low_rate.activation_energy
            Arr_low=[ A, b, E]
            return str(Arr_low).replace("\'", "")

    def build_falloff(j):
        falloff_str=str(',\n        falloff = Troe(' +
                        'A = ' + str(j[0]) +
                        ', T3 = ' + str(j[1]) +
                        ', T1 = ' + str(j[2]) +
                        ', T2 = ' + str(j[3]) +')       )\n\n')
        return falloff_str


    """-------------------------------------------------------------------------
    Write Title Block to file
    -------------------------------------------------------------------------"""
    section_break('CTI File converted from Solution Object')

    unit_string="units(length = \"cm\", time = \"s\"," +\
                        " quantity = \"mol\", act_energy = \"cal/mol\")"
    f.write(unit_string + '\n\n')


    """-------------------------------------------------------------------------
    Write Phase definition to file
    -------------------------------------------------------------------------"""

    element_names=eliminate( str(trimmed_solution.element_names), \
                                        ['[', ']', '\'', ','])

    species_names=wrap(
                        eliminate(  str(trimmed_solution.species_names), \
                                    ['[', ']', '\'', ','], \
                                    spaces='double')
                        )
    phase_string= Template('ideal_gas(name = \"$input_file_name_stripped\", \n' +
                    '     elements = \"$elements\", \n' +
                    '     species =""" $species""", \n' +
                    '     reactions = \"all\", \n' +
                    '     initial_state = state(temperature = $solution_T, \n \
                            pressure= $solution_P)   )       \n\n')


    f.write(phase_string.substitute(elements=element_names, \
                            species=species_names,\
                            input_file_name_stripped=input_file_name_stripped,\
                            solution_T=solution_T, solution_P=solution_P))

    """-------------------------------------------------------------------------
    Write Species to file
    -------------------------------------------------------------------------"""

    section_break('Species data')

    #write data for each species in the Solution object
    for i, name in enumerate(trimmed_solution.species_names):

        #physical Constant
        k=ct.boltzmann #joules/kelvin, boltzmann constant
        d=3.33564e-30 #1 debye = d coulomb-meters

        species=trimmed_solution.species(i)
        name=trimmed_solution.species(i).name
        nasa_coeffs=trimmed_solution.species(i).thermo.coeffs
        replace_list_1= {'{':'\"',       '}':'\"',       '\'':'',
                    ':  ':':',      '.0':"",         ',':'',       ' ': '  '}

        #build 7-coeff NASA polynomial array
        nasa_coeffs_1=[]
        for j, k in enumerate(nasa_coeffs):

                coeff="{:.9e}".format(nasa_coeffs[j+8])
                nasa_coeffs_1.append(coeff)
                if j == 6:
                    nasa_coeffs_1=wrap_nasa(eliminate(str(  nasa_coeffs_1), \
                                                    {'\'':""}))
                    break
        nasa_coeffs_2=[]
        for j, k in enumerate(nasa_coeffs):

                coeff="{:.9e}".format(nasa_coeffs[j+1])
                nasa_coeffs_2.append(coeff)
                if j == 6:
                    nasa_coeffs_2=wrap_nasa(eliminate(str(  nasa_coeffs_2), \
                                                    {'\'':""}))
                    break

        #Species attributes from trimmed solution object
        composition = replace_multiple(str(species.composition), replace_list_1)
        nasa_range_1 = str([ species.thermo.min_temp, nasa_coeffs[0] ])
        nasa_range_2 = str([ nasa_coeffs[0], species.thermo.max_temp ])
        transport_geometry = species.transport.geometry
        diameter = str(species.transport.diameter*(10**10))
        well_depth = str(species.transport.well_depth/k)
        polar = str(species.transport.polarizability*10**30)
        rot_relax = str(species.transport.rotational_relaxation)
        dipole=str(species.transport.dipole/d)
        #string template for each species
        species_string=Template('species(name = "$name",\n' +
                        '    atoms = $composition, \n' +
                        '    thermo = (\n' +
                        '       NASA(   $nasa_range_1, $nasa_coeffs_1  ),\n' +
                        '       NASA(   $nasa_range_2, $nasa_coeffs_2  )\n' +
                        '               ),\n'
                        '    transport = gas_transport(\n' +
                        '                   geom = \"$transport_geometry\",\n' +
                        '                   diam = $diameter, \n' +
                        '                   well_depth = $well_depth, \n' +
                        '                   polar = $polar, \n' +
                        '                   rot_relax = $rot_relax, \n' +
                        '                   dipole= $dipole) \n' +
                        '        )\n\n')
        #write string template
        f.write(species_string.substitute(name=name, composition=composition, \
                    nasa_range_1=nasa_range_1, nasa_coeffs_1=nasa_coeffs_1,\
                    nasa_range_2=nasa_range_2, nasa_coeffs_2=nasa_coeffs_2,\
                    transport_geometry=transport_geometry, diameter=diameter,\
                    well_depth=well_depth, polar=polar, rot_relax=rot_relax,\
                    dipole=dipole))

    """-------------------------------------------------------------------------
    Write reactions to file
    -------------------------------------------------------------------------"""

    section_break('Reaction Data')

    #write data for each reaction in the Solution Object
    for n, i in enumerate(trimmed_solution.reaction_equations()):

        equation_string=trimmed_solution.reaction_equation(n)
        equation_object=trimmed_solution.reaction(n)
        equation_type=type(equation_object).__name__
        m=str(n+1)

        #Case if a ThreeBody Reaction
        if equation_type == 'ThreeBodyReaction':

            #trimms efficiencies list
            efficiencies=equation_object.efficiencies
            trimmed_efficiencies=equation_object.efficiencies
            for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]

            Arr=build_Arr(equation_object)
            replace_list_2={"{":"\"",      "\'":"",     ": ":":", \
                                    ",":" ",     "}":"\"" }
            Efficiencies_string=replace_multiple(str(trimmed_efficiencies), replace_list_2)

            reaction_string = Template('#  Reaction $m\n' +
                        'three_body_reaction( \"$equation_string\",  $Arr,\n' +
                        '       efficiencies = $Efficiencies) \n\n')

            f.write(reaction_string.substitute( m=m, \
                            equation_string=equation_string, Arr=Arr, \
                            Efficiencies=Efficiencies_string))

        #Case if an ElementaryReaction
        if equation_type == 'ElementaryReaction':
            Arr=build_Arr(equation_object)

            reaction_string=Template('#  Reaction $m\n'+
                                'reaction( \"$equation_string\", $Arr)\n\n')
            f.write(reaction_string.substitute(m=m, \
                            equation_string=equation_string, Arr=Arr))

        #Case if a FalloffReaction
        if equation_type == 'FalloffReaction':
            #trimms efficiencies list
            efficiencies=equation_object.efficiencies
            trimmed_efficiencies=equation_object.efficiencies
            for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]

            kf=build_mod_Arr(equation_object, 'high')
            kf0=build_mod_Arr(equation_object, 'low')
            replace_list_2={"{":"\"",      "\'":"",     ": ":":", \
                                    ",":" ",     "}":"\"" }
            Efficiencies_string=replace_multiple(str(trimmed_efficiencies), replace_list_2)

            reaction_string = Template('#  Reaction $m\n' +
                        'falloff_reaction( \"$equation_string\",\n' +
                        '        kf = $kf,\n' +
                        '        kf0   = $kf0,\n' +
                        '        efficiencies = $Efficiencies')
            f.write(reaction_string.substitute( m=m, \
                            equation_string=equation_string, kf=kf, kf0=kf0,\
                            Efficiencies=Efficiencies_string))
            j=equation_object.falloff.parameters
            #If optional Arrhenius data included:
            try:
                falloff_str=build_falloff(j)
                f.write(falloff_str)
            except (IndexError):
                f.write('\n           )\n\n')
                pass

    f.close()
    cw='atom ' + output_file_name
    os.system(cw)
