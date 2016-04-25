"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required"""

from identify_file_extension import readin
import os
import textwrap
from string import Template


def write(input_file):



    input_file=os.path.abspath('Input_Data_Files/' + input_file)
    output_file_name=os.path.abspath('Output_Data_Files/'+'gri30_converted.cti')
    os.system('rm -r ' +  output_file_name)
    f=open(output_file_name, 'w+')


    #Read In trimmed Solution to Cantera Object
    Cantera_Solutions=readin(input_file, [])
    trimmed_solution=Cantera_Solutions[1]

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
        output_string= textwrap.fill(input_string, width=51, \
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
    phase_string= Template('ideal_gas(name = \"gri30\", \n \
        elements = \"$elements\", \n \
        species =""" $species""", \n\
        reactions = \"all\", \n \
        initial_state = state(temperature = 300.0, \n \
                                pressure= OneAtm)        \n\
        )\n\n')


    f.write(phase_string.substitute(elements=element_names, \
                                    species=species_names))

    """-------------------------------------------------------------------------
    Write Species to file
    -------------------------------------------------------------------------"""

    section_break('Species data')

    for i, name in enumerate(trimmed_solution.species_names):
        species=trimmed_solution.species(i)
        name=trimmed_solution.species(i).name

        nasa_coeffs=trimmed_solution.species(i).thermo.coeffs
        replace_list= {'{':'\"',       '}':'\"',       '\'':'',
                    ':  ':':',      '.0':"",         ',':'',       ' ': '  '}

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
        composition=replace_multiple(str(species.composition), replace_list )
        nasa_range_1=str([ species.thermo.min_temp, nasa_coeffs[0] ])
        nasa_range_2=str([ nasa_coeffs[0], species.thermo.max_temp ])
        transport_geometry=species.transport.geometry
        diameter= str(species.transport.diameter*(10**10))
        well_depth = str(species.transport.well_depth)
        polar = str(species.transport.polarizability*10**30)
        rot_relax = str(species.transport.rotational_relaxation)


        species_string=Template('species(name = "$name",\n\
            atoms = $composition, \n\
            thermo= (\n\
            NASA(   $nasa_range_1, $nasa_coeffs_1  ),\n\
            NASA(   $nasa_range_2, $nasa_coeffs_2  ),\n\
            transport = gas_transport(\n\
                    geom = \"$transport_geometry\", \n\
                    diam = $diameter, \n\
                    well_depth = $well_depth, \n\
                    polar = $polar, \n\
                    rot_relax = $rot_relax \n\
                )\n')

        f.write(species_string.substitute(name=name, composition=composition, \
                    nasa_range_1=nasa_range_1, nasa_coeffs_1=nasa_coeffs_1,\
                    nasa_range_2=nasa_range_2, nasa_coeffs_2=nasa_coeffs_2,\
                    transport_geometry=transport_geometry, diameter=diameter,\
                    well_depth=well_depth, polar=polar, rot_relax=rot_relax))

    """-----------------------------------------------------------------------------
    Write reactions to file
    -----------------------------------------------------------------------------"""

    section_break('Reaction Data')

    for n, i in enumerate(trimmed_solution.reaction_equations()):
        equation_string=trimmed_solution.reaction_equation(n)
        equation_object=trimmed_solution.reaction(n)
        equation_type=type(equation_object).__name__
        m=str(n+1)
        if equation_type == 'ThreeBodyReaction':
            Arr=[str("{:.5E}".format(equation_object.rate.pre_exponential_factor)).replace("\'", ""), equation_object.rate.temperature_exponent, equation_object.rate.activation_energy]
            Efficiencies=str(equation_object.efficiencies).replace("{", "\"").replace("\'", "").replace(": ", ":").replace(",", " ").replace("}", "\"")
            f.write('#  ' + 'Reaction' + ' ' + m + '\n')
            f.write('three_body_reaction(  ' + '\"'+ equation_string + '\",   ' + str(Arr).replace("\'", "")  +',')
            f.write('\n          ' + 'efficiencies = '+ Efficiencies + ')' + '\n\n')

        if equation_type == 'ElementaryReaction':
            Arr=["{:.5E}".format(equation_object.rate.pre_exponential_factor), equation_object.rate.temperature_exponent, equation_object.rate.activation_energy]
            f.write('#  ' + 'Reaction' + ' ' + m + '\n')
            f.write('reaction( ' + '\"'+ equation_string + '\",   ' + str(Arr).replace("\'", "")  + '),'+ '\n\n')

        if equation_type == 'FalloffReaction':
            Efficiencies=str(equation_object.efficiencies).replace("{", "\"").replace("\'", "").replace(": ", ":").replace(",", " ").replace("}", "\"")
            kf=["{:.5E}".format(equation_object.high_rate.pre_exponential_factor), equation_object.high_rate.temperature_exponent, equation_object.high_rate.activation_energy]
            kf0=["{:.5E}".format(equation_object.low_rate.pre_exponential_factor), equation_object.low_rate.temperature_exponent, equation_object.low_rate.activation_energy]
            j=equation_object.falloff.parameters
            f.write('#  ' + 'Reaction' + ' ' + m + '\n')
            f.write('falloff_reaction(  ' + '\"'+ equation_string + '\",   ' +'\n'  )
            f.write('           kf = ' + str(kf).replace("\'", "")+ ',\n')
            f.write('           kf0 = ' + str(kf0).replace("\'", "") +',')
            try:
                falloff_str=str('\n' + '           falloff = Troe(A = ' + str(j[0]) + ', T3 = ' + str(j[1]) + ', T1 = ' + str(j[2]) + ', T2 = ' + str(j[3]) +'),')
                f.write(falloff_str)
            except (IndexError):
                pass
            f.write('\n          ' + ' efficiencies = '+ Efficiencies + ')' + '\n\n')

    f.close()
    cw='atom ' + output_file_name
    os.system(cw)
write('gri30.cti')
