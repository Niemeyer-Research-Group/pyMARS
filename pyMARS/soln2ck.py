"""writes a solution object to a chemkin inp file

currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required"""

#from identify_file_extension import readin
import os
import textwrap
from string import Template
import cantera as ct


def write(solution):
    """Function to write cantera solution object to inp file.

    Parameters
    ----------
    cantera Solution

    Returns
    -------
        Trimmed Mechanism file
    """
    trimmed_solution=solution
    input_file_name_stripped=trimmed_solution.name
    cwd= os.getcwd()
    output_file_name=os.path.join(cwd, 'pym_' + input_file_name_stripped + '.inp')

    f=open(output_file_name, 'w+')


    #Read In trimmed Solution to Cantera Object

    solution_T=trimmed_solution.T
    solution_P=trimmed_solution.P
    """-------------------------------------------------------------------------
    Work Functions
    -------------------------------------------------------------------------"""
    c=4184.0 #number of calories in 1000 Joules of energy
    def eliminate(input_string, char_to_replace, spaces='single'):
        for char in char_to_replace:
                    input_string= input_string.replace(char, "")
        if spaces == 'double':
                    input_string=input_string.replace(" ", "  ")
        return input_string

    def wrap(input_string):
        output_string= textwrap.fill(input_string, width=60) #subsequent_indent= '                        '
        return output_string
    def wrap_nasa(input_string):
        output_string= textwrap.fill(input_string, width=50, \
                                    subsequent_indent= '                ')
        return output_string

    def section_break(title):
        f.write('!'+ "-"*75 + '\n')
        f.write('!  ' + title +'\n')
        f.write('!'+ "-"*75 + '\n')

    def replace_multiple(input_string, replace_list):
        for a, b in replace_list.items():
            input_string= input_string.replace(a, b)
        return input_string

    def build_Arr(equation_object, equation_type):
            coeff_sum=sum(equation_object.reactants.values())
            if equation_type == 'ElementaryReaction':
                if coeff_sum == 1:
                    A=str("{:.5E}".format(equation_object.rate.pre_exponential_factor))
                if coeff_sum == 2:
                    A=str("{:.5E}".format(equation_object.rate.pre_exponential_factor*10**3))
                if coeff_sum == 3:
                    A=str("{:.5E}".format(equation_object.rate.pre_exponential_factor*10**6))
                #if equation_object.duplicate is True:

            if equation_type =='ThreeBodyReaction':
                if coeff_sum == 1:
                    A=str("{:.5E}".format(equation_object.rate.pre_exponential_factor*10**3))
                if coeff_sum == 2:
                    A=str("{:.5E}".format(equation_object.rate.pre_exponential_factor*10**6))

            if equation_type !='ElementaryReaction' and equation_type != 'ThreeBodyReaction':
                A=str("{:.5E}".format(equation_object.rate.pre_exponential_factor)) #*10**6
            b=equation_object.rate.temperature_exponent
            E=equation_object.rate.activation_energy/c
            Arr=[ A, b, E]
            return str(Arr).replace("\'", "")

    def build_mod_Arr(equation_object, t_range):
        if t_range =='high':
            A=str("{:.5E}".format(equation_object.high_rate.pre_exponential_factor))  #*10**3
            b=equation_object.high_rate.temperature_exponent
            E=equation_object.high_rate.activation_energy/c
            Arr_high=[ A, b, E]
            return str(Arr_high).replace("\'", "")
        if t_range == 'low':
            A=str("{:.5E}".format(equation_object.low_rate.pre_exponential_factor))   #*10**6
            b=equation_object.low_rate.temperature_exponent
            E=equation_object.low_rate.activation_energy/c
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
    section_break('Chemkin File converted from Solution Object by pyMARS')


    """-------------------------------------------------------------------------
    Write Phase definition to file
    -------------------------------------------------------------------------"""

    element_names=eliminate( str(trimmed_solution.element_names).upper(), \
                                        ['[', ']', '\'', ','])

    element_string=Template('ELEMENTS\n'    +
                            '$element_names\n' +
                            'END\n')
    f.write(element_string.substitute(element_names=element_names))

    species_names=wrap(
                        eliminate(str(trimmed_solution.species_names).upper(), \
                                    ['[', ']', '\'', ','], \
                                    spaces='double')
                        )
    species_string=Template('SPECIES\n' +
                    '$species_names\n'+
                    'END\n')

    f.write(species_string.substitute(species_names=species_names))


    """-------------------------------------------------------------------------
    Write Species to file
    -------------------------------------------------------------------------"""

    section_break('Species data')

    #write data for each species in the Solution object
    for i, name in enumerate(trimmed_solution.species_names):

        #physical Constant
        boltzmann=ct.boltzmann #joules/kelvin, boltzmann constant
        d=3.33564e-30 #1 debye = d coulomb-meters

        species=trimmed_solution.species(i)
        name=str(trimmed_solution.species(i).name).upper()
        nasa_coeffs=trimmed_solution.species(i).thermo.coeffs

        #Species attributes from trimmed solution object
        n_molecules = len(species.composition.keys())
        t_low='{0:.3f}'.format(species.thermo.min_temp)
        t_max='{0:.3f}'.format(species.thermo.max_temp)



        temp_range= str(t_low) + '  ' + str(t_max) + '  1000.000'
        species_comp=''
        for ind, atom in enumerate(species.composition):
            species_comp += atom.upper()
            species_comp += '   '
            species_comp += str(int(species.composition[atom]))
        #species_comp += '   00'*(4-n_molecules)

        if type(species.transport).__name__ == 'GasTransportData':
            species_phase= 'G'
        else:
            print 'Species phase not found. Assumed to be Gas'
            species_phase='G'



        line_1 = '{:<18}'.format(name) + \
                    '{:<6}'.format('date') + \
                '{:<20}'.format(species_comp) + \
                '{:<4}'.format(species_phase) + \
                '{:<31}'.format(temp_range) + \
                '{:<1}'.format('1') + \
                            '\n'
        f.write(line_1)


        def build_nasa(nasa_coeffs, row):
            line_coeffs=''
            lines=[[1,2,3,4,5], [6,7,8,9,10], [11,12,13,14]]
            line_index=lines[row-2]
            for ix, c in enumerate(nasa_coeffs):
                if ix in line_index:
                    if c >= 0:
                        line_coeffs += ' '
                    line_coeffs += str('{:.8e}'.format(c))
            return line_coeffs

        line_2_coeffs=build_nasa(nasa_coeffs, 2)
        line_2 = line_2_coeffs  + '    2\n'
        f.write(line_2)

        line_3_coeffs=build_nasa(nasa_coeffs, 3)
        line_3 = line_3_coeffs + '    3\n'
        f.write(line_3)

        line_4_coeffs=build_nasa(nasa_coeffs, 4)
        line_4 = line_4_coeffs + '                   4\n'
        f.write(line_4)


        """
        if species.transport.dipole != 0:
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
        if species.transport.dipole == 0:
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
                            '                   rot_relax = $rot_relax) \n' +
                            '        )\n\n')
            #write string template
            f.write(species_string.substitute(name=name, composition=composition, \
                        nasa_range_1=nasa_range_1, nasa_coeffs_1=nasa_coeffs_1,\
                        nasa_range_2=nasa_range_2, nasa_coeffs_2=nasa_coeffs_2,\
                        transport_geometry=transport_geometry, diameter=diameter,\
                        well_depth=well_depth, polar=polar, rot_relax=rot_relax,\
                        ))
            """


    """-------------------------------------------------------------------------
    Write reactions to file
    -------------------------------------------------------------------------"""

    section_break('Reaction Data')

    #write data for each reaction in the Solution Object
    for n, i in enumerate(trimmed_solution.reaction_equations()):
        equation_string=str(trimmed_solution.reaction_equation(n)).upper()
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

            Arr=build_Arr(equation_object, equation_type)
            replace_list_2={"{":"\"",      "\'":"",     ": ":":", \
                                    ",":" ",     "}":"\"" }
            Efficiencies_string=replace_multiple(str(trimmed_efficiencies).upper(),\
                        replace_list_2)

            reaction_string = Template('#  Reaction $m\n' +
                        'three_body_reaction( \"$equation_string\",  $Arr,\n' +
                        '       efficiencies = $Efficiencies) \n\n')

            f.write(reaction_string.substitute( m=m, \
                            equation_string=equation_string, Arr=Arr, \
                            Efficiencies=Efficiencies_string))
        #Case if an elementary Reaction
        if equation_type == 'ElementaryReaction':
            Arr=build_Arr(equation_object, equation_type)
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
            Efficiencies_string=replace_multiple(str(trimmed_efficiencies).upper(),\
                                        replace_list_2)

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
    return output_file_name
        #print ('trimmed mechanism file:  ' + output_file_name)
    f.close()

A=ct.Solution('gri301.cti')
write(A)
os.system('atom pym_gri30.inp')
