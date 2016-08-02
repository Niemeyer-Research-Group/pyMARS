"""writes a solution object to a chemkin inp file

currently only works for Elementary, Falloff and ThreeBody Reactions
Cantera development version 2.3.0a2 required"""

#from identify_file_extension import readin
import os
import textwrap
from string import Template
import cantera as ct
from tests.test_mechanism_from_solution import test
import ck2cti


def write(solution):
    """Function to write cantera solution object to inp file.

    Parameters
    ----------
    cantera Solution

    Returns
    -------
        Trimmed Mechanism file (.inp)
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
                    A=str("{:.3E}".format(equation_object.rate.pre_exponential_factor))
                if coeff_sum == 2:
                    A=str("{:.3E}".format(equation_object.rate.pre_exponential_factor*10**3))
                if coeff_sum == 3:
                    A=str("{:.3E}".format(equation_object.rate.pre_exponential_factor*10**6))
                #if equation_object.duplicate is True:

            if equation_type =='ThreeBodyReaction':
                if coeff_sum == 1:
                    A=str("{:.3E}".format(equation_object.rate.pre_exponential_factor*10**3))
                if coeff_sum == 2:
                    A=str("{:.3E}".format(equation_object.rate.pre_exponential_factor*10**6))

            if equation_type !='ElementaryReaction' and equation_type != 'ThreeBodyReaction':
                A=str("{:.3E}".format(equation_object.rate.pre_exponential_factor)) #*10**6
            b='{:.3f}'.format(equation_object.rate.temperature_exponent)
            E='{:.2f}'.format(equation_object.rate.activation_energy/c)
            Arr=[ A, b, E]
            return Arr

    def build_mod_Arr(equation_object, t_range):
        if t_range =='high':
            A=str("{:.3E}".format(equation_object.high_rate.pre_exponential_factor))  #*10**3
            b='{:.3f}'.format(equation_object.high_rate.temperature_exponent)
            E='{:.2f}'.format(equation_object.high_rate.activation_energy/c)
            Arr_high=[ A, b, E]
            return Arr_high
        if t_range == 'low':
            A=str("{:.3E}".format(equation_object.low_rate.pre_exponential_factor))   #*10**6
            b='{:.3f}'.format(equation_object.low_rate.temperature_exponent)
            E='{:.2f}'.format(equation_object.low_rate.activation_energy/c)
            Arr_low=[ A, b, E]
            return Arr_low

    def build_falloff(j):
        falloff_str=str(',\n        falloff = Troe(' +
                        'A = ' + str(j[0]) +
                        ', T3 = ' + str(j[1]) +
                        ', T1 = ' + str(j[2]) +
                        ', T2 = ' + str(j[3]) +')       )\n\n')
        return falloff_str

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
    f.write('THERMO ALL' +'\n')
    f.write('   300.000  1000.000  5000.000' +'\n')

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
        t_mid='{0:.3f}'.format(species.thermo.coeffs[0])


        temp_range= str(t_low) + '  ' + str(t_max) + '  ' + t_mid
        species_comp=''
        for ind, atom in enumerate(species.composition):
            species_comp += '{:<4}'.format(atom.upper())
            species_comp += str(int(species.composition[atom]))
        #species_comp += '   00'*(4-n_molecules)

        if type(species.transport).__name__ == 'GasTransportData':
            species_phase= 'G'
        else:
            print 'Species phase not found. Assumed to be Gas'
            species_phase='G'



        line_1 = '{:<18}'.format(name) + \
                    '{:<6}'.format('    ') + \
                '{:<20}'.format(species_comp) + \
                '{:<4}'.format(species_phase) + \
                '{:<31}'.format(temp_range) + \
                '{:<1}'.format('1') + \
                            '\n'
        f.write(line_1)

        line_2_coeffs=build_nasa(nasa_coeffs, 2)
        line_2 = line_2_coeffs  + '    2\n'
        f.write(line_2)

        line_3_coeffs=build_nasa(nasa_coeffs, 3)
        line_3 = line_3_coeffs + '    3\n'
        f.write(line_3)

        line_4_coeffs=build_nasa(nasa_coeffs, 4)
        line_4 = line_4_coeffs + '                   4\n'
        f.write(line_4)

    f.write('END\n')

    """-------------------------------------------------------------------------
    Write reactions to file
    -------------------------------------------------------------------------"""
    section_break('Reaction Data')
    f.write('REACTIONS\n')
    #write data for each reaction in the Solution Object
    for n, i in enumerate(trimmed_solution.reaction_equations()):
        equation_string=str(trimmed_solution.reaction_equation(n)).upper()
        equation_string=eliminate(equation_string, ' ', 'single')
        equation_object=trimmed_solution.reaction(n)
        equation_type=type(equation_object).__name__
        m=str(n+1)

        #Case if a ThreeBody Reaction
        if equation_type == 'ThreeBodyReaction':
            arrhenius=build_Arr(equation_object, equation_type)
            main_line= '{:<41}'.format(equation_string) + \
                        '{:>9}'.format(arrhenius[0])+\
                        '{:>9}'.format(arrhenius[1])+\
                        '{:>11}'.format(arrhenius[2])+\
                        '\n'
            f.write(main_line)

            #trimms efficiencies list
            efficiencies=equation_object.efficiencies
            trimmed_efficiencies=equation_object.efficiencies
            for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]

            replace_list_2={'{':'', '}':'/', '\'':'', ':':'/', ',':'/'}
            Efficiencies_string=replace_multiple(str(trimmed_efficiencies).upper(),\
            replace_list_2)

            secondary_line= str(Efficiencies_string) + '\n'
            if bool(efficiencies) is True:
                f.write(secondary_line)

        #Case if an elementary Reaction
        if equation_type == 'ElementaryReaction':
            arrhenius=build_Arr(equation_object, equation_type)
            main_line= '{:<41}'.format(equation_string) + \
                        '{:>9}'.format(arrhenius[0])+\
                        '{:>9}'.format(arrhenius[1])+\
                        '{:>11}'.format(arrhenius[2])+\
                        '\n'
            f.write(main_line)

        #Case if a FalloffReaction
        if equation_type == 'FalloffReaction':
            arr_high=build_mod_Arr(equation_object, 'high')
            main_line= '{:<41}'.format(equation_string) + \
                        '{:>9}'.format(arr_high[0])+\
                        '{:>9}'.format(arr_high[1])+\
                        '{:>11}'.format(arr_high[2])+\
                        '\n'
            f.write(main_line)

            arr_low=build_mod_Arr(equation_object, 'low')
            second_line = '     LOW  /' + \
                            '  ' + arr_low[0] +\
                            '  ' + arr_low[1] +\
                            '  ' + arr_low[2] + '/'+ '\n'
            f.write(second_line)
            j=equation_object.falloff.parameters
            #If optional Arrhenius data included:
            try:
                third_line = '     TROE/' + \
                                '   ' + str(j[0]) +\
                                '  ' + str(j[1]) +\
                                '  ' + str(j[2]) +\
                                '  ' + str(j[3]) +' /' + '\n'
                f.write(third_line)
            except (IndexError):
                pass

            #trimms efficiencies list
            efficiencies=equation_object.efficiencies
            trimmed_efficiencies=equation_object.efficiencies
            for s in efficiencies:
                    if s not in trimmed_solution.species_names:
                        del trimmed_efficiencies[s]

            replace_list_2={'{':'', '}':'/', '\'':'', ':':'/', ',':'/'}
            Efficiencies_string=replace_multiple(str(trimmed_efficiencies).upper(),\
            replace_list_2)

            fourth_line= str(Efficiencies_string) + '\n'
            if bool(efficiencies) is True:
                f.write(fourth_line)

        #dupluicate option
        if equation_object.duplicate is True:
            duplicate_line=' DUPLICATE' +'\n'
            f.write(duplicate_line)
    f.write('END')
    f.close()
    """-------------------------------------------------------------------------
    Test mechanism file
    -------------------------------------------------------------------------"""
    original_solution=solution

    parser = ck2cti.Parser()
    outName='test_file.cti'
    parser.convertMech(output_file_name, outName=outName)
    new_solution=ct.Solution(outName)
    test(original_solution, new_solution)
    os.remove(outName)
    return output_file_name



#A=ct.Solution('h2air_highT.cti')
#write(A)
