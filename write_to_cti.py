"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff and ThreeBody Reactions
based on Cantera version 2.3.0a2"""

from identify_file_extension import readin



import os
import textwrap

filename= 'gri30.cti'
A=readin(filename, [])
initial=A[0]
final=A[1]


file_path= os.path.relpath('Output_Data_Files/test_file.cti')
os.system('rm -r test_file.cti')
f=open('test_file.cti', 'w+')


"""-----------------------------------------------------------------------------
Write Title Block to file
-----------------------------------------------------------------------------"""
f.write('#'+ "-"*80 + '\n')
f.write('#  CTI File converted from Solution Object\n')
f.write('#'+ "-"*80 + '\n\n')

units_string="units(length = \"cm\", time = \"s\", quantity = \"mol\", act_energy = \"cal/mol\")"
f.write(units_string + '\n\n')

"""-----------------------------------------------------------------------------
Write Phase definition to file
-----------------------------------------------------------------------------"""





"""-----------------------------------------------------------------------------
Write Species to file
-----------------------------------------------------------------------------"""

f.write('#'+ "-"*80 + '\n')
f.write('#  Species Data\n')
f.write('#'+ "-"*80 + '\n\n')

for i, name in enumerate(final.species_names):
    species=final.species(i)
    name=final.species(i).name
    composition= str(species.composition).replace("{", "\"").replace("}", "\"").replace("\'", "").replace(": ", ":").replace(".0", "").replace(",", "").replace(" ", "  ")
    nasa_coeffs=species.thermo.coeffs
    nasa_range_1=str([species.thermo.min_temp, nasa_coeffs[0]])
    nasa_range_2=str([nasa_coeffs[0], species.thermo.max_temp])

    f.write("species(name = \"" + name + "\"," + '\n')
    f.write("   atoms = "  + composition + ', \n' )
    f.write("   thermo = (" + '\n')

    f.write('       NASA( ' + nasa_range_1 + ', [' + str(nasa_coeffs[1]) + ', ' + str(nasa_coeffs[2]) + ',\n')
    f.write('           ' + str(nasa_coeffs[3]) + ', ' + str(nasa_coeffs[4]) + ', ' + str(nasa_coeffs[5]) + ',' + '\n')
    f.write('           ' + str(nasa_coeffs[6]) + ', ' + str(nasa_coeffs[7]) + ']  ),\n')

    f.write('       NASA( ' + nasa_range_2 + ', [' + str(nasa_coeffs[8]) + ', ' + str(nasa_coeffs[9]) + ',\n')
    f.write('           ' + str(nasa_coeffs[10]) + ', ' + str(nasa_coeffs[11]) + ', ' + str(nasa_coeffs[12]) + ',' + '\n')
    f.write('           ' + str(nasa_coeffs[13]) + ', ' + str(nasa_coeffs[14]) + ']  )\n           ), \n  ')


    geom_string= '                      geom = '+ '\"' + species.transport.geometry + '\"' + ',\n'
    diam_string= '                      diam = ' + str(species.transport.diameter*(10**10)) + ',\n'
    well_depth_string='                      well_depth = ' + str(species.transport.well_depth*(10**22)) + ',\n'
    polar_string='                      polar = ' + str(species.transport.polarizability*10**30) + ',\n'
    rot_relax_string='                       rot_relax= ' + str(species.transport.rotational_relaxation)
    f.write('     transport = gas_transport( \n')
    f.write(geom_string)
    f.write(diam_string)
    f.write(well_depth_string)
    f.write(polar_string)
    f.write(rot_relax_string)
    f.write('         )        \n        )\n\n')


"""-----------------------------------------------------------------------------
Write reactions to file
-----------------------------------------------------------------------------"""

f.write('#'+ "-"*80 + '\n')
f.write('#  Reaction Data\n')
f.write('#'+ "-"*80 + '\n\n')

for n, i in enumerate(final.reaction_equations()):
    equation_string=final.reaction_equation(n)
    equation_object=final.reaction(n)
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

os.system('atom test_file.cti')
