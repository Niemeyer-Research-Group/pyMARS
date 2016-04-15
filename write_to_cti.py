"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff and ThreeBody Reactions
based on Cantera version 2.3.0a2"""

from identify_file_extension import readin

A=readin('gri30.cti', [])
initial=A[0]
final=A[1]

import os
os.system('rm -r testfile.cti')
f=open('testfile.cti', 'w+')


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
    f.write("species(name = \"" + name + "\"," + '\n')
    f.write("   atoms = "  + composition + '\n' )
    f.write("   thermo = (" + '\n')
    f.write(str(nasa_coeffs)+ '\n\n')








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
        Arr=["{:.2E}".format(equation_object.rate.pre_exponential_factor), equation_object.rate.temperature_exponent, equation_object.rate.activation_energy]
        Efficiencies=str(equation_object.efficiencies).replace("{", "\"").replace("\'", "").replace(": ", ":").replace(",", " ").replace("}", "\"")
        f.write('#  ' + 'Reaction' + ' ' + m + '\n')
        f.write('three_body_reaction(  ' + '\"'+ equation_string + '\",   ' + str(Arr))
        f.write('\n          ' + 'efficiencies = '+ Efficiencies + ')' + '\n\n')

    if equation_type == 'ElementaryReaction':
        Arr=["{:.2E}".format(equation_object.rate.pre_exponential_factor), equation_object.rate.temperature_exponent, equation_object.rate.activation_energy]
        f.write('#  ' + 'Reaction' + ' ' + m + '\n')
        f.write('reaction( ' + '\"'+ equation_string + '\",   ' + str(Arr)  + ')'+ '\n\n')

    if equation_type == 'FalloffReaction':
        Efficiencies=str(equation_object.efficiencies).replace("{", "\"").replace("\'", "").replace(": ", ":").replace(",", " ").replace("}", "\"")
        kf=["{:.2E}".format(equation_object.high_rate.pre_exponential_factor), equation_object.high_rate.temperature_exponent, equation_object.high_rate.activation_energy]
        kf0=["{:.2E}".format(equation_object.low_rate.pre_exponential_factor), equation_object.low_rate.temperature_exponent, equation_object.low_rate.activation_energy]
        j=equation_object.falloff.parameters
        f.write('#  ' + 'Reaction' + ' ' + m + '\n')
        f.write('falloff_reaction(  ' + '\"'+ equation_string + '\",   ' +'\n'  )
        f.write('           kf = ' + str(kf).replace("\'", "")+ '\n')
        f.write('           kf0 = ' + str(kf0).replace("\'", ""))
        try:
            falloff_str=str('\n' + '           falloff = Troe(A = ' + str(j[0]) + ', T3 = ' + str(j[1]) + ', T1 = ' + str(j[2]) + ', T2 = ' + str(j[3]))
            f.write(falloff_str)
        except (IndexError):
            pass
        f.write('\n          ' + ' efficiencies = '+ Efficiencies + ')' + '\n\n')

os.system('atom testfile.cti')
