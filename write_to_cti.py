"""writes a solution object to a cantera cti file.

currently only works for Elementary, Falloff and ThreeBody Reactions"""

from identify_file_extension import readin
import re

A=readin('gri30.cti', [])
initial=A[0]
final=A[1]

"['H2', 'CO2']"
import cPickle
import os
os.system('rm -r testfile.cti')
f=open('testfile.cti', 'w+')


'Intro Block for Reactions'
f.write('#'+ "-"*80 + '\n')
f.write('#  Reaction Data\n')
f.write('#'+ "-"*80 + '\n\n')

for n, i in enumerate(final.reaction_equations()):
    equation_string=final.reaction_equation(n)
    equation_object=final.reaction(n)
    equation_type=type(equation_object).__name__
    if equation_type == 'ThreeBodyReaction':
        Arr=["{:.2E}".format(equation_object.rate.pre_exponential_factor), equation_object.rate.temperature_exponent, equation_object.rate.activation_energy]
        m=str(n+1)
        Efficiencies=str(equation_object.efficiencies).replace("{", "\"").replace("\'", "").replace(": ", ":").replace(",", " ").replace("}", "\"")
        f.write('#  ' + 'Reaction' + ' ' + m + '\n')
        f.write('three_body_reaction(  ' + '\"'+ equation_string + '\",   ' + str(Arr))
        f.write('\n          ' + 'efficiencies = '+ Efficiencies + ')' + '\n\n')

    if equation_type == 'ElementaryReaction':
        Arr=["{:.2E}".format(equation_object.rate.pre_exponential_factor), equation_object.rate.temperature_exponent, equation_object.rate.activation_energy]
        m=str(n+1)
        f.write('#  ' + 'Reaction' + ' ' + m + '\n')
        f.write('reaction( ' + '\"'+ equation_string + '\",   ' + str(Arr)  + ')'+ '\n\n')

    if equation_type == 'FalloffReaction':
        m=str(n+1)
        Efficiencies=str(equation_object.efficiencies).replace("{", "\"").replace("\'", "").replace(": ", ":").replace(",", " ")
        f.write('#  ' + 'Reaction' + ' ' + m + '\n')
        f.write('falloff_reaction(  ' + '\"'+ equation_string + '\",   '   + ')')
        f.write('\n          ' + 'efficiencies = '+ Efficiencies + ')' + '\n\n')

os.system('atom testfile.cti')
