from identify_file_extension import readin

A=readin('gri30.cti', ['H2, CO2, H4'])
initial=A[0]
final=A[1]

import cPickle
import os
os.system('rm -r testfile.cti')
f=open('testfile.cti', 'w+')

f.write('#'+ "-"*80 + '\n')
f.write('#  Reaction Data\n')
f.write('#'+ "-"*80 + '\n\n')

for n, i in enumerate(final.reaction_equations()):
    equation=final.reaction_equation(n)
    m=str(n+1)
    f.write('#  ' + 'Reaction' + ' ' + m + '\n')
    f.write('reaction(  ' + '\"'+ equation + '\",' + ')'+ '\n\n')

os.system('atom testfile.cti')

























"""
import xml.etree.cElementTree as ET

ctml= ET.Element("ctml")
phase=ET.SubElement(ctml, "phase", dim="3", id="gri30")
elementarray=ET.SubElement(phase, "elementarray", datasrc="elements.xml").text =" O H C N Ar"

tree=ET.ElementTree(ctml)
tree.write('tree.xml')
"""









"""
x=final.X
tad=final.T


csv_file='yourfile.csv'
with open(csv_file, 'w') as outfile:
    writer= csv.writer(outfile)
    writer.writerow(['T (K)'] + final.species_names)
    writer.writerow([tad]+list(x[i]))
"""
