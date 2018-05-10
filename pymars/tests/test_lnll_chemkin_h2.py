"""Tests LNLL h2 mechanism files
Complete Cantera mechanism file is left in output folder"""


"""
file:
    convert_chemkin_file.py

test parameters:
    Ck2cti works, and produces a mechanism file that can be used
    by Cantera to create a Soltuion
"""

import os
import cantera as ct
from pyMARS import convert_chemkin_file



def test_is_solution():
    #setup
    os.chdir('Data_Files')
    h2_package=['h2_v1b_mech.dat', 'h2_v1a_therm.dat', 'h2_v1a_tran.dat'  ]
    for i, files in enumerate(h2_package):
        os.system('cp ' + h2_package[i] +' ../Input_Data_Files')
    os.chdir('../')

    #test cti file
    file_name=convert_chemkin_file.convert('h2_v1b_mech.dat', h2_package[1], h2_package[2])
    solution=ct.Solution(file_name)
    assert solution.__qualname__ == 'Solution'

    #teardown
    os.chdir('Input_Data_Files')
    for i, files in enumerate(h2_package):
            os.system('rm ' + h2_package[i])
    return file_name

def test_create_trimmed_model(file_name):
    trimmed_solution_objects=os.system('python __main__.py --file=' + file_name)
