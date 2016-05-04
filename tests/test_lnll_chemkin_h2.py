"""Tests LNLL h2 mechanism files
Complete Cantera mechanism file is left in output folder"""


"""
file:
    convert_chemkin_file.py

test parameters:
    Ck2cti works, and produces a mechanism file that can be used
    by Cantera to create a Soltuion
"""


import cantera as ct
from convert_chemkin_file import convert
import os

def test_convert_is_solution():
    file_name=convert('h2_v1b_mech.dat', 'h2_v1a_therm.dat', 'h2_v1a_tran.dat')
    solution=ct.Solution(file_name)
    assert solution.__qualname__ == 'Solution'
