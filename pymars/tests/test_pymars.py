#tests to verify that various input species get included in retained or safe species

import os
import pkg_resources
import ruamel_yaml as yaml

from ..pymars import parse_inputs

def relative_location(file):
    file_path = os.path.join(file)
    return pkg_resources.resource_filename(__name__, file_path)


def test_species_check1():
    input_file = relative_location(os.path.join('assets', 'example1_input_file.yaml'))
    with open(input_file, 'r') as the_file:
        input_dict = yaml.safe_load(the_file)
    inputs = parse_inputs(input_dict)
    assert 'H2' in inputs.safe_species

def test_species_check2():
    input_file = relative_location(os.path.join('assets', 'example2_input_file.yaml'))
    with open(input_file, 'r') as the_file:
        input_dict = yaml.safe_load(the_file)
    inputs = parse_inputs(input_dict)
    assert 'AR' in inputs.safe_species

def test_species_check3():
    input_file = relative_location(os.path.join('assets', 'example3_input_file.yaml'))
    with open(input_file, 'r') as the_file:
        input_dict = yaml.safe_load(the_file)
    inputs = parse_inputs(input_dict)
    assert 'AR' in inputs.safe_species

def test_species_check4():
    input_file = relative_location(os.path.join('assets', 'example4_input_file.yaml'))
    with open(input_file, 'r') as the_file:
        input_dict = yaml.safe_load(the_file)
    inputs = parse_inputs(input_dict)
    assert 'AR' in inputs.safe_species

def test_species_check5():
    input_file = relative_location(os.path.join('assets', 'example5_input_file.yaml'))
    with open(input_file, 'r') as the_file:
        input_dict = yaml.safe_load(the_file)
    inputs = parse_inputs(input_dict)
    assert 'AR' in inputs.safe_species
