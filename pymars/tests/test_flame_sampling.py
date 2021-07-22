""" Tests the sampling module in pyMARS """

import os
import pkg_resources

import pytest
import numpy as np
import cantera as ct
import yaml

from pymars.sampling import parse_flame_inputs, InputLaminarFlame

# Taken from http://stackoverflow.com/a/22726782/1569494
try:
    from tempfile import TemporaryDirectory
except ImportError:
    from contextlib import contextmanager
    import shutil
    import tempfile
    import errno

    @contextmanager
    def TemporaryDirectory():
        name = tempfile.mkdtemp()
        try:
            yield name
        finally:
            try:
                shutil.rmtree(name)
            except OSError as e:
                # Reraise unless ENOENT: No such file or directory
                # (ok if directory has already been deleted)
                if e.errno != errno.ENOENT:
                    raise

def relative_location(file):
    file_path = os.path.join(file)
    return pkg_resources.resource_filename(__name__, file_path)


class TestCheckInputs:
    def test_good_example(self):
        """Tests (good) sample input file.
        """
        inputs = [
            {'kind': 'constant volume', 'pressure': 1.0, 'temperature': 1000.0,
             'fuel': {'CH4': 1.0}, 'oxidizer': {'O2': 1.0, 'N2': 3.76}, 'equivalence-ratio': 1.0, 'width': 0.014},
            {'kind': 'constant volume', 'pressure': 1.0, 'temperature': 1200.0,
             'fuel': {'CH4': 1.0}, 'oxidizer': {'O2': 1.0, 'N2': 3.76}, 'equivalence-ratio': 1.0}
             ]
        
        conditions = parse_flame_inputs('gri30.cti', inputs)
        for item in conditions:
            assert type(item) == InputLaminarFlame
    
    def test_good_example_alternate(self):
        """Tests correct sample input file with alternate values.
        """
        inputs = [
            {'kind': 'constant pressure', 'pressure': 1.0, 'temperature': 1000.0,
             'width': 0.1, 'reactants': {'CH4': 1.0, 'O2': 1.0, 'N2': 3.76}},
             ]
        
        conditions = parse_flame_inputs('gri30.cti', inputs)
        for item in conditions:
            assert type(item) == InputLaminarFlame
    
    @pytest.mark.parametrize('key', [
        'kind', 'pressure', 'temperature', 'fuel', 'oxidizer', 'equivalence-ratio'
        ])
    def test_missing_keys(self, key):
        """Tests correct errors for missing required keys.
        """
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000.0,
            'fuel': {'CH4': 1.0},
            'oxidizer': {'O2': 1.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        del case[0][key]

        with pytest.raises(AssertionError):
            parse_flame_inputs('gri30.cti', case)
    
    def test_bad_fuel_oxidizer_value(self):
        """Tests correct errors for improper value.
        """
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'width': 0.1,
            'fuel': {'CH4': 0.0},
            'oxidizer': {'O2': 1.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        with pytest.raises(AssertionError):
            parse_flame_inputs('gri30.cti', case)
        
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'width': 0.1,
            'fuel': {'CH4': 1.0},
            'oxidizer': {'O2': 0.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        with pytest.raises(AssertionError):
            parse_flame_inputs('gri30.cti', case)
    
    def test_bad_species(self):
        """Tests raising error for species not in model.
        """
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'width': 0.1,
            'fuel': {'C4H10': 0.0},
            'oxidizer': {'O2': 1.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        with pytest.raises(AssertionError):
            parse_flame_inputs('gri30.cti', case)
        
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'width': 0.1,
            'fuel': {'CH4': 1.0},
            'oxidizer': {'O2': 0.0, 'HE': 3.76},
            'equivalence-ratio': 1.0
            }]
        with pytest.raises(AssertionError):
            parse_flame_inputs('gri30.cti', case)

        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'width': 0.1,
            'reactants': {'C4H10': 1.0, 'O2': 1.0, 'N2': 3.76}
            }]
        with pytest.raises(AssertionError):
            parse_flame_inputs('gri30.cti', case)

    def test_bad_width(self):
        """Tests raising error for species not in model.
        """
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'width': -2,
            'fuel': {'C4H10': 0.0},
            'oxidizer': {'O2': 1.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        with pytest.raises(AssertionError):
            parse_flame_inputs('gri30.cti', case)