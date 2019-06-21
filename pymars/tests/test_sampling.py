""" Tests the sampling module in pyMARS """

import os
import pkg_resources

import pytest
import numpy as np
import cantera as ct
import yaml

from ..sampling import SamplingInputs, check_inputs

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
        inputs = SamplingInputs(
            input_ignition=relative_location(os.path.join('inputfiles', 'example_input_file.yaml')),
            )
        
        assert check_inputs(inputs)
    
    def test_good_example_alternate(self):
        """Tests correct sample input file with alternate values.
        """
        inputs = SamplingInputs(
            input_ignition=relative_location(os.path.join('inputfiles', 'example_input_file2.yaml')),
            )
        
        assert check_inputs(inputs)
    
    @pytest.mark.parametrize('key', [
        'kind', 'pressure', 'temperature', 'fuel', 'oxidizer', 'equivalence-ratio'
        ])
    def test_missing_keys(self, key):
        """Tests correct errors for missing required keys.
        """
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'end-time': 10.0,
            'fuel': {'CH4': 1.0},
            'oxidizer': {'O2': 1.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        del case[0][key]

        with TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, 'file.yaml')
            with open(filename, 'w') as the_file:
                yaml.dump(case, the_file)
            
            inputs = SamplingInputs(input_ignition=filename)

            with pytest.raises(KeyError):
                check_inputs(inputs)
    
    def test_bad_fuel_oxidizer_value(self):
        """Tests correct errors for improper value.
        """
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'end-time': 10.0,
            'fuel': {'CH4': 0.0},
            'oxidizer': {'O2': 1.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        with TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, 'file.yaml')
            with open(filename, 'w') as the_file:
                yaml.dump(case, the_file)
            
            inputs = SamplingInputs(input_ignition=filename)

            with pytest.raises(ValueError):
                check_inputs(inputs)
        
        case = [{
            'kind': 'constant volume',
            'pressure': 1.0,
            'temperature': 1000,
            'end-time': 10.0,
            'fuel': {'CH4': 1.0},
            'oxidizer': {'O2': 0.0, 'N2': 3.76},
            'equivalence-ratio': 1.0
            }]
        with TemporaryDirectory() as temp_dir:
            filename = os.path.join(temp_dir, 'file.yaml')
            with open(filename, 'w') as the_file:
                yaml.dump(case, the_file)
            
            inputs = SamplingInputs(input_ignition=filename)

            with pytest.raises(ValueError):
                check_inputs(inputs)
